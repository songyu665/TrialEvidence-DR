import argparse
import gzip
import hashlib
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.control import RunController
from trialevidence.paths import project_path

BASE_URL = "https://clinicaltrials.gov/api/v2/studies"
ADVANCED_FILTER = "(AREA[Phase]PHASE2 OR AREA[Phase]PHASE3) AND AREA[HasResults]true"
FIELDS = [
    "NCTId", "OverallStatus", "WhyStopped", "Phase", "Condition", "InterventionType",
    "InterventionName", "ArmGroupType", "ArmGroupLabel", "ArmGroupDescription",
    "ArmGroupInterventionName", "OutcomeMeasureType", "OutcomeMeasureTitle",
    "OutcomeMeasureDescription", "OutcomeMeasurePopulationDescription", "OutcomeMeasureTimeFrame",
    "OutcomeMeasureReportingStatus", "OutcomeMeasureParamType", "OutcomeMeasureUnitOfMeasure",
    "OutcomeAnalysisParamType", "OutcomeAnalysisParamValue", "OutcomeAnalysisPValue",
    "OutcomeAnalysisCILowerLimit", "OutcomeAnalysisCIUpperLimit", "OutcomeAnalysisTestedNonInferiority",
    "OutcomeAnalysisNonInferiorityType", "OutcomeAnalysisNonInferiorityComment",
    "OutcomeAnalysisGroupId", "OutcomeAnalysisGroupDescription", "OutcomeAnalysisStatisticalMethod",
    "OutcomeAnalysisStatisticalComment", "OutcomeAnalysisPValueComment",
    "OutcomeAnalysisEstimateComment", "OutcomeAnalysisOtherAnalysisDescription", "ResultsFirstPostDate",
]


def request_page(page_token: str | None, page_size: int) -> bytes:
    parameters = {
        "format": "json", "filter.overallStatus": "COMPLETED", "filter.advanced": ADVANCED_FILTER,
        "pageSize": str(page_size), "fields": ",".join(FIELDS),
    }
    if page_token:
        parameters["pageToken"] = page_token
    else:
        parameters["countTotal"] = "true"
    url = f"{BASE_URL}?{urllib.parse.urlencode(parameters)}"
    request = urllib.request.Request(url, headers={"User-Agent": "TrialEvidence-DR-E0/0.1"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read()
        except urllib.error.HTTPError as error:
            if error.code not in {429, 500, 502, 503, 504} or attempt == 4:
                raise
            retry_after = error.headers.get("Retry-After")
            time.sleep(float(retry_after) if retry_after else 2**attempt)
        except (urllib.error.URLError, TimeoutError):
            if attempt == 4:
                raise
            time.sleep(2**attempt)
    raise RuntimeError("Unreachable retry state")


def atomic_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("wb") as stream:
        stream.write(payload)
        stream.flush()
    temporary.replace(path)


def atomic_json(path: Path, payload: dict) -> None:
    atomic_bytes(path, (json.dumps(payload, ensure_ascii=False, indent=2) + "\n").encode("utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-name", default="e0_trials_completed_results_phase23_with_comments")
    parser.add_argument("--page-size", type=int, default=1000)
    parser.add_argument("--max-pages", type=int)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    controller = RunController(args.run_name)
    controller.install_signal_handlers()
    output_dir = project_path("data", "raw", "clinicaltrials", "completed_results_phase23_with_comments")
    checkpoint = output_dir / "download_checkpoint.json"
    if args.resume:
        state = controller.load_checkpoint()
        controller.clear_requests()
    else:
        if checkpoint.exists():
            raise FileExistsError("Download checkpoint exists; use --resume")
        state = {
            "next_token": None, "pages_downloaded": 0, "studies_downloaded": 0,
            "total_count": None, "page_size": args.page_size,
            "status_filter": "COMPLETED", "advanced_filter": ADVANCED_FILTER, "fields": FIELDS,
        }
    if state["page_size"] != args.page_size:
        raise ValueError("Cannot change page size when resuming")

    complete = False
    while True:
        reason = controller.stop_reason
        if reason:
            controller.save_checkpoint(state)
            raise SystemExit(85)
        if args.max_pages is not None and state["pages_downloaded"] >= args.max_pages:
            break
        if state["pages_downloaded"] > 0 and not state["next_token"]:
            complete = True
            break
        page_number = state["pages_downloaded"] + 1
        page_path = output_dir / f"page_{page_number:06d}.json.gz"
        if page_path.exists():
            with gzip.open(page_path, "rb") as stream:
                payload = json.load(stream)
        else:
            body = request_page(state["next_token"], args.page_size)
            payload = json.loads(body)
            if not isinstance(payload.get("studies"), list):
                raise ValueError("Response lacks studies array")
            atomic_bytes(page_path, gzip.compress(body, mtime=0))
            atomic_json(
                output_dir / f"page_{page_number:06d}.manifest.json",
                {
                    "page": page_number, "studies": len(payload["studies"]),
                    "sha256_uncompressed_json": hashlib.sha256(body).hexdigest(),
                    "downloaded_at": datetime.now(timezone.utc).isoformat(),
                },
            )
        studies = payload["studies"]
        state["pages_downloaded"] = page_number
        state["studies_downloaded"] += len(studies)
        state["next_token"] = payload.get("nextPageToken")
        if state["total_count"] is None:
            state["total_count"] = payload.get("totalCount")
        controller.save_checkpoint(state)
        atomic_json(checkpoint, state)
        print(json.dumps({"page": page_number, "studies": len(studies), "next_page": bool(state["next_token"])}), flush=True)
        if not state["next_token"]:
            complete = True
            break
    summary = {
        "finished_at": datetime.now(timezone.utc).isoformat(),
        "complete_query_download": complete, "max_pages": args.max_pages, **state,
    }
    atomic_json(output_dir / "download_summary.json", summary)
    controller.mark_completed(summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

import argparse
import gzip
import hashlib
import json
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.control import RunController
from trialevidence.paths import project_path

BASE_URL = "https://clinicaltrials.gov/api/v2/studies"
ADVANCED_FILTER = "(AREA[Phase]PHASE2 OR AREA[Phase]PHASE3)"


def request_page(page_token: str | None, page_size: int) -> tuple[bytes, str]:
    parameters = {
        "format": "json",
        "filter.overallStatus": "TERMINATED",
        "filter.advanced": ADVANCED_FILTER,
        "pageSize": str(page_size),
    }
    if page_token:
        parameters["pageToken"] = page_token
    else:
        parameters["countTotal"] = "true"
    url = f"{BASE_URL}?{urllib.parse.urlencode(parameters)}"
    request = urllib.request.Request(url, headers={"User-Agent": "TrialEvidence-DR-E0/0.1"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                return response.read(), url
        except urllib.error.HTTPError as error:
            if error.code not in {429, 500, 502, 503, 504} or attempt == 4:
                raise
            retry_after = error.headers.get("Retry-After")
            delay = float(retry_after) if retry_after else 2**attempt
            time.sleep(delay)
        except (urllib.error.URLError, TimeoutError):
            if attempt == 4:
                raise
            time.sleep(2**attempt)
    raise RuntimeError("Unreachable retry state")


def atomic_bytes(path: Path, payload: bytes) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("wb") as stream:
        stream.write(payload)
        stream.flush()
    temporary.replace(path)
    return hashlib.sha256(payload).hexdigest()


def atomic_json(path: Path, payload: dict) -> None:
    atomic_bytes(path, (json.dumps(payload, ensure_ascii=False, indent=2) + "\n").encode("utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-name", default="e0_trials_terminated_phase23")
    parser.add_argument("--page-size", type=int, default=1000)
    parser.add_argument("--max-pages", type=int)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.page_size <= 1000:
        parser.error("page-size must be between 1 and 1000")

    controller = RunController(args.run_name)
    controller.install_signal_handlers()
    output_dir = project_path("data", "raw", "clinicaltrials", "terminated_phase23")
    checkpoint = output_dir / "download_checkpoint.json"

    if args.resume:
        state = controller.load_checkpoint()
        controller.clear_requests()
    else:
        if checkpoint.exists():
            raise FileExistsError("Download checkpoint exists; use --resume or choose another run name")
        state = {
            "next_token": None,
            "pages_downloaded": 0,
            "studies_downloaded": 0,
            "total_count": None,
            "page_size": args.page_size,
            "status_filter": "TERMINATED",
            "advanced_filter": ADVANCED_FILTER,
        }
    if state["page_size"] != args.page_size:
        raise ValueError("Cannot change page size when resuming")

    started = time.monotonic()
    complete = False
    while True:
        reason = controller.stop_reason
        if reason:
            controller.save_checkpoint(state)
            print(json.dumps({"paused": reason, "state": state}))
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
            url = None
        else:
            body, url = request_page(state["next_token"], args.page_size)
            payload = json.loads(body.decode("utf-8"))
            if not isinstance(payload.get("studies"), list):
                raise ValueError("ClinicalTrials.gov response lacks a studies array")
            digest = hashlib.sha256(body).hexdigest()
            compressed = gzip.compress(body, mtime=0)
            atomic_bytes(page_path, compressed)
            page_meta = {
                "page": page_number,
                "file": str(page_path.relative_to(project_path())),
                "sha256_uncompressed_json": digest,
                "studies": len(payload["studies"]),
                "request_url": url,
                "downloaded_at": datetime.now(timezone.utc).isoformat(),
            }
            atomic_json(output_dir / f"page_{page_number:06d}.manifest.json", page_meta)

        studies = payload["studies"]
        state["pages_downloaded"] = page_number
        state["studies_downloaded"] += len(studies)
        state["next_token"] = payload.get("nextPageToken")
        if state["total_count"] is None:
            state["total_count"] = payload.get("totalCount")
        controller.save_checkpoint(state)
        atomic_json(checkpoint, state)
        print(json.dumps({"page": page_number, "studies": len(studies), "next_page": bool(state["next_token"])}))
        if not state["next_token"]:
            complete = True
            break

    summary = {
        "finished_at": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": round(time.monotonic() - started, 3),
        "complete_query_download": complete,
        "max_pages": args.max_pages,
        **state,
    }
    atomic_json(output_dir / "download_summary.json", summary)
    controller.mark_completed(summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

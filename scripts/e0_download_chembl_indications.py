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

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/drug_indication.json?limit=1000"


def fetch(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "TrialEvidence-DR-E0/0.1", "Accept": "application/json"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
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
    parser.add_argument("--run-name", default="e0_chembl_drug_indications")
    parser.add_argument("--max-pages", type=int)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    controller = RunController(args.run_name)
    controller.install_signal_handlers()
    output_dir = project_path("data", "raw", "chembl", "drug_indication")
    checkpoint = output_dir / "download_checkpoint.json"

    if args.resume:
        state = controller.load_checkpoint()
        controller.clear_requests()
    else:
        if checkpoint.exists():
            raise FileExistsError("Download checkpoint exists; use --resume")
        state = {"next_url": BASE_URL, "pages_downloaded": 0, "records_downloaded": 0, "total_count": None}

    complete = False
    while True:
        reason = controller.stop_reason
        if reason:
            controller.save_checkpoint(state)
            print(json.dumps({"paused": reason, "state": state}))
            raise SystemExit(85)
        if not state["next_url"]:
            complete = True
            break
        if args.max_pages is not None and state["pages_downloaded"] >= args.max_pages:
            break

        page_number = state["pages_downloaded"] + 1
        request_url = state["next_url"]
        page_path = output_dir / f"page_{page_number:06d}.json.gz"
        if page_path.exists():
            with gzip.open(page_path, "rb") as stream:
                payload = json.load(stream)
        else:
            body = fetch(request_url)
            payload = json.loads(body)
            records = payload.get("drug_indications")
            if not isinstance(records, list):
                raise ValueError("ChEMBL response lacks drug_indications array")
            atomic_bytes(page_path, gzip.compress(body, mtime=0))
            atomic_json(
                output_dir / f"page_{page_number:06d}.manifest.json",
                {
                    "page": page_number,
                    "request_url": request_url,
                    "records": len(records),
                    "sha256_uncompressed_json": hashlib.sha256(body).hexdigest(),
                    "downloaded_at": datetime.now(timezone.utc).isoformat(),
                },
            )
        records = payload["drug_indications"]
        state["records_downloaded"] += len(records)
        state["pages_downloaded"] = page_number
        state["total_count"] = payload.get("page_meta", {}).get("total_count")
        next_path = payload.get("page_meta", {}).get("next")
        state["next_url"] = urllib.parse.urljoin("https://www.ebi.ac.uk", next_path) if next_path else None
        controller.save_checkpoint(state)
        atomic_json(checkpoint, state)
        print(json.dumps({"page": page_number, "records": len(records), "next_page": bool(state["next_url"])}), flush=True)

    summary = {
        "finished_at": datetime.now(timezone.utc).isoformat(),
        "complete_query_download": complete,
        **state,
    }
    atomic_json(output_dir / "download_summary.json", summary)
    controller.mark_completed(summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

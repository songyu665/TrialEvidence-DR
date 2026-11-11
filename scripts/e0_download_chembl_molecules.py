import argparse
import csv
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

API = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"


def atomic_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("wb") as stream:
        stream.write(payload)
        stream.flush()
    temporary.replace(path)


def atomic_json(path: Path, payload: dict) -> None:
    atomic_bytes(path, (json.dumps(payload, ensure_ascii=False, indent=2) + "\n").encode("utf-8"))


def fetch(ids: list[str]) -> bytes:
    url = API + "?" + urllib.parse.urlencode({"molecule_chembl_id__in": ",".join(ids), "limit": 1000})
    request = urllib.request.Request(url, headers={"User-Agent": "TrialEvidence-DR-E0/0.1", "Accept": "application/json"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read()
        except Exception:
            if attempt == 4:
                raise
            time.sleep(2**attempt)
    raise RuntimeError("Unreachable retry state")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-name", default="e0_chembl_molecule_dictionary")
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--max-batches", type=int)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    source = project_path("data", "interim", "chembl_approved_indication_candidates.csv")
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    with source.open(encoding="utf-8-sig", newline="") as stream:
        ids = sorted({row["parent_molecule_chembl_id"] for row in csv.DictReader(stream)})
    batches = [ids[index:index + args.batch_size] for index in range(0, len(ids), args.batch_size)]

    controller = RunController(args.run_name)
    controller.install_signal_handlers()
    output_dir = project_path("data", "raw", "chembl", "molecule_dictionary")
    checkpoint = output_dir / "download_checkpoint.json"
    if args.resume:
        state = controller.load_checkpoint()
        controller.clear_requests()
        if state["source_sha256"] != source_hash or state["batch_size"] != args.batch_size:
            raise ValueError("Candidate source or batch size changed; start a new run")
    else:
        if checkpoint.exists():
            raise FileExistsError("Molecule checkpoint exists; use --resume")
        state = {
            "source_sha256": source_hash,
            "batch_size": args.batch_size,
            "candidate_molecules": len(ids),
            "next_batch": 0,
            "records_downloaded": 0,
        }

    while state["next_batch"] < len(batches):
        reason = controller.stop_reason
        if reason:
            controller.save_checkpoint(state)
            raise SystemExit(85)
        if args.max_batches is not None and state["next_batch"] >= args.max_batches:
            break
        batch_index = state["next_batch"]
        requested_ids = batches[batch_index]
        body = fetch(requested_ids)
        payload = json.loads(body)
        records = payload.get("molecules", [])
        returned_ids = {record["molecule_chembl_id"] for record in records}
        if not returned_ids.issubset(requested_ids):
            raise ValueError("ChEMBL returned unrequested molecule IDs")
        number = batch_index + 1
        atomic_bytes(output_dir / f"batch_{number:06d}.json.gz", gzip.compress(body, mtime=0))
        atomic_json(
            output_dir / f"batch_{number:06d}.manifest.json",
            {
                "batch": number,
                "requested_ids": requested_ids,
                "returned_ids": sorted(returned_ids),
                "records": len(records),
                "sha256_uncompressed_json": hashlib.sha256(body).hexdigest(),
                "downloaded_at": datetime.now(timezone.utc).isoformat(),
            },
        )
        state["next_batch"] = number
        state["records_downloaded"] += len(records)
        controller.save_checkpoint(state)
        atomic_json(checkpoint, state)
        print(json.dumps({"batch": number, "batches_total": len(batches), "records": len(records)}), flush=True)

    complete = state["next_batch"] == len(batches)
    summary = {
        "finished_at": datetime.now(timezone.utc).isoformat(),
        "complete_query_download": complete,
        **state,
    }
    atomic_json(output_dir / "download_summary.json", summary)
    controller.mark_completed(summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

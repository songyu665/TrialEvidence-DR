import gzip
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path

DATASETS = {
    "clinicaltrials_terminated_phase23": {
        "directory": "data/raw/clinicaltrials/terminated_phase23",
        "array_key": "studies",
        "id_path": ("protocolSection", "identificationModule", "nctId"),
    },
    "clinicaltrials_completed_results_phase23": {
        "directory": "data/raw/clinicaltrials/completed_results_phase23",
        "array_key": "studies",
        "id_path": ("protocolSection", "identificationModule", "nctId"),
    },
    "clinicaltrials_completed_results_phase23_with_comments": {
        "directory": "data/raw/clinicaltrials/completed_results_phase23_with_comments",
        "array_key": "studies",
        "id_path": ("protocolSection", "identificationModule", "nctId"),
    },
    "chembl_drug_indication": {
        "directory": "data/raw/chembl/drug_indication",
        "array_key": "drug_indications",
        "id_path": ("drugind_id",),
    },
}


def nested(record: dict, path: tuple[str, ...]):
    value = record
    for key in path:
        value = value.get(key, {})
    return value if value != {} else None


def verify(name: str, config: dict) -> dict:
    directory = project_path(config["directory"])
    state = json.loads((directory / "download_checkpoint.json").read_text(encoding="utf-8"))
    count = 0
    ids = []
    for page in range(1, state["pages_downloaded"] + 1):
        data_path = directory / f"page_{page:06d}.json.gz"
        manifest_path = directory / f"page_{page:06d}.manifest.json"
        with gzip.open(data_path, "rb") as stream:
            body = stream.read()
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        digest = hashlib.sha256(body).hexdigest()
        if digest != manifest["sha256_uncompressed_json"]:
            raise ValueError(f"SHA-256 mismatch: {name} page {page}")
        payload = json.loads(body)
        records = payload[config["array_key"]]
        if len(records) != manifest["studies" if name.startswith("clinical") else "records"]:
            raise ValueError(f"Record count mismatch: {name} page {page}")
        count += len(records)
        ids.extend(nested(record, config["id_path"]) for record in records)
    expected_count = state["studies_downloaded" if name.startswith("clinical") else "records_downloaded"]
    if count != expected_count:
        raise ValueError(f"Total count mismatch for {name}")
    if any(identifier is None for identifier in ids):
        raise ValueError(f"Missing identifiers in {name}")
    duplicates = len(ids) - len(set(ids))
    if duplicates:
        raise ValueError(f"Duplicate identifiers in {name}: {duplicates}")
    return {
        "pages_verified": state["pages_downloaded"],
        "records_verified": count,
        "unique_ids_verified": len(set(ids)),
        "download_complete": not state.get("next_token") and not state.get("next_url"),
        "status": "passed",
    }


def main() -> None:
    results = {name: verify(name, config) for name, config in DATASETS.items()}
    report = {"verified_at": datetime.now(timezone.utc).isoformat(), "datasets": results}
    output = project_path("reports", "raw_download_integrity.json")
    with output.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

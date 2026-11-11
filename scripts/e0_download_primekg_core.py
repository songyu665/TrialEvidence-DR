import hashlib
import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path

DATASET_URL = "https://dataverse.harvard.edu/api/datasets/:persistentId?persistentId=doi:10.7910/DVN/IXA7BM"
ACCESS_URL = "https://dataverse.harvard.edu/api/access/datafile/"
WANTED = {"nodes.tab", "drug_features.tab", "kg_grouped_diseases.tab", "kg_grouped_diseases_bert_map.tab", "README.txt"}
USER_AGENT = "TrialEvidence-DR-E0/0.1"


def fetch(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read()


def atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("wb") as stream:
        stream.write(payload)
        stream.flush()
    temporary.replace(path)


def main() -> None:
    output_dir = project_path("data", "raw", "primekg")
    metadata_bytes = fetch(DATASET_URL)
    metadata = json.loads(metadata_bytes)
    atomic_write(output_dir / "dataset_metadata.json", metadata_bytes)
    version = metadata["data"]["latestVersion"]
    available = {item["dataFile"]["filename"]: item["dataFile"] for item in version["files"]}
    missing = WANTED - available.keys()
    if missing:
        raise KeyError(f"PrimeKG files missing from metadata: {sorted(missing)}")

    manifest = []
    for name in sorted(WANTED):
        item = available[name]
        payload = fetch(f"{ACCESS_URL}{item['id']}?format=original")
        digest = hashlib.md5(payload).hexdigest()
        expected = item.get("checksum", {}).get("value")
        if expected and digest != expected:
            raise ValueError(f"Checksum mismatch for {name}")
        atomic_write(output_dir / name, payload)
        manifest.append(
            {
                "name": name,
                "dataverse_file_id": item["id"],
                "bytes": len(payload),
                "metadata_tabular_size": item.get("filesize"),
                "md5": digest,
                "verified": True,
            }
        )
        print(json.dumps(manifest[-1]))

    summary = {
        "downloaded_at": datetime.now(timezone.utc).isoformat(),
        "dataset_persistent_url": metadata["data"].get("persistentUrl"),
        "dataset_version": f"{version.get('versionNumber')}.{version.get('versionMinorNumber')}",
        "dataset_release_time": version.get("releaseTime"),
        "files": manifest,
        "scope": "Node, drug-feature and source README files only; full graph edges not yet downloaded",
    }
    atomic_write(
        output_dir / "core_download_manifest.json",
        (json.dumps(summary, ensure_ascii=False, indent=2) + "\n").encode("utf-8"),
    )


if __name__ == "__main__":
    main()

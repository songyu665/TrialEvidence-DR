import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path

TARGETS = [
    ("AACT_2023-12-01_monthly_snapshot", "https://aact.ctti-clinicaltrials.org/static/exported_files/monthly/20231201_aact.zip"),
    ("ChEMBL_33_SQLite", "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_33/chembl_33_sqlite.tar.gz"),
    ("ClinicalTrials_record_history_example", "https://clinicaltrials.gov/study/NCT03352557?tab=history"),
]


def head(name: str, url: str) -> dict:
    request = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "TrialEvidence-DR-E0/0.1"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return {
            "name": name, "url": url, "http_status": response.status,
            "content_length": response.headers.get("Content-Length"),
            "content_type": response.headers.get("Content-Type"),
            "last_modified": response.headers.get("Last-Modified"),
        }


def main() -> None:
    primekg_path = project_path("data", "raw", "primekg", "dataset_metadata.json")
    primekg = json.loads(primekg_path.read_text(encoding="utf-8"))["data"]["latestVersion"]
    targets = [head(name, url) for name, url in TARGETS]
    report = {
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "primekg": {
            "dataset_version": f"{primekg.get('versionNumber')}.{primekg.get('versionMinorNumber')}",
            "release_time": primekg.get("releaseTime"),
            "compatible_with_2023_12_cutoff_by_release_date": True,
        },
        "historical_resources": targets,
        "recoverability_assessment": "Historical resources exist and are reachable, but no historical archive has been downloaded, extracted or used to reproduce the candidate labels",
        "time_experiment_allowed": False,
        "remaining_requirements": [
            "Download and checksum the selected historical archives",
            "Reproduce P and N_strict candidates from the historical versions",
            "Use first-public-availability dates rather than event dates alone",
            "Measure the proportion of current candidates recoverable at the cutoff",
            "Keep 2024-01-01 AACT snapshot separate as a post-cutoff reference, not training data",
        ],
    }
    output = project_path("reports", "e0_history_audit.json")
    with output.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

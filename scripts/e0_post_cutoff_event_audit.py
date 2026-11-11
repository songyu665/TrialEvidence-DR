import csv
import gzip
import json
import sqlite3
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path


def read_csv(name: str) -> list[dict]:
    with project_path("data", "interim", name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    historical = {
        (row["drug_node_index"], row["disease_node_index"])
        for row in read_csv("historical_chembl33_positive_pairs.csv")
    }
    current_rows = [row for row in read_csv("chembl_primekg_positive_candidates.csv") if row["mapping_status"] == "both_mapped_unique"]
    current = {}
    for row in current_rows:
        key = (row["drug_node_index"], row["disease_node_index"])
        if key not in current:
            current[key] = {
                "drug_node_index": key[0], "disease_node_index": key[1],
                "drug_node_name": row["drug_node_name"], "disease_node_name": row["disease_node_name"],
                "parent_molecule_chembl_ids": set(), "drugind_ids": set(),
                "reference_types": set(), "reference_urls": set(),
            }
        current[key]["parent_molecule_chembl_ids"].add(row["parent_molecule_chembl_id"])
        current[key]["drugind_ids"].update(json.loads(row["drugind_ids"]))
        current[key]["reference_types"].update(json.loads(row["reference_types"]))
        current[key]["reference_urls"].update(json.loads(row["reference_urls"]))

    molecule_dir = project_path("data", "raw", "chembl", "molecule_dictionary")
    summary = json.loads((molecule_dir / "download_summary.json").read_text())
    molecules = {}
    for batch in range(1, summary["next_batch"] + 1):
        with gzip.open(molecule_dir / f"batch_{batch:06d}.json.gz", "rt", encoding="utf-8") as stream:
            for record in json.load(stream)["molecules"]:
                molecules[record["molecule_chembl_id"]] = record

    database = project_path("data", "raw", "chembl", "historical", "chembl_33.db")
    con = sqlite3.connect(f"file:{database}?mode=ro", uri=True)
    historical_drugind_ids = {str(row[0]) for row in con.execute("SELECT drugind_id FROM drug_indication")}
    con.close()

    output_rows = []
    for key, row in current.items():
        if key in historical:
            continue
        current_ids = row["drugind_ids"]
        old_ids = current_ids & historical_drugind_ids
        new_ids = current_ids - historical_drugind_ids
        if old_ids and new_ids:
            category = "mixed_existing_and_new_indication_records"
        elif new_ids:
            category = "new_indication_records_after_2023_version"
        else:
            category = "existing_records_with_mapping_or_curation_difference"
        approvals = sorted({molecules[identifier].get("first_approval") for identifier in row["parent_molecule_chembl_ids"] if identifier in molecules and molecules[identifier].get("first_approval")})
        output_rows.append({
            "drug_node_index": key[0], "disease_node_index": key[1],
            "drug_node_name": row["drug_node_name"], "disease_node_name": row["disease_node_name"],
            "parent_molecule_chembl_ids": json.dumps(sorted(row["parent_molecule_chembl_ids"])),
            "first_approval_years": json.dumps(approvals),
            "current_drugind_ids": json.dumps(sorted(current_ids)),
            "new_drugind_ids_relative_to_chembl33": json.dumps(sorted(new_ids)),
            "existing_drugind_ids_in_chembl33": json.dumps(sorted(old_ids)),
            "reference_types": json.dumps(sorted(row["reference_types"])),
            "reference_urls": json.dumps(sorted(row["reference_urls"]), ensure_ascii=False),
            "event_category": category,
            "event_date_confirmation_status": "pending_original_regulatory_or_label_date_review",
        })
    output = project_path("data", "interim", "post_cutoff_positive_event_candidates.csv")
    with output.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "current_positive_pairs_not_in_historical_set": len(output_rows),
        "event_category_counts": dict(Counter(row["event_category"] for row in output_rows)),
        "with_drug_first_approval_after_2023": sum(any(year > 2023 for year in json.loads(row["first_approval_years"])) for row in output_rows),
        "event_date_confirmation_completed": False,
        "time_experiment_allowed": False,
        "warning": "A new ChEMBL record ID is not proof of a new approval or its public date; first_approval refers to the drug, not the indication",
    }
    report_path = project_path("reports", "e0_post_cutoff_event_audit.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

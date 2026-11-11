import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path


def read_csv(name: str) -> list[dict]:
    with project_path("data", "interim", name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    reviewed = read_csv("n_strict_conservative_rule_review.csv")
    passed = [row for row in reviewed if row["rule_review_status"] == "passed_conservative_rules"]
    grouped = {}
    for row in passed:
        key = (row["drug_node_index"], row["disease_node_index"])
        if key not in grouped:
            grouped[key] = {
                "drug_node_index": key[0], "disease_node_index": key[1],
                "drug_node_name": row["drug_node_name"], "disease_node_name": row["disease_node_name"],
                "trial_ids": set(), "raw_evidence": set(), "evidence_urls": set(),
                "drug_mapping_transforms": set(),
            }
        target = grouped[key]
        target["trial_ids"].add(row["trial_id"])
        target["raw_evidence"].add(row["raw_evidence"])
        target["evidence_urls"].add(row["evidence_url"])
        target["drug_mapping_transforms"].add(row["drug_mapping_transform"])
    positive = read_csv("chembl_primekg_positive_candidates.csv")
    positive_keys = {
        (row["drug_node_index"], row["disease_node_index"])
        for row in positive if row["mapping_status"] == "both_mapped_unique"
    }
    output_rows = []
    for key, row in grouped.items():
        for field in ("trial_ids", "raw_evidence", "evidence_urls", "drug_mapping_transforms"):
            row[field] = json.dumps(sorted(row[field]), ensure_ascii=False)
        row["trial_count"] = len(json.loads(row["trial_ids"]))
        row["label"] = "N_strict_rule_based_pending_domain_review"
        row["positive_overlap"] = key in positive_keys
        row["final_use_status"] = "U_pending_context_review" if key in positive_keys else "eligible_candidate_pending_domain_review"
        output_rows.append(row)
    output_rows.sort(key=lambda row: (row["disease_node_name"], row["drug_node_name"]))
    fields = [
        "drug_node_index", "disease_node_index", "drug_node_name", "disease_node_name",
        "trial_ids", "trial_count", "raw_evidence", "evidence_urls", "drug_mapping_transforms",
        "label", "positive_overlap", "final_use_status",
    ]
    output = project_path("data", "interim", "n_strict_rule_based_pairs.csv")
    temporary = output.with_suffix(".csv.tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output_rows)
    temporary.replace(output)

    nonoverlap = [row for row in output_rows if not row["positive_overlap"]]
    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "rule_passed_unique_pairs": len(output_rows),
        "positive_overlapping_pairs": len(output_rows) - len(nonoverlap),
        "nonoverlapping_rule_based_pairs": len(nonoverlap),
        "unique_nonoverlapping_drugs": len({row["drug_node_index"] for row in nonoverlap}),
        "unique_nonoverlapping_diseases": len({row["disease_node_index"] for row in nonoverlap}),
        "e0_minimum_200_independent_n_strict_met": len(nonoverlap) >= 200,
        "domain_review_completed": False,
        "time_history_reconstruction_completed": False,
        "decision": "Do not enter full model comparison; expand N_review adjudication or pivot to a small label-quality/bias study",
    }
    report_path = project_path("reports", "e0_rule_based_gate.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

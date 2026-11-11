import csv
import json
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


def write_csv(name: str, rows: list[dict], fields: list[str]) -> None:
    path = project_path("data", "interim", name)
    temporary = path.with_suffix(".csv.tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> None:
    positive_rows = [row for row in read_csv("chembl_primekg_positive_candidates.csv") if row["mapping_status"] == "both_mapped_unique"]
    negative_rows = read_csv("n_strict_pair_candidates_pending_review.csv")
    positives = {}
    for row in positive_rows:
        key = (row["drug_node_index"], row["disease_node_index"])
        if key not in positives:
            positives[key] = {
                "drug_node_index": key[0], "disease_node_index": key[1],
                "drug_node_name": row["drug_node_name"], "disease_node_name": row["disease_node_name"],
                "chembl_pair_keys": set(), "reference_urls": set(),
            }
        positives[key]["chembl_pair_keys"].add(row["pair_key"])
        positives[key]["reference_urls"].update(json.loads(row["reference_urls"]))
    negatives = {(row["drug_node_index"], row["disease_node_index"]): row for row in negative_rows}

    conflicts = []
    for key in sorted(positives.keys() & negatives.keys()):
        positive = positives[key]
        negative = negatives[key]
        conflicts.append({
            "drug_node_index": key[0], "disease_node_index": key[1],
            "drug_node_name": positive["drug_node_name"], "disease_node_name": positive["disease_node_name"],
            "chembl_pair_keys": json.dumps(sorted(positive["chembl_pair_keys"]), ensure_ascii=False),
            "positive_reference_urls": json.dumps(sorted(positive["reference_urls"]), ensure_ascii=False),
            "negative_trial_ids": negative["trial_ids"], "negative_raw_evidence": negative["raw_evidence"],
            "resolution": "U_pending_context_review",
        })
    conflict_fields = [
        "drug_node_index", "disease_node_index", "drug_node_name", "disease_node_name",
        "chembl_pair_keys", "positive_reference_urls", "negative_trial_ids",
        "negative_raw_evidence", "resolution",
    ]
    write_csv("positive_negative_conflict_review.csv", conflicts, conflict_fields)

    positive_by_disease = Counter(key[1] for key in positives)
    negative_by_disease = Counter(row["disease_node_index"] for row in negative_rows)
    conflict_by_disease = Counter(row["disease_node_index"] for row in conflicts)
    disease_names = {key[1]: value["disease_node_name"] for key, value in positives.items()}
    disease_names.update({row["disease_node_index"]: row["disease_node_name"] for row in negative_rows})
    coverage = []
    for disease in sorted(positive_by_disease.keys() & negative_by_disease.keys()):
        p_count = positive_by_disease[disease]
        n_count = negative_by_disease[disease]
        conflict_count = conflict_by_disease[disease]
        coverage.append({
            "disease_node_index": disease, "disease_node_name": disease_names[disease],
            "positive_candidate_pairs": p_count, "negative_candidate_pairs": n_count,
            "conflicting_pairs": conflict_count,
            "nonconflicting_negative_candidates": n_count - conflict_count,
            "exploratory_coverage_only": True,
        })
    coverage_fields = [
        "disease_node_index", "disease_node_name", "positive_candidate_pairs",
        "negative_candidate_pairs", "conflicting_pairs", "nonconflicting_negative_candidates",
        "exploratory_coverage_only",
    ]
    write_csv("e0_indication_coverage_candidates.csv", coverage, coverage_fields)

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "mapped_positive_pairs": len(positives),
        "negative_pair_candidates": len(negatives),
        "positive_negative_conflicting_pairs": len(conflicts),
        "nonconflicting_negative_pair_candidates": len(negatives) - len(conflicts),
        "indications_with_at_least_one_positive_and_negative_candidate": len(coverage),
        "indications_with_nonconflicting_negative_candidate_after_overlap_exclusion": sum(row["nonconflicting_negative_candidates"] > 0 for row in coverage),
        "coverage_basis": "At least one candidate of each type; final minimum P/N counts are not yet frozen",
        "e0_gate_status": "not_passed_pending_three_layer_review_and_context_conflict_resolution",
        "warning": "Overlapping P and N candidates are context conflicts for review, not automatic label errors or automatic exclusions without documentation",
    }
    report_path = project_path("reports", "e0_conflict_and_coverage.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

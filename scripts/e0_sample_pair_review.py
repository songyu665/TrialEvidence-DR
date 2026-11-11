import csv
import hashlib
import json
import random
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path

SEED = 20261110
SAMPLE_SIZE = 200
csv.field_size_limit(sys.maxsize)


def read_csv(name: str) -> list[dict]:
    with project_path("data", "interim", name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(name: str, rows: list[dict]) -> None:
    path = project_path("data", "interim", name)
    temporary = path.with_suffix(".csv.tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> None:
    source_name = "completed_mapped_n_review_candidates.csv"
    source = project_path("data", "interim", source_name)
    rows = read_csv(source_name)
    grouped = {}
    for row in rows:
        if row["mapping_status"] != "both_mapped_unique":
            continue
        key = (row["drug_node_index"], row["disease_node_index"])
        grouped.setdefault(key, []).append(row)
    high = []
    medium = []
    for key, pair_rows in grouped.items():
        priorities = {priority for row in pair_rows for priority in json.loads(row["review_priorities"])}
        if "high_text_supported_N_review" in priorities:
            high.append(key)
        else:
            medium.append(key)
    high.sort()
    medium.sort()
    selected = high + random.Random(SEED).sample(medium, SAMPLE_SIZE - len(high))
    if len(selected) != SAMPLE_SIZE or len(set(selected)) != SAMPLE_SIZE:
        raise ValueError("Unexpected review sample size")
    selected_set = set(selected)

    trial_info = {row["nct_id"]: row for row in read_csv("completed_primary_outcome_review_queue.csv")}
    pair_rows = []
    for key in selected:
        evidence = grouped[key]
        first = evidence[0]
        trial_ids = sorted({row["trial_id"] for row in evidence})
        pair_rows.append({
            "pair_id": f"{key[0]}|{key[1]}",
            "drug_node_index": key[0], "disease_node_index": key[1],
            "drug_node_name": first["drug_node_name"], "disease_node_name": first["disease_node_name"],
            "raw_intervention_names": json.dumps(sorted({row["raw_intervention_name"] for row in evidence}), ensure_ascii=False),
            "raw_conditions": json.dumps(sorted({row["raw_condition"] for row in evidence}), ensure_ascii=False),
            "trial_ids": json.dumps(trial_ids),
            "arm_groups": json.dumps({trial_id: json.loads(trial_info[trial_id]["arm_groups"]) for trial_id in trial_ids}, ensure_ascii=False),
            "review_priorities": json.dumps(sorted({priority for row in evidence for priority in json.loads(row["review_priorities"])})),
            "drug_mapping_rules": json.dumps(sorted({row["drug_mapping_rule"] for row in evidence})),
            "extraction_fidelity": "", "intervention_attribution": "", "context_aggregation": "",
            "statistical_semantics": "", "clinical_semantics": "", "final_label": "U_pending_review",
            "reviewer": "", "review_date": "", "review_notes": "",
        })
    write_csv("e0_pair_review_sample_200.csv", pair_rows)

    trial_to_pairs = {}
    for pair in pair_rows:
        key = tuple(pair["pair_id"].split("|"))
        for trial_id in json.loads(pair["trial_ids"]):
            trial_to_pairs.setdefault(trial_id, []).append(key)
    analysis_rows = []
    for row in read_csv("completed_primary_analysis_review.csv"):
        if row["review_priority"] not in {"high_text_supported_N_review", "medium_statistical_N_review"} or row["tested_noninferiority"] == "True":
            continue
        for key in trial_to_pairs.get(row["nct_id"], []):
            analysis_rows.append({
                "pair_id": f"{key[0]}|{key[1]}", "nct_id": row["nct_id"],
                "outcome_index": row["outcome_index"], "analysis_index": row["analysis_index"],
                "outcome_title": row["outcome_title"], "analysis_param_type": row["analysis_param_type"],
                "analysis_param_value": row["analysis_param_value"], "p_value": row["p_value"],
                "ci_lower_limit": row["ci_lower_limit"], "ci_upper_limit": row["ci_upper_limit"],
                "group_ids": row["group_ids"], "group_description": row["group_description"],
                "review_priority": row["review_priority"], "design_or_hypothetical_text": row["design_or_hypothetical_text"],
                "actual_negative_result_text": row["actual_negative_result_text"],
                "statistical_comment": row["statistical_comment"], "p_value_comment": row["p_value_comment"],
                "estimate_comment": row["estimate_comment"], "other_analysis_description": row["other_analysis_description"],
                "evidence_url": row["evidence_url"],
            })
    write_csv("e0_pair_review_sample_evidence.csv", analysis_rows)

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source_file": source_name,
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "sampling_unit": "unique PrimeKG drug-disease pair",
        "seed": SEED,
        "high_priority_pairs_all_included": len(high),
        "random_medium_priority_pairs": len(selected) - len(high),
        "total_pairs": len(selected),
        "evidence_analysis_rows": len(analysis_rows),
        "review_status": "not_reviewed",
        "warning": "This is a fixed review sample, not a labeled negative dataset; clinical and statistical semantics remain unassessed",
    }
    report_path = project_path("reports", "e0_pair_review_sample.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

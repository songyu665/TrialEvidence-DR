import csv
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path

FIELDS = [
    "drug_node_index", "disease_node_index", "drug_node_name", "disease_node_name",
    "disease_group_id", "disease_group_name", "trial_ids", "trial_count",
    "source_record_ids", "raw_conditions", "raw_interventions", "raw_evidence",
    "evidence_urls", "label", "review_status",
]


def read_disease_groups() -> dict[str, tuple[str, str]]:
    nodes_path = project_path("data", "raw", "primekg", "nodes.tab")
    with nodes_path.open(encoding="utf-8", newline="") as stream:
        sample = stream.read(4096)
        stream.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        index_to_id = {row["node_index"]: row["node_id"] for row in csv.DictReader(stream, dialect=dialect) if row.get("node_type") == "disease"}
    groups_path = project_path("data", "raw", "primekg", "kg_grouped_diseases_bert_map.tab")
    with groups_path.open(encoding="utf-8", newline="") as stream:
        sample = stream.read(4096)
        stream.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        by_node_id = {row["node_id"]: (row["group_id_bert"], row["group_name_bert"]) for row in csv.DictReader(stream, dialect=dialect)}
    return {index: by_node_id.get(node_id, (f"singleton:{index}", "ungrouped")) for index, node_id in index_to_id.items()}


def main() -> None:
    source = project_path("data", "interim", "mapped_pair_candidates_v2.csv")
    with source.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    disease_groups = read_disease_groups()
    selected = [
        row for row in rows
        if row["candidate_bucket"] == "N_strict_candidate_requires_three_layer_review"
        and row["active_intervention_candidate"] == "True"
        and row["combination_regimen_flag"] != "True"
        and row["mapping_status"] == "both_exact_normalized_unique"
    ]
    grouped = {}
    for row in selected:
        key = (row["drug_id"], row["graph_node_id"])
        group_id, group_name = disease_groups[key[1]]
        if key not in grouped:
            grouped[key] = {
                "drug_node_index": key[0],
                "disease_node_index": key[1],
                "drug_node_name": row["drug_node_name"],
                "disease_node_name": row["disease_node_name"],
                "disease_group_id": group_id,
                "disease_group_name": group_name,
                "trial_ids": set(),
                "source_record_ids": set(),
                "raw_conditions": set(),
                "raw_interventions": set(),
                "raw_evidence": set(),
                "evidence_urls": set(),
            }
        target = grouped[key]
        target["trial_ids"].add(row["trial_id"])
        target["source_record_ids"].add(row["record_id"])
        target["raw_conditions"].add(row["raw_condition"])
        target["raw_interventions"].add(row["raw_intervention_name"])
        target["raw_evidence"].add(row["raw_evidence"])
        target["evidence_urls"].add(row["evidence_url"])

    output_rows = []
    for row in grouped.values():
        row["trial_count"] = len(row["trial_ids"])
        for field in ("trial_ids", "source_record_ids", "raw_conditions", "raw_interventions", "raw_evidence", "evidence_urls"):
            row[field] = json.dumps(sorted(row[field]), ensure_ascii=False)
        row["label"] = "N_strict_pair_candidate_pending_review"
        row["review_status"] = "pending_extraction_attribution_and_context_review"
        output_rows.append(row)
    output_rows.sort(key=lambda row: (row["disease_node_name"], row["drug_node_name"]))

    output = project_path("data", "interim", "n_strict_pair_candidates_pending_review.csv")
    temporary = output.with_suffix(".csv.tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output_rows)
    temporary.replace(output)

    per_disease = Counter(row["disease_node_name"] for row in output_rows)
    per_group = Counter(row["disease_group_id"] for row in output_rows)
    nodes_per_group = {}
    for row in output_rows:
        nodes_per_group.setdefault(row["disease_group_id"], set()).add(row["disease_node_index"])
    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source_rows": len(selected),
        "unique_pair_candidates": len(output_rows),
        "unique_drugs": len({row["drug_node_index"] for row in output_rows}),
        "unique_diseases": len({row["disease_node_index"] for row in output_rows}),
        "unique_disease_groups": len(per_group),
        "disease_groups_with_multiple_mapped_nodes": sum(len(nodes) > 1 for nodes in nodes_per_group.values()),
        "pairs_with_multiple_trials": sum(row["trial_count"] > 1 for row in output_rows),
        "disease_candidate_count_distribution": dict(Counter(per_disease.values())),
        "disease_group_candidate_count_distribution": dict(Counter(per_group.values())),
        "e0_gate_status": "not_evaluated_pending_manual_review_positive_mapping_and_conflict_checks",
        "warning": "Groups are for entity-group audits and split sensitivity, not automatic clinical-label merging; candidate pairs remain unvalidated",
    }
    report_path = project_path("reports", "e0_aggregated_n_strict_candidates.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

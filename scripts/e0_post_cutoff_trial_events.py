import csv
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path

CUTOFF = "2023-12-01"
csv.field_size_limit(sys.maxsize)


def read_csv(name: str) -> list[dict]:
    with project_path("data", "interim", name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    queue = {row["nct_id"]: row for row in read_csv("completed_primary_outcome_review_queue.csv")}
    historical_negative = {(row["drug_node_index"], row["disease_node_index"]) for row in read_csv("historical_20231201_n_rule_pairs.csv")}
    current_negative = read_csv("n_strict_rule_based_pairs.csv")
    historical_positive = {(row["drug_node_index"], row["disease_node_index"]) for row in read_csv("historical_chembl33_positive_pairs.csv")}

    output_rows = []
    for row in read_csv("completed_mapped_n_review_candidates.csv"):
        trial = queue[row["trial_id"]]
        posted = trial["results_first_posted"]
        if not posted or posted <= CUTOFF:
            continue
        key = (row["drug_node_index"], row["disease_node_index"])
        output_rows.append({
            "event_type": "post_cutoff_completed_trial_result_N_review",
            "trial_id": row["trial_id"], "drug_node_index": key[0], "disease_node_index": key[1],
            "drug_node_name": row["drug_node_name"], "disease_node_name": row["disease_node_name"],
            "event_date": posted, "date_type": "results_first_posted_date",
            "review_priorities": row["review_priorities"], "evidence_url": row["evidence_url"],
            "pre_cutoff_positive_pair": key in historical_positive,
            "pre_cutoff_rule_negative_pair": key in historical_negative,
            "label": "N_review_pending_event_adjudication",
        })

    for row in current_negative:
        key = (row["drug_node_index"], row["disease_node_index"])
        if key in historical_negative:
            continue
        output_rows.append({
            "event_type": "later_futility_text_not_visible_in_2023_snapshot",
            "trial_id": ";".join(json.loads(row["trial_ids"])),
            "drug_node_index": key[0], "disease_node_index": key[1],
            "drug_node_name": row["drug_node_name"], "disease_node_name": row["disease_node_name"],
            "event_date": "", "date_type": "exact_field_publication_date_not_recovered",
            "review_priorities": "current_N_rule", "evidence_url": row["evidence_urls"],
            "pre_cutoff_positive_pair": key in historical_positive,
            "pre_cutoff_rule_negative_pair": False,
            "label": "N_review_pending_event_and_date_adjudication",
        })

    output = project_path("data", "interim", "post_cutoff_trial_event_candidates.csv")
    with output.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(), "cutoff": CUTOFF,
        "event_type_counts": dict(Counter(row["event_type"] for row in output_rows)),
        "unique_event_pairs": len({(row["drug_node_index"], row["disease_node_index"]) for row in output_rows}),
        "completed_result_events_by_year": dict(Counter(row["event_date"][:4] for row in output_rows if row["event_type"] == "post_cutoff_completed_trial_result_N_review")),
        "events_with_exact_results_first_posted_date": sum(bool(row["event_date"]) for row in output_rows),
        "events_without_exact_field_publication_date": sum(not row["event_date"] for row in output_rows),
        "event_adjudication_completed": False,
        "interpretation": "Post-cutoff evidence candidates, not confirmed positive or negative outcomes",
    }
    with project_path("reports", "e0_post_cutoff_trial_events.json").open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

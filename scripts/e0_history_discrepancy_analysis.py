import csv
import gzip
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path

CUTOFF = "2023-12-01"


def read_csv(name: str) -> list[dict]:
    with project_path("data", "interim", name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    current_pre = read_csv("n_strict_pair_candidates_pending_review.csv")
    historical_pre = read_csv("historical_20231201_n_strict_candidates.csv")
    current_rule = read_csv("n_strict_rule_based_pairs.csv")
    historical_rule = read_csv("historical_20231201_n_rule_pairs.csv")
    historical_pre_keys = {(row["drug_node_index"], row["disease_node_index"]) for row in historical_pre}
    historical_rule_keys = {(row["drug_node_index"], row["disease_node_index"]) for row in historical_rule}

    needed_trials = set()
    for row in current_pre + current_rule:
        needed_trials.update(json.loads(row["trial_ids"]))
    trial_dates = {}
    trial_dir = project_path("data", "raw", "clinicaltrials", "terminated_phase23")
    state = json.loads((trial_dir / "download_checkpoint.json").read_text())
    for page in range(1, state["pages_downloaded"] + 1):
        with gzip.open(trial_dir / f"page_{page:06d}.json.gz", "rt", encoding="utf-8") as stream:
            for study in json.load(stream)["studies"]:
                protocol = study["protocolSection"]
                identifier = protocol["identificationModule"]["nctId"]
                if identifier in needed_trials:
                    status = protocol.get("statusModule", {})
                    trial_dates[identifier] = {
                        "study_first_posted": status.get("studyFirstPostDateStruct", {}).get("date"),
                        "results_first_posted": status.get("resultsFirstPostDateStruct", {}).get("date"),
                        "last_update_posted": status.get("lastUpdatePostDateStruct", {}).get("date"),
                    }

    def discrepancy(rows: list[dict], historical_keys: set, level: str):
        output = []
        for row in rows:
            key = (row["drug_node_index"], row["disease_node_index"])
            if key in historical_keys:
                continue
            trials = json.loads(row["trial_ids"])
            dates = [trial_dates[trial] for trial in trials]
            first_dates = sorted({item["study_first_posted"] for item in dates if item["study_first_posted"]})
            if first_dates and all(date > CUTOFF for date in first_dates):
                category = "all_trial_records_first_posted_after_snapshot"
            elif any(date and date <= CUTOFF for date in first_dates):
                category = "existed_by_snapshot_but_text_mapping_or_rule_changed"
            else:
                category = "first_posted_date_missing"
            output.append({
                "level": level, "drug_node_index": key[0], "disease_node_index": key[1],
                "drug_node_name": row["drug_node_name"], "disease_node_name": row["disease_node_name"],
                "trial_ids": json.dumps(trials), "study_first_posted_dates": json.dumps(first_dates),
                "last_update_posted_dates": json.dumps(sorted({item["last_update_posted"] for item in dates if item["last_update_posted"]})),
                "discrepancy_category": category,
            })
        return output

    output_rows = discrepancy(current_pre, historical_pre_keys, "pre_review_pair") + discrepancy(current_rule, historical_rule_keys, "rule_passed_pair")
    output = project_path("data", "interim", "history_candidate_discrepancies.csv")
    with output.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "snapshot_date": CUTOFF,
        "category_counts_by_level": {
            level: dict(Counter(row["discrepancy_category"] for row in output_rows if row["level"] == level))
            for level in ("pre_review_pair", "rule_passed_pair")
        },
        "warning": "Date categories explain potential sources of discrepancy but do not prove when whyStopped or a mapping became publicly available",
    }
    report_path = project_path("reports", "e0_history_discrepancies.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    for row in output_rows:
        if row["level"] == "rule_passed_pair":
            print(json.dumps(row, ensure_ascii=False))


if __name__ == "__main__":
    main()

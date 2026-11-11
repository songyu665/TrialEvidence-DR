import csv
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path

REVIEW_BUCKETS = {
    "N_strict_candidate_requires_three_layer_review",
    "conflicting_reason_candidates_exclude_pending_review",
    "N_review_candidate_not_automatic_negative",
}
ACTIVE_TYPES = {"DRUG", "BIOLOGICAL", "COMBINATION_PRODUCT"}
NONACTIVE = re.compile(r"\bplacebos?\b|\bvehicle\b|\bsham\b|normal saline|standard of care|best supportive care|physician'?s choice|standard treatment|best medical therapy", re.I)
COMBINATION = re.compile(r"\s*\+\s*|\bcombined with\b|\bin combination with\b|\bplus\b|\s+and\s+|;", re.I)

FIELDS = [
    "record_id", "trial_id", "drug_id", "disease_id", "raw_intervention_name",
    "raw_condition", "intervention_type", "raw_evidence", "intervention_attribution",
    "indication_match", "comparator_type", "trial_purpose", "stop_reason",
    "result_available", "candidate_bucket", "active_intervention_candidate",
    "combination_regimen_flag", "canonical_context", "graph_node_id", "mapping_status",
    "aggregation_rule", "label", "exclusion_reason", "reviewer", "review_date",
    "evidence_url", "review_notes",
]


def main() -> None:
    source = project_path("data", "interim", "trial_label_candidates.csv")
    with source.open(encoding="utf-8-sig", newline="") as stream:
        trials = list(csv.DictReader(stream))
    rows = []
    skipped = Counter()
    for trial in trials:
        if trial["candidate_bucket"] not in REVIEW_BUCKETS:
            continue
        interventions = json.loads(trial["interventions"])
        conditions = json.loads(trial["conditions"])
        if not interventions or not conditions:
            skipped["missing_intervention_or_condition"] += 1
            continue
        for intervention_index, intervention in enumerate(interventions):
            if intervention.get("type") not in ACTIVE_TYPES:
                skipped["non_drug_intervention"] += 1
                continue
            name = (intervention.get("name") or "").strip()
            if not name:
                skipped["missing_intervention_name"] += 1
                continue
            active = not NONACTIVE.search(name)
            combination = bool(COMBINATION.search(name))
            for condition_index, condition in enumerate(conditions):
                condition = condition.strip()
                if not condition:
                    continue
                rows.append(
                    {
                        "record_id": f"{trial['nct_id']}:I{intervention_index}:C{condition_index}",
                        "trial_id": trial["nct_id"],
                        "drug_id": "",
                        "disease_id": "",
                        "raw_intervention_name": name,
                        "raw_condition": condition,
                        "intervention_type": intervention.get("type"),
                        "raw_evidence": trial["why_stopped"],
                        "intervention_attribution": "pending",
                        "indication_match": "pending",
                        "comparator_type": "",
                        "trial_purpose": "",
                        "stop_reason": trial["why_stopped"],
                        "result_available": trial["results_available"],
                        "candidate_bucket": trial["candidate_bucket"],
                        "active_intervention_candidate": active,
                        "combination_regimen_flag": combination,
                        "canonical_context": "",
                        "graph_node_id": "",
                        "mapping_status": "unmapped",
                        "aggregation_rule": "pending",
                        "label": "U_pending_review",
                        "exclusion_reason": "" if active else "nonactive_comparator_or_placebo",
                        "reviewer": "",
                        "review_date": "",
                        "evidence_url": f"https://clinicaltrials.gov/study/{trial['nct_id']}",
                        "review_notes": "",
                    }
                )

    output = project_path("data", "interim", "label_adjudication_queue.csv")
    temporary = output.with_suffix(".csv.tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(output)

    download_summary_path = project_path("data", "raw", "clinicaltrials", "terminated_phase23", "download_summary.json")
    with download_summary_path.open(encoding="utf-8") as stream:
        download_complete = json.load(stream)["complete_query_download"]
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source_file": str(source.relative_to(project_path())),
        "source_download_is_complete": download_complete,
        "review_rows": len(rows),
        "unique_trials": len({row["trial_id"] for row in rows}),
        "unique_raw_intervention_names": len({row["raw_intervention_name"].casefold() for row in rows}),
        "unique_raw_conditions": len({row["raw_condition"].casefold() for row in rows}),
        "active_intervention_rows": sum(row["active_intervention_candidate"] for row in rows),
        "combination_regimen_rows": sum(row["combination_regimen_flag"] for row in rows),
        "skipped": dict(skipped),
        "warning": "Rows are unadjudicated candidate pairs, not confirmed N_strict labels or independent experimental units",
    }
    report = project_path("reports", "e0_adjudication_queue_summary.json")
    with report.open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

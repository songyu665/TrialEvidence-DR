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


def read_csv(name: str) -> list[dict]:
    with project_path("data", "interim", name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    rows = read_csv("mapped_pair_candidates_v2.csv")
    selected = [
        row for row in rows
        if row["candidate_bucket"] == "N_strict_candidate_requires_three_layer_review"
        and row["active_intervention_candidate"] == "True"
        and row["combination_regimen_flag"] != "True"
        and row["mapping_status"] == "both_exact_normalized_unique"
    ]
    trial_ids = {row["trial_id"] for row in selected}
    raw_trials = {}
    trial_dir = project_path("data", "raw", "clinicaltrials", "terminated_phase23")
    state = json.loads((trial_dir / "download_checkpoint.json").read_text())
    for page in range(1, state["pages_downloaded"] + 1):
        with gzip.open(trial_dir / f"page_{page:06d}.json.gz", "rt", encoding="utf-8") as stream:
            for study in json.load(stream)["studies"]:
                identifier = study["protocolSection"]["identificationModule"]["nctId"]
                if identifier in trial_ids:
                    raw_trials[identifier] = study

    all_rows_by_trial = {}
    selected_by_trial = {}
    for row in rows:
        if row["candidate_bucket"] == "N_strict_candidate_requires_three_layer_review":
            all_rows_by_trial.setdefault(row["trial_id"], []).append(row)
    for row in selected:
        selected_by_trial.setdefault(row["trial_id"], []).append(row)
    output_rows = []
    for trial_id, selected_rows in sorted(selected_by_trial.items()):
        study = raw_trials[trial_id]
        protocol = study.get("protocolSection", {})
        arm_groups = protocol.get("armsInterventionsModule", {}).get("armGroups", [])
        arm_types_by_intervention = {}
        for arm in arm_groups:
            for name in arm.get("interventionNames", []):
                key = name.split(":", 1)[-1].strip().casefold()
                arm_types_by_intervention.setdefault(key, set()).add(arm.get("type", ""))
        active_rows = [row for row in all_rows_by_trial[trial_id] if row["active_intervention_candidate"] == "True"]
        active_drugs = {row["drug_id"] for row in active_rows}
        active_names = {row["raw_intervention_name"].casefold() for row in active_rows}
        disease_nodes = {row["graph_node_id"] for row in active_rows}
        active_comparator_present = any(arm.get("type") == "ACTIVE_COMPARATOR" for arm in arm_groups)

        pair_groups = {}
        for row in selected_rows:
            key = (row["drug_id"], row["graph_node_id"])
            pair_groups.setdefault(key, []).append(row)
        for (drug, disease), pair_rows in sorted(pair_groups.items()):
            reasons = []
            if len(active_drugs) != 1:
                reasons.append("multiple_or_unresolved_active_interventions")
            if len(disease_nodes) != 1:
                reasons.append("multiple_or_unresolved_disease_contexts")
            if any(row["disease_mapping_status"] != "exact_normalized_unique" for row in active_rows):
                reasons.append("not_all_trial_conditions_uniquely_mapped")
            if active_comparator_present:
                reasons.append("active_comparator_requires_clinical_review")
            target_types = set()
            for name in {row["raw_intervention_name"].split(":", 1)[-1].strip().casefold() for row in pair_rows}:
                target_types.update(arm_types_by_intervention.get(name, set()))
            if not target_types:
                reasons.append("arm_attribution_missing")
            elif target_types != {"EXPERIMENTAL"}:
                reasons.append("target_not_exclusively_experimental")
            first = pair_rows[0]
            passed = not reasons
            output_rows.append({
                "trial_id": trial_id,
                "drug_node_index": drug,
                "disease_node_index": disease,
                "drug_node_name": first["drug_node_name"],
                "disease_node_name": first["disease_node_name"],
                "raw_intervention_name": first["raw_intervention_name"],
                "raw_condition": first["raw_condition"],
                "raw_evidence": first["raw_evidence"],
                "drug_mapping_transform": first.get("drug_mapping_transform", ""),
                "arm_types": json.dumps(sorted(target_types)),
                "active_intervention_names": json.dumps(sorted(active_names), ensure_ascii=False),
                "rule_review_status": "passed_conservative_rules" if passed else "excluded_or_needs_review",
                "exclusion_reasons": json.dumps(reasons),
                "label": "N_strict_rule_based_pending_domain_review" if passed else "U_pending_review",
                "reviewer": "automated_conservative_rule_v1",
                "review_date": datetime.now(timezone.utc).date().isoformat(),
                "evidence_url": first["evidence_url"],
            })

    fields = list(output_rows[0])
    output = project_path("data", "interim", "n_strict_conservative_rule_review.csv")
    temporary = output.with_suffix(".csv.tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output_rows)
    temporary.replace(output)

    passed = [row for row in output_rows if row["rule_review_status"] == "passed_conservative_rules"]
    reason_counts = Counter(reason for row in output_rows for reason in json.loads(row["exclusion_reasons"]))
    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "input_trial_pair_records": len(output_rows),
        "passed_trial_pair_records": len(passed),
        "passed_unique_pair_candidates": len({(row["drug_node_index"], row["disease_node_index"]) for row in passed}),
        "passed_unique_drugs": len({row["drug_node_index"] for row in passed}),
        "passed_unique_diseases": len({row["disease_node_index"] for row in passed}),
        "exclusion_reason_counts": dict(reason_counts),
        "review_scope": "Text signal, structured intervention attribution and unique context aggregation rules",
        "clinical_validation": "not_performed_no_domain_reviewer",
        "warning": "Rule-passed pairs are frozen rule-based evidence candidates without clinical adjudication; they are not clinical gold-standard negatives",
    }
    report_path = project_path("reports", "e0_conservative_rule_review.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

import csv
import json
import sys
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from trialevidence.paths import project_path
from e0_build_adjudication_queue import ACTIVE_TYPES, NONACTIVE
from e0_map_chembl_primekg import read_nodes, unique_matches

csv.field_size_limit(sys.maxsize)


def read_csv(name: str) -> list[dict]:
    with project_path("data", "interim", name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def rows_from_zip(z: zipfile.ZipFile, name: str):
    with z.open(name) as stream:
        wrapper = (line.decode("utf-8-sig" if index == 0 else "utf-8") for index, line in enumerate(stream))
        yield from csv.DictReader(wrapper, delimiter="|")


def main() -> None:
    candidates = read_csv("historical_20231201_n_strict_candidates.csv")
    candidate_trials = {trial for row in candidates for trial in json.loads(row["trial_ids"])}
    known = {"drug": {}, "disease": {}}
    for row in read_csv("observed_trial_name_mapping.csv"):
        if row["node_index"]:
            known[row["node_type"]][row["raw_name"]] = (row["node_index"], row["node_name"])
    nodes = read_nodes()

    def map_name(raw: str, node_type: str):
        if raw in known[node_type]:
            return known[node_type][raw]
        matches = unique_matches({raw}, nodes[node_type])
        return (matches[0]["node_index"], matches[0]["node_name"]) if len(matches) == 1 else (None, None)

    trials = {identifier: {"interventions": [], "conditions": [], "groups": []} for identifier in candidate_trials}
    intervention_by_id = {}
    with zipfile.ZipFile(project_path("data", "raw", "aact", "20231201_aact.zip")) as z:
        for row in rows_from_zip(z, "interventions.txt"):
            if row["nct_id"] in candidate_trials:
                item = {"id": row["id"], "type": row.get("intervention_type", "").upper(), "name": row.get("name", "")}
                trials[row["nct_id"]]["interventions"].append(item)
                intervention_by_id[row["id"]] = item
        for row in rows_from_zip(z, "conditions.txt"):
            if row["nct_id"] in candidate_trials:
                trials[row["nct_id"]]["conditions"].append(row.get("name", ""))
        for row in rows_from_zip(z, "design_groups.txt"):
            if row["nct_id"] in candidate_trials:
                trials[row["nct_id"]]["groups"].append({"id": row["id"], "type": row.get("group_type", "")})
        for row in rows_from_zip(z, "design_group_interventions.txt"):
            if row["nct_id"] in candidate_trials and row["intervention_id"] in intervention_by_id:
                for group in trials[row["nct_id"]]["groups"]:
                    if group["id"] == row["design_group_id"]:
                        group.setdefault("intervention_ids", set()).add(row["intervention_id"])

    reviewed = []
    for nct_id in sorted(candidate_trials):
        trial = trials[nct_id]
        active = [item for item in trial["interventions"] if item["type"] in ACTIVE_TYPES and item["name"] and not NONACTIVE.search(item["name"])]
        active_drugs = {map_name(item["name"], "drug")[0] for item in active}
        disease_nodes = {map_name(name, "disease")[0] for name in trial["conditions"]}
        active_comparator = any(group["type"] == "Active Comparator" for group in trial["groups"])
        types_by_intervention = {}
        for group in trial["groups"]:
            for intervention_id in group.get("intervention_ids", set()):
                types_by_intervention.setdefault(intervention_id, set()).add(group["type"])
        trial_candidates = [row for row in candidates if nct_id in json.loads(row["trial_ids"])]
        for candidate in trial_candidates:
            target_types = set()
            for item in active:
                if map_name(item["name"], "drug")[0] == candidate["drug_node_index"]:
                    target_types.update(types_by_intervention.get(item["id"], set()))
            reasons = []
            if len(active_drugs) != 1:
                reasons.append("multiple_or_unresolved_active_interventions")
            if len(disease_nodes) != 1 or None in disease_nodes:
                reasons.append("multiple_or_unresolved_disease_contexts")
            if active_comparator:
                reasons.append("active_comparator_requires_clinical_review")
            if not target_types:
                reasons.append("arm_attribution_missing")
            elif target_types != {"Experimental"}:
                reasons.append("target_not_exclusively_experimental")
            reviewed.append({
                "trial_id": nct_id, "drug_node_index": candidate["drug_node_index"],
                "disease_node_index": candidate["disease_node_index"],
                "drug_node_name": candidate["drug_node_name"], "disease_node_name": candidate["disease_node_name"],
                "raw_evidence": candidate["raw_evidence"], "rule_review_status": "passed_conservative_rules" if not reasons else "excluded_or_needs_review",
                "exclusion_reasons": json.dumps(reasons), "snapshot": "AACT_2023-12-01",
            })

    passed = [row for row in reviewed if row["rule_review_status"] == "passed_conservative_rules"]
    passed_pairs = {}
    for row in passed:
        key = (row["drug_node_index"], row["disease_node_index"])
        if key not in passed_pairs:
            passed_pairs[key] = {
                "drug_node_index": key[0], "disease_node_index": key[1],
                "drug_node_name": row["drug_node_name"], "disease_node_name": row["disease_node_name"],
                "trial_ids": set(), "raw_evidence": set(), "snapshot": "AACT_2023-12-01",
                "label": "historical_N_rule_pending_domain_review",
            }
        passed_pairs[key]["trial_ids"].add(row["trial_id"])
        passed_pairs[key]["raw_evidence"].update(json.loads(row["raw_evidence"]))
    output_rows = []
    for row in passed_pairs.values():
        row["trial_ids"] = json.dumps(sorted(row["trial_ids"]))
        row["raw_evidence"] = json.dumps(sorted(row["raw_evidence"]), ensure_ascii=False)
        output_rows.append(row)
    output = project_path("data", "interim", "historical_20231201_n_rule_pairs.csv")
    with output.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    current_pairs = {(row["drug_node_index"], row["disease_node_index"]) for row in read_csv("n_strict_rule_based_pairs.csv")}
    historical_pairs = set(passed_pairs)
    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "historical_mapped_pair_candidates": len(candidates),
        "historical_trial_pair_records_reviewed": len(reviewed),
        "historical_rule_passed_unique_pairs": len(historical_pairs),
        "current_rule_pairs_recoverable_in_history": len(current_pairs & historical_pairs),
        "current_rule_pairs_total": len(current_pairs),
        "historical_rule_pair_recoverability_fraction": len(current_pairs & historical_pairs) / len(current_pairs),
        "exclusion_reason_counts": dict(Counter(reason for row in reviewed for reason in json.loads(row["exclusion_reasons"]))),
        "time_experiment_allowed": False,
        "warning": "Historical rule labels are not clinically adjudicated; positive-label reconstruction and information-availability audits remain separate",
    }
    report_path = project_path("reports", "e0_aact_20231201_conservative_review.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

import csv
import gzip
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from trialevidence.paths import project_path
from e0_build_adjudication_queue import ACTIVE_TYPES, NONACTIVE
from e0_map_chembl_primekg import read_nodes, unique_matches

SEED = 20261110
DATASETS = {"terminated": "terminated_phase23", "completed_results": "completed_results_phase23_with_comments"}


def read_csv(name: str) -> list[dict]:
    with project_path("data", "interim", name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    known = {"drug": {}, "disease": {}}
    for row in read_csv("observed_trial_name_mapping.csv"):
        if row["node_index"]:
            known[row["node_type"]][row["raw_name"]] = (row["node_index"], row["node_name"])
    nodes = read_nodes()

    def map_name(raw: str, node_type: str):
        if raw in known[node_type]:
            return known[node_type][raw][0]
        matches = unique_matches({raw}, nodes[node_type])
        return matches[0]["node_index"] if len(matches) == 1 else None

    memberships = []
    global_counts = {"drug": Counter(), "disease": Counter()}
    for dataset, directory_name in DATASETS.items():
        directory = project_path("data", "raw", "clinicaltrials", directory_name)
        state = json.loads((directory / "download_checkpoint.json").read_text())
        for page in range(1, state["pages_downloaded"] + 1):
            with gzip.open(directory / f"page_{page:06d}.json.gz", "rt", encoding="utf-8") as stream:
                studies = json.load(stream)["studies"]
            for study in studies:
                protocol = study.get("protocolSection", {})
                nct_id = protocol["identificationModule"]["nctId"]
                interventions = protocol.get("armsInterventionsModule", {}).get("interventions", [])
                conditions = protocol.get("conditionsModule", {}).get("conditions", [])
                drug_nodes = set()
                disease_nodes = set()
                for intervention in interventions:
                    name = (intervention.get("name") or "").strip()
                    if intervention.get("type") not in ACTIVE_TYPES or not name or NONACTIVE.search(name):
                        continue
                    mapped = map_name(name, "drug")
                    if mapped:
                        drug_nodes.add(mapped)
                for condition in conditions:
                    mapped = map_name(condition, "disease")
                    if mapped:
                        disease_nodes.add(mapped)
                memberships.append({"nct_id": nct_id, "drug_nodes": drug_nodes, "disease_nodes": disease_nodes})
                for node in drug_nodes:
                    global_counts["drug"][node] += 1
                for node in disease_nodes:
                    global_counts["disease"][node] += 1
    membership_by_trial = {row["nct_id"]: row for row in memberships}

    positive_keys = {
        (int(row["drug_node_index"]), int(row["disease_node_index"]))
        for row in read_csv("chembl_primekg_positive_candidates.csv")
        if row["mapping_status"] == "both_mapped_unique"
    }
    source_trials = defaultdict(set)
    rule_keys = set()
    for row in read_csv("n_strict_rule_based_pairs.csv"):
        key = (int(row["drug_node_index"]), int(row["disease_node_index"]))
        rule_keys.add(key)
        source_trials[key].update(json.loads(row["trial_ids"]))
    for row in read_csv("completed_mapped_n_review_candidates.csv"):
        if row["mapping_status"] == "both_mapped_unique":
            key = (int(row["drug_node_index"]), int(row["disease_node_index"]))
            if key not in positive_keys and key not in rule_keys:
                source_trials[key].add(row["trial_id"])

    original_rows = read_csv("selection_bias_pair_features.csv")
    adjusted_rows = []
    for row in original_rows:
        key = (int(row["drug_node_index"]), int(row["disease_node_index"]))
        excluded = source_trials.get(key, set())
        for node_type, index in (("drug", key[0]), ("disease", key[1])):
            original = int(row[f"{node_type}_observed_trials"])
            source_count = sum(str(index) in membership_by_trial[trial][f"{node_type}_nodes"] for trial in excluded if trial in membership_by_trial)
            adjusted = original - source_count
            if adjusted < 0:
                raise ValueError("Negative leave-source-out count")
            row[f"{node_type}_observed_trials_original"] = original
            row[f"{node_type}_observed_trials"] = adjusted
            row[f"{node_type}_label_source_trials_removed"] = source_count
        row["leave_label_source_out"] = True
        adjusted_rows.append(row)

    output = project_path("data", "interim", "selection_bias_pair_features_leave_source_out.csv")
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(adjusted_rows[0]))
        writer.writeheader()
        writer.writerows(adjusted_rows)

    feature_sets = {
        "degrees_only": ["drug_auxiliary_degree", "disease_auxiliary_degree"],
        "trials_only_adjusted": ["drug_observed_trials", "disease_observed_trials"],
        "all_features_adjusted": ["drug_auxiliary_degree", "disease_auxiliary_degree", "drug_observed_trials", "disease_observed_trials"],
    }
    result_rows = []
    for task in ("P", "N_rule", "N_review"):
        subset = [row for row in adjusted_rows if row["evidence_set"] in {task, "U"}]
        y = np.array([row["evidence_set"] == task for row in subset], dtype=int)
        groups = np.array([row["disease_node_index"] for row in subset])
        for feature_name, features in feature_sets.items():
            x = np.log1p(np.array([[float(row[feature]) for feature in features] for row in subset]))
            oof = np.zeros(len(y))
            for train, test in StratifiedGroupKFold(5, shuffle=True, random_state=SEED).split(x, y, groups):
                model = make_pipeline(StandardScaler(), LogisticRegression(class_weight="balanced", max_iter=2000, random_state=SEED))
                model.fit(x[train], y[train])
                oof[test] = model.predict_proba(x[test])[:, 1]
            result_rows.append({
                "task": f"{task}_vs_U", "feature_set": feature_name,
                "n": len(y), "positive_prevalence": float(y.mean()),
                "oof_auprc": average_precision_score(y, oof), "oof_auroc": roc_auc_score(y, oof),
                "prediction_target": "evidence status after removing each pair's label-source trials",
            })
            print(json.dumps(result_rows[-1]), flush=True)

    with project_path("reports", "tables", "f3_leave_label_source_out_results.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(result_rows[0]))
        writer.writeheader()
        writer.writerows(result_rows)

    with project_path("reports", "tables", "f3_evidence_status_model_results.csv").open(encoding="utf-8", newline="") as stream:
        original_model_rows = list(csv.DictReader(stream))
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.6))
    for axis, task in zip(axes, ("P", "N_rule", "N_review")):
        original = next(row for row in original_model_rows if row["task"] == f"{task}_vs_U" and row["group_mode"] == "disease" and row["feature_set"] == "all_features" and row["model"] == "logistic")
        adjusted = next(row for row in result_rows if row["task"] == f"{task}_vs_U" and row["feature_set"] == "all_features_adjusted")
        axis.bar([0, 1], [float(original["oof_auprc"]), adjusted["oof_auprc"]], color=["#B8B8B8", "#3B6E8F"])
        axis.set_xticks([0, 1], ["Original", "Leave-source-out"])
        axis.set_ylim(0, 1)
        axis.set_title(f"{task} vs U")
        axis.set_ylabel("Disease-grouped AUPRC")
    fig.suptitle("Effect of removing label-source trial records")
    fig.text(0.01, -0.04, "Sensitivity analysis for same-source information cycling; not independent external validation.", fontsize=8)
    base = project_path("reports", "figures", "f3_leave_label_source_out")
    plt.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "pairs": len(adjusted_rows),
        "pairs_with_any_removed_source_trial": sum(int(row["drug_label_source_trials_removed"]) + int(row["disease_label_source_trials_removed"]) > 0 for row in adjusted_rows),
        "source_trial_mapping": "N_rule and N_review pairs; P uses regulatory/label selection and was not adjusted",
        "scope": "Leave-label-source-trials-out sensitivity, not independent external validation",
    }
    with project_path("reports", "f3_leave_label_source_out.json").open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

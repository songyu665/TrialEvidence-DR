import csv
import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path

SEED = 20261110
BOOTSTRAPS = 200
FEATURE_SETS = {
    "degrees_only": ["drug_auxiliary_degree", "disease_auxiliary_degree"],
    "trials_only": ["drug_observed_trials", "disease_observed_trials"],
    "all_features": ["drug_auxiliary_degree", "disease_auxiliary_degree", "drug_observed_trials", "disease_observed_trials"],
}
TASKS = [("P", "U"), ("N_rule", "U"), ("N_review", "U")]


def bootstrap_auprc(y: np.ndarray, scores: np.ndarray, groups: np.ndarray) -> list[float]:
    grouped = defaultdict(list)
    for index, group in enumerate(groups):
        grouped[group].append(index)
    group_list = list(grouped)
    indices = [np.array(grouped[group]) for group in group_list]
    rng = np.random.default_rng(SEED)
    estimates = []
    for _ in range(BOOTSTRAPS):
        sampled = rng.integers(0, len(group_list), len(group_list))
        selected = np.concatenate([indices[index] for index in sampled])
        if len(np.unique(y[selected])) < 2:
            continue
        estimates.append(average_precision_score(y[selected], scores[selected]))
    return [float(np.percentile(estimates, 2.5)), float(np.percentile(estimates, 97.5))]


def main() -> None:
    with project_path("data", "interim", "selection_bias_pair_features.csv").open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    output_rows = []
    for task_index, (positive, negative) in enumerate(TASKS):
        subset = [row for row in rows if row["evidence_set"] in {positive, negative}]
        y = np.array([row["evidence_set"] == positive for row in subset], dtype=int)
        for group_mode in ("disease", "drug"):
            groups = np.array([row[f"{group_mode}_node_index"] for row in subset])
            for feature_name, features in FEATURE_SETS.items():
                x = np.log1p(np.array([[float(row[feature]) for feature in features] for row in subset]))
                for model_name in ("logistic", "random_forest"):
                    oof = np.zeros(len(y))
                    fold_auprc = []
                    fold_auroc = []
                    for fold, (train, test) in enumerate(StratifiedGroupKFold(5, shuffle=True, random_state=SEED).split(x, y, groups)):
                        if model_name == "logistic":
                            model = make_pipeline(StandardScaler(), LogisticRegression(class_weight="balanced", max_iter=2000, random_state=SEED))
                        else:
                            model = RandomForestClassifier(n_estimators=200, class_weight="balanced", min_samples_leaf=3, n_jobs=2, random_state=SEED + fold)
                        model.fit(x[train], y[train])
                        scores = model.predict_proba(x[test])[:, 1]
                        oof[test] = scores
                        fold_auprc.append(average_precision_score(y[test], scores))
                        fold_auroc.append(roc_auc_score(y[test], scores))
                    result = {
                        "task": f"{positive}_vs_{negative}", "group_mode": group_mode,
                        "feature_set": feature_name, "model": model_name,
                        "n": len(y), "positive_prevalence": float(y.mean()),
                        "oof_auprc": average_precision_score(y, oof),
                        "oof_auroc": roc_auc_score(y, oof),
                        "fold_auprc_mean": float(np.mean(fold_auprc)),
                        "fold_auprc_min": float(np.min(fold_auprc)), "fold_auprc_max": float(np.max(fold_auprc)),
                        "fold_auroc_mean": float(np.mean(fold_auroc)),
                        "auprc_bootstrap_low_high": json.dumps(bootstrap_auprc(y, oof, groups)),
                        "prediction_target": "evidence-observation status, not therapeutic efficacy",
                    }
                    output_rows.append(result)
                    print(json.dumps({k: result[k] for k in ("task", "group_mode", "feature_set", "model", "oof_auprc", "oof_auroc")}), flush=True)

    output = project_path("reports", "tables", "f3_evidence_status_model_results.csv")
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    disease_rows = [row for row in output_rows if row["group_mode"] == "disease" and row["feature_set"] == "all_features"]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.8))
    for axis, (positive, _) in zip(axes, TASKS):
        task_rows = [row for row in disease_rows if row["task"] == f"{positive}_vs_U"]
        labels = [row["model"] for row in task_rows]
        values = [row["oof_auprc"] for row in task_rows]
        prevalence = task_rows[0]["positive_prevalence"]
        axis.bar(np.arange(2), values, color=["#3B6E8F", "#6B9E78"])
        axis.axhline(prevalence, color="black", linestyle="--", linewidth=1, label="Prevalence")
        axis.set_xticks(np.arange(2), labels)
        axis.set_ylim(0, 1)
        axis.set_title(f"{positive} vs U")
        axis.set_ylabel("Out-of-fold AUPRC")
    fig.suptitle("Predicting evidence-observation status from selection-bias features")
    fig.text(0.01, -0.04, "Disease-grouped five-fold evaluation. This is not drug-efficacy prediction.", fontsize=8)
    base = project_path("reports", "figures", "f3_evidence_status_auprc")
    plt.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "result_rows": len(output_rows), "bootstrap_iterations": BOOTSTRAPS,
        "features": FEATURE_SETS, "evaluation": "Five-fold stratified group CV by disease, with drug-group sensitivity",
        "scope": "Prediction of evidence-observation status only; not drug efficacy, clinical utility, or biological association",
    }
    report_path = project_path("reports", "f3_evidence_status_models.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

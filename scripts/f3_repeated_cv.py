import csv
import json
import sys
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

from trialevidence.paths import project_path

SEEDS = [20261110 + index for index in range(10)]
FEATURES = ["drug_auxiliary_degree", "disease_auxiliary_degree", "drug_observed_trials", "disease_observed_trials"]


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    variants = {
        "original": read_csv(project_path("data", "interim", "selection_bias_pair_features.csv")),
        "leave_label_source_out": read_csv(project_path("data", "interim", "selection_bias_pair_features_leave_source_out.csv")),
    }
    result_rows = []
    for variant, rows in variants.items():
        for task in ("P", "N_rule", "N_review"):
            subset = [row for row in rows if row["evidence_set"] in {task, "U"}]
            y = np.array([row["evidence_set"] == task for row in subset], dtype=int)
            groups = np.array([row["disease_node_index"] for row in subset])
            x = np.log1p(np.array([[float(row[feature]) for feature in FEATURES] for row in subset]))
            for seed in SEEDS:
                oof = np.zeros(len(y))
                for train, test in StratifiedGroupKFold(5, shuffle=True, random_state=seed).split(x, y, groups):
                    model = make_pipeline(StandardScaler(), LogisticRegression(class_weight="balanced", max_iter=2000, random_state=seed))
                    model.fit(x[train], y[train])
                    oof[test] = model.predict_proba(x[test])[:, 1]
                result_rows.append({
                    "variant": variant, "task": f"{task}_vs_U", "cv_seed": seed,
                    "n": len(y), "positive_prevalence": float(y.mean()),
                    "oof_auprc": average_precision_score(y, oof), "oof_auroc": roc_auc_score(y, oof),
                })
            print(json.dumps({"variant": variant, "task": task, "completed_seeds": len(SEEDS)}), flush=True)

    output = project_path("reports", "tables", "f3_repeated_cv_results.csv")
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(result_rows[0]))
        writer.writeheader()
        writer.writerows(result_rows)
    summary_rows = []
    for variant in variants:
        for task in ("P", "N_rule", "N_review"):
            values = np.array([row["oof_auprc"] for row in result_rows if row["variant"] == variant and row["task"] == f"{task}_vs_U"])
            summary_rows.append({
                "variant": variant, "task": f"{task}_vs_U", "repeats": len(values),
                "auprc_mean": float(values.mean()), "auprc_sd": float(values.std(ddof=1)),
                "auprc_min": float(values.min()), "auprc_max": float(values.max()),
            })
    with project_path("reports", "tables", "f3_repeated_cv_summary.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)

    fig, axes = plt.subplots(1, 3, figsize=(11, 4.6))
    for axis, task in zip(axes, ("P", "N_rule", "N_review")):
        data = [[row["oof_auprc"] for row in result_rows if row["variant"] == variant and row["task"] == f"{task}_vs_U"] for variant in variants]
        axis.boxplot(data, tick_labels=["Original", "Leave-source-out"])
        axis.set_title(f"{task} vs U")
        axis.set_ylabel("Out-of-fold AUPRC")
        axis.tick_params(axis="x", rotation=20)
    fig.suptitle("Evidence-status AUPRC across 10 retrained group-CV repetitions")
    fig.text(0.01, -0.04, "Distributions reflect retraining and split variability; U and evidence sources remain fixed.", fontsize=8)
    base = project_path("reports", "figures", "f3_repeated_cv")
    plt.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(), "seeds": SEEDS,
        "result_rows": len(result_rows),
        "scope": "Repeated disease-grouped retraining distributions, not confidence intervals or independent replication",
    }
    with project_path("reports", "f3_repeated_cv.json").open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

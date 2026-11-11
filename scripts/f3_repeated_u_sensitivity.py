import csv
import json
import random
import sys
from collections import Counter
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
from f3_selection_bias import LABEL_RELATIONS
from f3_evidence_status_models import FEATURE_SETS

SEED = 20261110
SEEDS = [SEED + index for index in range(10)]


def read_csv(name: str) -> list[dict]:
    with project_path("data", "interim", name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    auxiliary_degree = Counter()
    with project_path("data", "raw", "primekg", "edges.csv").open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            if row["display_relation"].casefold() not in LABEL_RELATIONS:
                auxiliary_degree[int(row["x_index"])] += 1
                auxiliary_degree[int(row["y_index"])] += 1
    trial_counts = {
        (row["node_type"], int(row["node_index"])): int(row["observed_trial_count"])
        for row in read_csv("entity_observed_trial_counts.csv")
    }
    positive = {(int(row["drug_node_index"]), int(row["disease_node_index"])) for row in read_csv("chembl_primekg_positive_candidates.csv") if row["mapping_status"] == "both_mapped_unique"}
    negative = {(int(row["drug_node_index"]), int(row["disease_node_index"])) for row in read_csv("n_strict_rule_based_pairs.csv")}
    review = {(int(row["drug_node_index"]), int(row["disease_node_index"])) for row in read_csv("completed_mapped_n_review_candidates.csv") if row["mapping_status"] == "both_mapped_unique"}
    negative -= positive
    review -= positive | negative
    occupied = positive | negative | review

    with project_path("data", "raw", "primekg", "nodes.tab").open(encoding="utf-8", newline="") as stream:
        sample = stream.read(4096)
        stream.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        nodes = list(csv.DictReader(stream, dialect=dialect))
    drugs = sorted(int(row["node_index"]) for row in nodes if row["node_type"] == "drug")
    diseases = sorted(int(row["node_index"]) for row in nodes if row["node_type"] == "disease")

    output_rows = []
    for seed in SEEDS:
        rng = random.Random(seed)
        unknown = set()
        while len(unknown) < 10000:
            pair = (rng.choice(drugs), rng.choice(diseases))
            if pair not in occupied:
                unknown.add(pair)
        for positive_name, positive_pairs in (("P", positive), ("N_rule", negative), ("N_review", review)):
            pairs = sorted(positive_pairs | unknown)
            y = np.array([pair in positive_pairs for pair in pairs], dtype=int)
            groups = np.array([pair[1] for pair in pairs])
            for feature_name, features in FEATURE_SETS.items():
                x = np.log1p(np.array([[
                    (auxiliary_degree[pair[0]] if feature == "drug_auxiliary_degree" else
                     auxiliary_degree[pair[1]] if feature == "disease_auxiliary_degree" else
                     trial_counts.get(("drug", pair[0]), 0) if feature == "drug_observed_trials" else
                     trial_counts.get(("disease", pair[1]), 0))
                    for feature in features
                ] for pair in pairs]))
                oof = np.zeros(len(y))
                for train, test in StratifiedGroupKFold(5, shuffle=True, random_state=SEED).split(x, y, groups):
                    model = make_pipeline(StandardScaler(), LogisticRegression(class_weight="balanced", max_iter=2000, random_state=SEED))
                    model.fit(x[train], y[train])
                    oof[test] = model.predict_proba(x[test])[:, 1]
                result = {
                    "task": f"{positive_name}_vs_U", "u_seed": seed, "feature_set": feature_name,
                    "n": len(y), "positive_prevalence": float(y.mean()),
                    "oof_auprc": average_precision_score(y, oof), "oof_auroc": roc_auc_score(y, oof),
                }
                output_rows.append(result)
        print(json.dumps({"completed_u_seed": seed}), flush=True)

    output = project_path("reports", "tables", "f3_repeated_u_results.csv")
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    summary_rows = []
    for task in ("P_vs_U", "N_rule_vs_U", "N_review_vs_U"):
        for feature_name in FEATURE_SETS:
            selected = [row for row in output_rows if row["task"] == task and row["feature_set"] == feature_name]
            values = np.array([row["oof_auprc"] for row in selected])
            summary_rows.append({
                "task": task, "feature_set": feature_name, "u_samples": len(values),
                "auprc_mean": float(values.mean()), "auprc_sd": float(values.std(ddof=1)),
                "auprc_min": float(values.min()), "auprc_max": float(values.max()),
                "auprc_q2_5": float(np.percentile(values, 2.5)), "auprc_q97_5": float(np.percentile(values, 97.5)),
            })
    with project_path("reports", "tables", "f3_repeated_u_summary.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)

    fig, axes = plt.subplots(1, 3, figsize=(11, 4.6))
    for axis, task in zip(axes, ("P_vs_U", "N_rule_vs_U", "N_review_vs_U")):
        data = [[row["oof_auprc"] for row in output_rows if row["task"] == task and row["feature_set"] == feature] for feature in FEATURE_SETS]
        axis.boxplot(data, tick_labels=list(FEATURE_SETS))
        axis.set_title(task)
        axis.set_ylabel("Out-of-fold AUPRC")
        axis.tick_params(axis="x", rotation=20)
    fig.suptitle("Evidence-status AUPRC across 10 random U draws")
    fig.text(0.01, -0.04, "Disease-grouped logistic regression. U is unknown, not a biological negative class.", fontsize=8)
    base = project_path("reports", "figures", "f3_repeated_u_sensitivity")
    plt.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(), "u_seeds": SEEDS,
        "result_rows": len(output_rows), "scope": "Sensitivity to U sampling; not a solution to same-source label/feature construction",
    }
    with project_path("reports", "f3_repeated_u_sensitivity.json").open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

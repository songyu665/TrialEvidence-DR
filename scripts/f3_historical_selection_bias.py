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

SEED = 20261110
FEATURE_SETS = {
    "degrees_only": ["drug_auxiliary_degree", "disease_auxiliary_degree"],
    "trials_only": ["drug_observed_trials", "disease_observed_trials"],
    "all_features": ["drug_auxiliary_degree", "disease_auxiliary_degree", "drug_observed_trials", "disease_observed_trials"],
}


def read_csv(name: str) -> list[dict]:
    with project_path("data", "interim", name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    from scipy.stats import rankdata
    ranks = rankdata(np.concatenate([a, b]))
    auc = (ranks[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b))
    return float(2 * auc - 1)


def main() -> None:
    auxiliary_degree = Counter()
    with project_path("data", "raw", "primekg", "edges.csv").open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            if row["display_relation"].casefold() not in LABEL_RELATIONS:
                auxiliary_degree[int(row["x_index"])] += 1
                auxiliary_degree[int(row["y_index"])] += 1
    trial_counts = {
        (row["node_type"], int(row["node_index"])): int(row["observed_trial_count"])
        for row in read_csv("entity_observed_trial_counts_aact20231201.csv")
    }
    positive = {(int(row["drug_node_index"]), int(row["disease_node_index"])) for row in read_csv("historical_chembl33_positive_pairs.csv")}
    negative = {(int(row["drug_node_index"]), int(row["disease_node_index"])) for row in read_csv("historical_20231201_n_rule_pairs.csv")}
    negative -= positive

    with project_path("data", "raw", "primekg", "nodes.tab").open(encoding="utf-8", newline="") as stream:
        sample = stream.read(4096)
        stream.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        nodes = list(csv.DictReader(stream, dialect=dialect))
    drugs = sorted(int(row["node_index"]) for row in nodes if row["node_type"] == "drug")
    diseases = sorted(int(row["node_index"]) for row in nodes if row["node_type"] == "disease")
    rng = random.Random(SEED)
    unknown = set()
    while len(unknown) < 10000:
        pair = (rng.choice(drugs), rng.choice(diseases))
        if pair not in positive | negative:
            unknown.add(pair)

    sets = {"P": positive, "N_rule": negative, "U": unknown}
    pair_rows = []
    for name, pairs in sets.items():
        for drug, disease in sorted(pairs):
            pair_rows.append({
                "evidence_set": name, "drug_node_index": drug, "disease_node_index": disease,
                "drug_auxiliary_degree": auxiliary_degree[drug], "disease_auxiliary_degree": auxiliary_degree[disease],
                "drug_observed_trials": trial_counts.get(("drug", drug), 0),
                "disease_observed_trials": trial_counts.get(("disease", disease), 0),
            })
    with project_path("data", "interim", "historical_selection_bias_pair_features.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(pair_rows[0]))
        writer.writeheader()
        writer.writerows(pair_rows)

    summary_rows = []
    model_rows = []
    for set_name in ("P", "N_rule"):
        a_rows = [row for row in pair_rows if row["evidence_set"] == set_name]
        b_rows = [row for row in pair_rows if row["evidence_set"] == "U"]
        for feature in FEATURE_SETS["all_features"]:
            a = np.array([row[feature] for row in a_rows])
            b = np.array([row[feature] for row in b_rows])
            summary_rows.append({
                "comparison": f"{set_name}_vs_U", "feature": feature,
                "median_a": float(np.median(a)), "median_b": float(np.median(b)),
                "cliffs_delta": cliffs_delta(a, b),
            })
        y = np.array([row["evidence_set"] == set_name for row in a_rows + b_rows], dtype=int)
        groups = np.array([row["disease_node_index"] for row in a_rows + b_rows])
        for feature_name, features in FEATURE_SETS.items():
            x = np.log1p(np.array([[row[feature] for feature in features] for row in a_rows + b_rows]))
            oof = np.zeros(len(y))
            for train, test in StratifiedGroupKFold(5, shuffle=True, random_state=SEED).split(x, y, groups):
                model = make_pipeline(StandardScaler(), LogisticRegression(class_weight="balanced", max_iter=2000, random_state=SEED))
                model.fit(x[train], y[train])
                oof[test] = model.predict_proba(x[test])[:, 1]
            model_rows.append({
                "task": f"{set_name}_vs_U", "feature_set": feature_name,
                "n": len(y), "positive_prevalence": float(y.mean()),
                "oof_auprc": average_precision_score(y, oof), "oof_auroc": roc_auc_score(y, oof),
                "prediction_target": "historical evidence-observation status, not efficacy",
            })

    for name, rows in (("historical_selection_bias_summary.csv", summary_rows), ("historical_evidence_status_models.csv", model_rows)):
        with project_path("reports", "tables", name).open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))
    for axis, task in zip(axes, ("P_vs_U", "N_rule_vs_U")):
        selected = [row for row in model_rows if row["task"] == task]
        labels = [row["feature_set"] for row in selected]
        values = [row["oof_auprc"] for row in selected]
        axis.bar(np.arange(3), values, color="#3B6E8F")
        axis.axhline(selected[0]["positive_prevalence"], color="black", linestyle="--", linewidth=1)
        axis.set_xticks(np.arange(3), labels, rotation=20, ha="right")
        axis.set_ylim(0, 1)
        axis.set_title(task)
        axis.set_ylabel("Historical out-of-fold AUPRC")
    fig.suptitle("Historical learnability of evidence-observation status")
    fig.text(0.01, -0.04, "PrimeKG 2022 + AACT 2023-12-01. Disease-grouped logistic regression; not efficacy prediction.", fontsize=8)
    base = project_path("reports", "figures", "f3_historical_evidence_status_auprc")
    plt.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "set_sizes": {name: len(pairs) for name, pairs in sets.items()},
        "historical_sources": "PrimeKG 2022, AACT 2023-12-01, ChEMBL 33",
        "evaluation": "Disease-grouped five-fold stratified logistic regression",
        "scope": "Historical sensitivity analysis, not independent external validation or efficacy prediction",
    }
    with project_path("reports", "f3_historical_selection_bias.json").open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps({"summary_rows": summary_rows, "model_rows": model_rows, "report": report}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

import csv
import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path

SEED = 20261110
BOOTSTRAPS = 500
FEATURES = [
    "drug_auxiliary_degree", "disease_auxiliary_degree",
    "drug_observed_trials", "disease_observed_trials",
]
COMPARISONS = [("P", "U"), ("N_rule", "U"), ("N_review", "U")]


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    combined = np.concatenate([a, b])
    ranks = rankdata(combined)
    auc = (ranks[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b))
    return float(2 * auc - 1)


def main() -> None:
    with project_path("data", "interim", "selection_bias_pair_features.csv").open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        for feature in FEATURES:
            row[feature] = float(row[feature])
    by_set = {name: [row for row in rows if row["evidence_set"] == name] for name in ("P", "N_rule", "N_review", "U")}
    rng = np.random.default_rng(SEED)
    output_rows = []
    for set_a, set_b in COMPARISONS:
        for feature in FEATURES:
            a = np.array([row[feature] for row in by_set[set_a]])
            b = np.array([row[feature] for row in by_set[set_b]])
            point = cliffs_delta(a, b)
            result = {
                "comparison": f"{set_a}_vs_{set_b}", "feature": feature,
                "n_a": len(a), "n_b": len(b),
                "median_a": float(np.median(a)), "median_b": float(np.median(b)),
                "median_difference": float(np.median(a) - np.median(b)),
                "cliffs_delta": point,
            }
            for cluster in ("disease_node_index", "drug_node_index"):
                a_groups = defaultdict(list)
                b_groups = defaultdict(list)
                for row in by_set[set_a]:
                    a_groups[row[cluster]].append(row[feature])
                for row in by_set[set_b]:
                    b_groups[row[cluster]].append(row[feature])
                a_group_list = list(a_groups)
                b_group_list = list(b_groups)
                a_values = [np.array(a_groups[group]) for group in a_group_list]
                b_values = [np.array(b_groups[group]) for group in b_group_list]
                differences = []
                for _ in range(BOOTSTRAPS):
                    a_indices = rng.integers(0, len(a_group_list), len(a_group_list))
                    b_indices = rng.integers(0, len(b_group_list), len(b_group_list))
                    a_sample = np.concatenate([a_values[index] for index in a_indices])
                    b_sample = np.concatenate([b_values[index] for index in b_indices])
                    differences.append(np.median(a_sample) - np.median(b_sample))
                result[f"{cluster}_bootstrap_low"] = float(np.percentile(differences, 2.5))
                result[f"{cluster}_bootstrap_high"] = float(np.percentile(differences, 97.5))
            output_rows.append(result)

    output = project_path("reports", "tables", "f3_selection_bias_effects.csv")
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    fig, axis = plt.subplots(figsize=(8, 8))
    ordered = [(comparison, feature) for comparison in COMPARISONS for feature in FEATURES]
    positions = np.arange(len(ordered))
    labels = []
    for position, (comparison, feature) in zip(positions, ordered):
        row = next(item for item in output_rows if item["comparison"] == f"{comparison[0]}_vs_{comparison[1]}" and item["feature"] == feature)
        labels.append(f"{comparison[0]} vs U | {feature.replace('_', ' ')}")
        axis.scatter(row["cliffs_delta"], position, color="#3B6E8F")
    axis.axvline(0, color="black", linewidth=0.8)
    axis.set_yticks(positions, labels)
    axis.set_xlabel("Cliff's delta")
    axis.invert_yaxis()
    fig.suptitle("Effect sizes for evidence-set selection differences")
    fig.text(0.01, -0.03, "Positive values indicate larger values in P, N_rule or N_review than U. Cluster intervals are in the CSV.", fontsize=8)
    base = project_path("reports", "figures", "f3_selection_bias_effect_sizes")
    plt.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "comparisons": [f"{a}_vs_{b}" for a, b in COMPARISONS],
        "features": FEATURES,
        "cliffs_delta_scale": "-1 to 1; positive means higher values in the evidence set than U",
        "bootstrap_iterations": BOOTSTRAPS,
        "bootstrap_unit": "Separate drug-entity and disease-entity cluster percentile intervals for median differences",
        "seed": SEED,
        "scope": "Exploratory selection-bias effect sizes; no causal inference and no claim of independent biological replication",
    }
    report_path = project_path("reports", "f3_selection_bias_effects.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps({"rows": len(output_rows), "report": report}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

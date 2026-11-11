import csv
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from trialevidence.paths import project_path
from f3_selection_bias import LABEL_RELATIONS
from f3_future_event_baseline import tied_ranks


def read_csv(name: str) -> list[dict]:
    with project_path("data", "interim", name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    trial_rows = read_csv("entity_observed_trial_counts_aact20231201.csv")
    trial_counts = {(row["node_type"], row["node_index"]): int(row["observed_trial_count"]) for row in trial_rows}
    drugs = {row["node_index"] for row in trial_rows if row["node_type"] == "drug"}
    diseases = {row["node_index"] for row in trial_rows if row["node_type"] == "disease"}
    for name in ("historical_chembl33_positive_pairs.csv", "historical_20231201_n_rule_pairs.csv"):
        for row in read_csv(name):
            drugs.add(row["drug_node_index"])
            diseases.add(row["disease_node_index"])
    candidate_pairs = {(drug, disease) for drug in drugs for disease in diseases}

    auxiliary_degree = Counter()
    with project_path("data", "raw", "primekg", "edges.csv").open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            if row["display_relation"].casefold() not in LABEL_RELATIONS:
                auxiliary_degree[int(row["x_index"])] += 1
                auxiliary_degree[int(row["y_index"])] += 1
    score_maps = {
        "drug_observed_trials": {pair: trial_counts.get(("drug", pair[0]), 0) for pair in candidate_pairs},
        "disease_observed_trials": {pair: trial_counts.get(("disease", pair[1]), 0) for pair in candidate_pairs},
        "total_observed_trials": {pair: trial_counts.get(("drug", pair[0]), 0) + trial_counts.get(("disease", pair[1]), 0) for pair in candidate_pairs},
        "total_auxiliary_degree": {pair: auxiliary_degree[int(pair[0])] + auxiliary_degree[int(pair[1])] for pair in candidate_pairs},
    }
    event_pairs = {(row["drug_node_index"], row["disease_node_index"]) for row in read_csv("post_cutoff_trial_event_candidates.csv")}
    evaluable = event_pairs & candidate_pairs
    result_rows = []
    for name, scores in score_maps.items():
        ranks = tied_ranks(candidate_pairs, scores)
        event_ranks = np.array([ranks[pair] for pair in evaluable])
        result_rows.append({
            "score": name, "candidate_pairs": len(candidate_pairs),
            "event_pairs_total": len(event_pairs), "event_pairs_evaluable_in_candidate_set": len(evaluable),
            "hits_at_10": float(np.mean(event_ranks <= 10)), "hits_at_100": float(np.mean(event_ranks <= 100)),
            "hits_at_1000": float(np.mean(event_ranks <= 1000)), "mrr": float(np.mean(1 / event_ranks)),
            "median_event_rank": float(np.median(event_ranks)),
        })
    with project_path("reports", "tables", "f3_future_event_expanded_baseline.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(result_rows[0]))
        writer.writeheader()
        writer.writerows(result_rows)

    with project_path("reports", "tables", "f3_future_event_baseline.csv").open(encoding="utf-8", newline="") as stream:
        original_rows = list(csv.DictReader(stream))
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))
    for axis, score in zip(axes, ("total_observed_trials", "total_auxiliary_degree")):
        old = next(row for row in original_rows if row["score"] == score)
        new = next(row for row in result_rows if row["score"] == score)
        axis.bar([0, 1], [float(old["hits_at_100"]), new["hits_at_100"]], color=["#B8B8B8", "#3B6E8F"])
        axis.set_xticks([0, 1], ["Observed pairs", "Entity Cartesian set"])
        axis.set_ylabel("Hits@100")
        axis.set_title(score.replace("_", " "))
    fig.suptitle("Effect of candidate-frame expansion on later-event ranking")
    fig.text(0.01, -0.04, "Same events and pre-cutoff features. Candidate definition sensitivity, not efficacy prediction.", fontsize=8)
    base = project_path("reports", "figures", "f3_future_event_candidate_frame")
    plt.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "historical_drugs": len(drugs), "historical_diseases": len(diseases),
        "candidate_pairs": len(candidate_pairs), "event_pairs_total": len(event_pairs),
        "event_pairs_evaluable": len(evaluable), "event_pairs_outside_candidate_set": len(event_pairs - candidate_pairs),
        "candidate_definition": "Cartesian product of drugs and diseases observed in historical AACT or historical P/N",
        "interpretation": "Candidate-frame sensitivity analysis, not therapeutic efficacy or complete prospective validation",
    }
    with project_path("reports", "f3_future_event_candidate_frame.json").open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps({"report": report, "results": result_rows}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

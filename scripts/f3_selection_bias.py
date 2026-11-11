import csv
import hashlib
import json
import random
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path

SEED = 20261110
U_SIZE = 10000
LABEL_RELATIONS = {"indication", "contraindication", "off-label use", "off_label_use"}


def read_csv(name: str) -> list[dict]:
    with project_path("data", "interim", name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    edge_path = project_path("data", "raw", "primekg", "edges.csv")
    checksum_path = project_path("data", "raw", "primekg", "edges.csv.sha256")
    if not edge_path.exists() or not checksum_path.exists():
        raise FileNotFoundError("PrimeKG edges.csv has not completed verified download")
    metadata = json.loads(project_path("data", "raw", "primekg", "dataset_metadata.json").read_text())
    files = metadata["data"]["latestVersion"]["files"]
    expected_md5 = next(item["dataFile"]["checksum"]["value"] for item in files if item["dataFile"]["filename"] == "edges.csv")
    hasher = hashlib.md5()
    with edge_path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            hasher.update(chunk)
    if hasher.hexdigest() != expected_md5:
        raise ValueError("PrimeKG edges.csv does not match Dataverse original-file MD5")
    total_degree = Counter()
    auxiliary_degree = Counter()
    relations = Counter()
    with edge_path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            x = int(row["x_index"])
            y = int(row["y_index"])
            total_degree[x] += 1
            total_degree[y] += 1
            relation = row["display_relation"].casefold()
            relations[relation] += 1
            if relation not in LABEL_RELATIONS:
                auxiliary_degree[x] += 1
                auxiliary_degree[y] += 1

    trial_counts = {
        (row["node_type"], int(row["node_index"])): int(row["observed_trial_count"])
        for row in read_csv("entity_observed_trial_counts.csv")
    }
    positive = {
        (int(row["drug_node_index"]), int(row["disease_node_index"]))
        for row in read_csv("chembl_primekg_positive_candidates.csv")
        if row["mapping_status"] == "both_mapped_unique"
    }
    negative = {
        (int(row["drug_node_index"]), int(row["disease_node_index"]))
        for row in read_csv("n_strict_rule_based_pairs.csv")
    }
    review = {
        (int(row["drug_node_index"]), int(row["disease_node_index"]))
        for row in read_csv("completed_mapped_n_review_candidates.csv")
        if row["mapping_status"] == "both_mapped_unique"
    }
    negative -= positive
    review -= positive | negative

    with project_path("data", "raw", "primekg", "nodes.tab").open(encoding="utf-8", newline="") as stream:
        sample = stream.read(4096)
        stream.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        nodes = list(csv.DictReader(stream, dialect=dialect))
    drugs = sorted(int(row["node_index"]) for row in nodes if row["node_type"] == "drug")
    diseases = sorted(int(row["node_index"]) for row in nodes if row["node_type"] == "disease")
    occupied = positive | negative | review
    rng = random.Random(SEED)
    unknown = set()
    while len(unknown) < U_SIZE:
        pair = (rng.choice(drugs), rng.choice(diseases))
        if pair not in occupied:
            unknown.add(pair)

    sets = {"P": positive, "N_rule": negative, "N_review": review, "U": unknown}
    pair_rows = []
    for set_name, pairs in sets.items():
        for drug, disease in sorted(pairs):
            pair_rows.append({
                "evidence_set": set_name, "drug_node_index": drug, "disease_node_index": disease,
                "drug_total_degree": total_degree[drug], "disease_total_degree": total_degree[disease],
                "drug_auxiliary_degree": auxiliary_degree[drug], "disease_auxiliary_degree": auxiliary_degree[disease],
                "drug_observed_trials": trial_counts.get(("drug", drug), 0),
                "disease_observed_trials": trial_counts.get(("disease", disease), 0),
            })
    pair_output = project_path("data", "interim", "selection_bias_pair_features.csv")
    with pair_output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(pair_rows[0]))
        writer.writeheader()
        writer.writerows(pair_rows)

    features = [
        "drug_total_degree", "disease_total_degree", "drug_auxiliary_degree", "disease_auxiliary_degree",
        "drug_observed_trials", "disease_observed_trials",
    ]
    summary_rows = []
    for set_name in sets:
        rows = [row for row in pair_rows if row["evidence_set"] == set_name]
        for feature in features:
            values = [row[feature] for row in rows]
            summary_rows.append({
                "evidence_set": set_name, "feature": feature, "n": len(values),
                "median": float(np.median(values)),
                "q1": float(np.percentile(values, 25)), "q3": float(np.percentile(values, 75)),
                "mean": float(np.mean(values)), "zero_fraction": float(np.mean(np.array(values) == 0)),
            })
    summary_path = project_path("reports", "tables", "f3_selection_bias_summary.csv")
    with summary_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)

    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    for axis, feature in zip(axes.flat, features):
        data = [np.log1p([row[feature] for row in pair_rows if row["evidence_set"] == set_name]) for set_name in sets]
        axis.boxplot(data, tick_labels=list(sets), showfliers=False)
        axis.set_title(feature.replace("_", " "))
        axis.set_ylabel("log(1 + value)")
    fig.suptitle("Selection-bias features across mutually exclusive evidence sets")
    fig.text(0.01, -0.02, "U is unknown, not negative. Trial counts are limited to downloaded query collections and mapping coverage.", fontsize=8)
    figure_base = project_path("reports", "figures", "f3_selection_bias_features")
    plt.savefig(figure_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(figure_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "set_sizes": {name: len(pairs) for name, pairs in sets.items()},
        "u_sampling_seed": SEED,
        "u_sampling_population": "PrimeKG drug-disease pairs excluding P, N_rule and N_review",
        "label_edges_excluded_from_auxiliary_degree": sorted(LABEL_RELATIONS),
        "observed_relations": dict(relations),
        "analysis_scope": "Descriptive distributions only; no significance test or causal claim",
        "limitations": [
            "Trial counts cover only the two downloaded ClinicalTrials.gov query collections",
            "Name-mapping failures can undercount observed research activity",
            "Small N_rule set supports exploratory description only",
        ],
    }
    report_path = project_path("reports", "f3_selection_bias.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

import bz2
import csv
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
from e0_map_chembl_primekg import normalize

LABEL_EDGES = {"treats", "palliates"}
SEED = 20261110


def read_csv(name: str) -> list[dict]:
    with project_path("data", "interim", name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    with bz2.open(project_path("data", "raw", "hetionet", "hetionet-v1.0.json.bz2"), "rt", encoding="utf-8") as stream:
        graph = json.load(stream)
    name_to_keys = {"Compound": defaultdict(list), "Disease": defaultdict(list)}
    for node in graph["nodes"]:
        if node["kind"] in name_to_keys:
            name_to_keys[node["kind"]][normalize(node["name"])].append((node["kind"], node["identifier"]))
    degree = Counter()
    auxiliary_degree = Counter()
    for edge in graph["edges"]:
        source = tuple(edge["source_id"])
        target = tuple(edge["target_id"])
        degree[source] += 1
        degree[target] += 1
        if edge["kind"].casefold() not in LABEL_EDGES:
            auxiliary_degree[source] += 1
            auxiliary_degree[target] += 1

    with project_path("data", "raw", "primekg", "nodes.tab").open(encoding="utf-8", newline="") as stream:
        sample = stream.read(4096)
        stream.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        primekg_nodes = list(csv.DictReader(stream, dialect=dialect))
    primekg = {int(row["node_index"]): row for row in primekg_nodes}

    pair_rows = []
    for row in read_csv("selection_bias_pair_features.csv"):
        drug_index = int(row["drug_node_index"])
        disease_index = int(row["disease_node_index"])
        drug_matches = name_to_keys["Compound"].get(normalize(primekg[drug_index]["node_name"]), [])
        disease_matches = name_to_keys["Disease"].get(normalize(primekg[disease_index]["node_name"]), [])
        if len(drug_matches) != 1 or len(disease_matches) != 1:
            continue
        drug = drug_matches[0]
        disease = disease_matches[0]
        pair_rows.append({
            "evidence_set": row["evidence_set"], "drug_node_index": drug_index, "disease_node_index": disease_index,
            "hetionet_compound_id": drug[1], "hetionet_disease_id": disease[1],
            "hetionet_drug_total_degree": degree[drug], "hetionet_disease_total_degree": degree[disease],
            "hetionet_drug_auxiliary_degree": auxiliary_degree[drug],
            "hetionet_disease_auxiliary_degree": auxiliary_degree[disease],
        })
    with project_path("data", "interim", "hetionet_pair_features.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(pair_rows[0]))
        writer.writeheader()
        writer.writerows(pair_rows)

    summary_rows = []
    model_rows = []
    for set_name in ("P", "N_rule", "N_review"):
        selected = [row for row in pair_rows if row["evidence_set"] == set_name]
        unknown = [row for row in pair_rows if row["evidence_set"] == "U"]
        for feature in ("hetionet_drug_auxiliary_degree", "hetionet_disease_auxiliary_degree"):
            a = np.array([row[feature] for row in selected])
            b = np.array([row[feature] for row in unknown])
            summary_rows.append({
                "comparison": f"{set_name}_vs_U", "feature": feature,
                "n_a": len(a), "n_b": len(b), "median_a": float(np.median(a)), "median_b": float(np.median(b)),
            })
        subset = selected + unknown
        y = np.array([row["evidence_set"] == set_name for row in subset], dtype=int)
        groups = np.array([row["disease_node_index"] for row in subset])
        x = np.log1p(np.array([[row["hetionet_drug_auxiliary_degree"], row["hetionet_disease_auxiliary_degree"]] for row in subset]))
        oof = np.zeros(len(y))
        for train, test in StratifiedGroupKFold(5, shuffle=True, random_state=SEED).split(x, y, groups):
            model = make_pipeline(StandardScaler(), LogisticRegression(class_weight="balanced", max_iter=2000, random_state=SEED))
            model.fit(x[train], y[train])
            oof[test] = model.predict_proba(x[test])[:, 1]
        model_rows.append({
            "task": f"{set_name}_vs_U", "n": len(y), "positive_prevalence": float(y.mean()),
            "oof_auprc": average_precision_score(y, oof), "oof_auroc": roc_auc_score(y, oof),
            "features": "Hetionet auxiliary degrees only",
            "interpretation": "Not interpretable as validation because only 12 U pairs were mapped",
        })

    for name, rows in (("f3_hetionet_summary.csv", summary_rows), ("f3_hetionet_models.csv", model_rows)):
        with project_path("reports", "tables", name).open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))
    for axis, feature in zip(axes, ("hetionet_drug_auxiliary_degree", "hetionet_disease_auxiliary_degree")):
        data = [np.log1p([row[feature] for row in pair_rows if row["evidence_set"] == name]) for name in ("P", "N_rule", "N_review", "U")]
        axis.boxplot(data, tick_labels=["P", "N_rule", "N_review", "U"], showfliers=False)
        axis.set_title(feature.replace("_", " "))
        axis.set_ylabel("log(1 + Hetionet degree)")
    fig.suptitle("Evidence-set degree distributions in Hetionet")
    fig.text(0.01, -0.04, "Exact-name-matched subset; treatment and palliation edges excluded. Shared labels, different graph.", fontsize=8)
    base = project_path("reports", "figures", "f3_hetionet_degrees")
    plt.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()

    original_counts = Counter(row["evidence_set"] for row in read_csv("selection_bias_pair_features.csv"))
    mapped_counts = Counter(row["evidence_set"] for row in pair_rows)
    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "hetionet_nodes": len(graph["nodes"]), "hetionet_edges": len(graph["edges"]),
        "original_pair_counts": dict(original_counts),
        "mapped_pair_counts": dict(mapped_counts),
        "mapped_pair_fractions": {name: mapped_counts[name] / original_counts[name] for name in original_counts},
        "mapping": "Unique exact normalized names between PrimeKG and Hetionet; no fuzzy matching",
        "label_edges_excluded": sorted(LABEL_EDGES),
        "validation_status": "not_informative_due_to_extreme_mapping_selection_only_12_U_pairs",
        "model_interpretation": "AUPRC values are not used as validation evidence because the mapped subset is highly selected",
        "scope": "Graph-structure sensitivity analysis over shared evidence labels, not independent label validation",
    }
    with project_path("reports", "f3_hetionet_validation.json").open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps({"report": report, "model_rows": model_rows}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

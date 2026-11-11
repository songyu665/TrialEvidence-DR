import csv
import json
import sys
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from trialevidence.paths import project_path
from e0_build_adjudication_queue import ACTIVE_TYPES, NONACTIVE
from e0_map_chembl_primekg import read_nodes, unique_matches
from f3_selection_bias import LABEL_RELATIONS

CUTOFF = "2023-12-01"
csv.field_size_limit(sys.maxsize)


def read_csv(name: str) -> list[dict]:
    with project_path("data", "interim", name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def rows_from_zip(z: zipfile.ZipFile, name: str):
    with z.open(name) as stream:
        wrapper = (line.decode("utf-8-sig" if index == 0 else "utf-8") for index, line in enumerate(stream))
        yield from csv.DictReader(wrapper, delimiter="|")


def tied_ranks(pairs: list, scores: dict) -> dict:
    ordered = sorted(pairs, key=lambda pair: (-scores[pair], pair))
    ranks = {}
    start = 0
    while start < len(ordered):
        end = start + 1
        score = scores[ordered[start]]
        while end < len(ordered) and scores[ordered[end]] == score:
            end += 1
        average = (start + 1 + end) / 2
        for pair in ordered[start:end]:
            ranks[pair] = average
        start = end
    return ranks


def main() -> None:
    known = {"drug": {}, "disease": {}}
    for row in read_csv("observed_trial_name_mapping.csv"):
        if row["node_index"]:
            known[row["node_type"]][row["raw_name"]] = row["node_index"]
    nodes = read_nodes()

    def map_name(raw: str, node_type: str):
        if raw in known[node_type]:
            return known[node_type][raw]
        matches = unique_matches({raw}, nodes[node_type])
        return matches[0]["node_index"] if len(matches) == 1 else None

    selected = {}
    archive = project_path("data", "raw", "aact", "20231201_aact.zip")
    with zipfile.ZipFile(archive) as z:
        for row in rows_from_zip(z, "studies.txt"):
            phase = row.get("phase", "")
            if row.get("study_type") != "Interventional" or not ("Phase 2" in phase or "Phase 3" in phase):
                continue
            category = None
            if row.get("overall_status") == "Terminated":
                category = "terminated"
            elif row.get("overall_status") == "Completed" and row.get("results_first_posted_date") and row["results_first_posted_date"] <= CUTOFF:
                category = "completed_results"
            if category:
                selected[row["nct_id"]] = {"interventions": [], "conditions": []}
        for row in rows_from_zip(z, "interventions.txt"):
            if row["nct_id"] in selected:
                selected[row["nct_id"]]["interventions"].append({"type": row.get("intervention_type", "").upper(), "name": row.get("name", "")})
        for row in rows_from_zip(z, "conditions.txt"):
            if row["nct_id"] in selected:
                selected[row["nct_id"]]["conditions"].append(row.get("name", ""))

    candidate_pairs = set()
    for trial in selected.values():
        drugs = set()
        diseases = set()
        for item in trial["interventions"]:
            if item["type"] not in ACTIVE_TYPES or not item["name"] or NONACTIVE.search(item["name"]):
                continue
            mapped = map_name(item["name"], "drug")
            if mapped:
                drugs.add(mapped)
        for name in trial["conditions"]:
            mapped = map_name(name, "disease")
            if mapped:
                diseases.add(mapped)
        candidate_pairs.update((drug, disease) for drug in drugs for disease in diseases)
    candidate_pairs |= {(row["drug_node_index"], row["disease_node_index"]) for row in read_csv("historical_chembl33_positive_pairs.csv")}
    candidate_pairs |= {(row["drug_node_index"], row["disease_node_index"]) for row in read_csv("historical_20231201_n_rule_pairs.csv")}

    auxiliary_degree = Counter()
    with project_path("data", "raw", "primekg", "edges.csv").open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            if row["display_relation"].casefold() not in LABEL_RELATIONS:
                auxiliary_degree[int(row["x_index"])] += 1
                auxiliary_degree[int(row["y_index"])] += 1
    trials = {
        (row["node_type"], row["node_index"]): int(row["observed_trial_count"])
        for row in read_csv("entity_observed_trial_counts_aact20231201.csv")
    }

    score_maps = {
        "drug_observed_trials": {pair: trials.get(("drug", pair[0]), 0) for pair in candidate_pairs},
        "disease_observed_trials": {pair: trials.get(("disease", pair[1]), 0) for pair in candidate_pairs},
        "total_observed_trials": {pair: trials.get(("drug", pair[0]), 0) + trials.get(("disease", pair[1]), 0) for pair in candidate_pairs},
        "total_auxiliary_degree": {pair: auxiliary_degree[int(pair[0])] + auxiliary_degree[int(pair[1])] for pair in candidate_pairs},
    }

    event_pairs = {(row["drug_node_index"], row["disease_node_index"]) for row in read_csv("post_cutoff_trial_event_candidates.csv")}
    evaluable = event_pairs & candidate_pairs
    result_rows = []
    for score_name, scores in score_maps.items():
        ranks = tied_ranks(candidate_pairs, scores)
        event_ranks = np.array([ranks[pair] for pair in evaluable])
        result_rows.append({
            "score": score_name, "candidate_pairs": len(candidate_pairs),
            "event_pairs_total": len(event_pairs), "event_pairs_evaluable_in_candidate_set": len(evaluable),
            "hits_at_10": float(np.mean(event_ranks <= 10)), "hits_at_100": float(np.mean(event_ranks <= 100)),
            "hits_at_1000": float(np.mean(event_ranks <= 1000)), "mrr": float(np.mean(1 / event_ranks)),
            "median_event_rank": float(np.median(event_ranks)),
        })
    output = project_path("reports", "tables", "f3_future_event_baseline.csv")
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(result_rows[0]))
        writer.writeheader()
        writer.writerows(result_rows)

    fig, ax = plt.subplots(figsize=(8, 4.8))
    x = np.arange(len(result_rows))
    ax.bar(x, [row["hits_at_100"] for row in result_rows], color="#3B6E8F")
    short_labels = {"drug_observed_trials": "Drug trials", "disease_observed_trials": "Disease trials", "total_observed_trials": "Total trials", "total_auxiliary_degree": "Auxiliary degree"}
    ax.set_xticks(x, [short_labels[row["score"]] for row in result_rows], rotation=15, ha="right")
    ax.set_ylabel("Hits@100 for post-cutoff evidence events")
    ax.set_title("Pre-cutoff graph and research-attention ranking baselines")
    fig.text(0.01, -0.04, "Events are later evidence observations, not clinical efficacy. Ties receive average ranks.", fontsize=8)
    base = project_path("reports", "figures", "f3_future_event_baseline")
    plt.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(), "candidate_pairs": len(candidate_pairs),
        "event_pairs_total": len(event_pairs), "event_pairs_evaluable": len(evaluable),
        "event_pairs_outside_candidate_set": len(event_pairs - candidate_pairs),
        "ranking": "Descending pre-cutoff feature; average ranks assigned to ties",
        "interpretation": "Prediction of later evidence observation, not therapeutic efficacy; no Precision@K because non-events are not negatives",
    }
    with project_path("reports", "f3_future_event_baseline.json").open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps({"report": report, "results": result_rows}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

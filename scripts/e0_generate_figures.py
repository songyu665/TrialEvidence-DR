import csv
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path


def report(name: str) -> dict:
    return json.loads(project_path("reports", name).read_text(encoding="utf-8"))


def write_csv(name: str, rows: list[dict]) -> None:
    path = project_path("reports", "tables", name)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_figure(name: str) -> None:
    path = project_path("reports", "figures", name)
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()


def main() -> None:
    plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans", "axes.spines.top": False, "axes.spines.right": False})
    integrity = report("raw_download_integrity.json")["datasets"]
    trials = report("e0_trial_sample_audit.json")
    negative = report("e0_aggregated_n_strict_candidates.json")
    review = report("e0_conservative_rule_review.json")
    gate = report("e0_rule_based_gate.json")
    ai = report("e0_ai_assisted_review_200_audited.json")["strata"]

    integrity_rows = []
    for name in ("clinicaltrials_terminated_phase23", "clinicaltrials_completed_results_phase23_with_comments", "chembl_drug_indication"):
        item = integrity[name]
        integrity_rows.append({"dataset": name, "records": item["records_verified"], "complete": item["download_complete"], "integrity": item["status"]})
    write_csv("e0_data_integrity.csv", integrity_rows)

    attrition = [
        ("Terminated phase II/III records", trials["records_audited"], "trial"),
        ("Explicit futility/lack-of-efficacy text", trials["candidate_buckets"]["N_strict_candidate_requires_three_layer_review"], "trial"),
        ("Unique mapped candidates", negative["unique_pair_candidates"], "pair"),
        ("Passed conservative trial-pair records", review["passed_trial_pair_records"], "trial-pair"),
        ("Unique rule-based candidates", review["passed_unique_pair_candidates"], "pair"),
        ("Nonoverlapping rule-based candidates", gate["nonoverlapping_rule_based_pairs"], "pair"),
    ]
    write_csv("e0_candidate_attrition.csv", [{"stage": stage, "count": count, "unit": unit} for stage, count, unit in attrition])
    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    y = np.arange(len(attrition))
    values = [item[1] for item in attrition]
    ax.barh(y, values, color="#3B6E8F")
    ax.set_yticks(y, [item[0] for item in attrition])
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("Count (log scale)")
    ax.set_title("Clinical-trial negative-evidence candidate attrition")
    for index, (_, value, unit) in enumerate(attrition):
        ax.text(value * 1.12, index, f"{value} {unit}s", va="center")
    ax.set_xlim(1, max(values) * 3)
    fig.text(0.01, -0.05, "Units change across stages; bars are not independent samples or model results.", fontsize=8)
    save_figure("e0_candidate_attrition")

    with project_path("data", "interim", "n_strict_pair_candidates_pending_review.csv").open(encoding="utf-8-sig", newline="") as stream:
        pairs = list(csv.DictReader(stream))
    disease_counts = Counter(row["disease_node_name"] for row in pairs)
    top = disease_counts.most_common(15)
    write_csv("e0_top_disease_candidates.csv", [{"disease": name, "unique_pair_candidates": count} for name, count in top])
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    names = [item[0] for item in reversed(top)]
    values = [item[1] for item in reversed(top)]
    ax.barh(np.arange(len(top)), values, color="#6B9E78")
    ax.set_yticks(np.arange(len(top)), names)
    ax.set_xlabel("Unique mapped pair candidates")
    ax.set_title("Disease coverage before conservative semantic review")
    for index, value in enumerate(values):
        ax.text(value + 0.2, index, str(value), va="center")
    ax.set_xlim(0, max(values) + 3)
    save_figure("e0_disease_coverage")

    decisions = ["supports_negative_candidate", "does_not_support", "insufficient"]
    labels = ["Support", "Do not support", "Insufficient"]
    colors = ["#2E7D5B", "#B45A4A", "#B8B8B8"]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    left = np.zeros(2)
    for decision, label, color in zip(decisions, labels, colors):
        values = np.array([ai["high"]["decisions"][decision], ai["medium"]["decisions"][decision]])
        ax.barh([0, 1], values, left=left, label=label, color=color)
        left += values
    ax.set_yticks([0, 1], ["High text priority\n(census, n=8)", "Medium statistical priority\n(random sample, n=192)"])
    ax.set_xlabel("Pairs")
    ax.set_title("AI-assisted screening outcomes by review stratum")
    ax.legend(loc="upper right")
    fig.text(0.01, -0.05, "AI-assisted screening only; strata are not pooled and clinical expert review remains pending.", fontsize=8)
    save_figure("e0_ai_review")

    input_names = [
        "raw_download_integrity.json", "e0_trial_sample_audit.json", "e0_aggregated_n_strict_candidates.json",
        "e0_conservative_rule_review.json", "e0_rule_based_gate.json", "e0_ai_assisted_review_200_audited.json",
    ]
    manifest = {
        "inputs": {name: hashlib.sha256(project_path("reports", name).read_bytes()).hexdigest() for name in input_names},
        "outputs": ["e0_candidate_attrition", "e0_disease_coverage", "e0_ai_review"],
        "scope": "Descriptive E0 audit outputs; no model performance or clinical efficacy claims",
    }
    manifest_path = project_path("reports", "figures", "e0_figure_manifest.json")
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

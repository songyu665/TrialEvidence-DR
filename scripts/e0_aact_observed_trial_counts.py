import csv
import json
import sys
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from trialevidence.paths import project_path
from e0_build_adjudication_queue import ACTIVE_TYPES, NONACTIVE
from e0_map_chembl_primekg import read_nodes, unique_matches

SNAPSHOT = "2023-12-01"
csv.field_size_limit(sys.maxsize)


def rows_from_zip(z: zipfile.ZipFile, name: str):
    with z.open(name) as stream:
        wrapper = (line.decode("utf-8-sig" if index == 0 else "utf-8") for index, line in enumerate(stream))
        yield from csv.DictReader(wrapper, delimiter="|")


def main() -> None:
    known = {"drug": {}, "disease": {}}
    with project_path("data", "interim", "observed_trial_name_mapping.csv").open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            if row["node_index"]:
                known[row["node_type"]][row["raw_name"]] = (row["node_index"], row["node_name"])
    nodes = read_nodes()
    node_names = {
        node_type: {row["node_index"]: row["node_name"] for matches in nodes[node_type].values() for row in matches}
        for node_type in ("drug", "disease")
    }

    def map_name(raw: str, node_type: str):
        if raw in known[node_type]:
            return known[node_type][raw]
        matches = unique_matches({raw}, nodes[node_type])
        return (matches[0]["node_index"], matches[0]["node_name"]) if len(matches) == 1 else (None, None)

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
            elif row.get("overall_status") == "Completed" and row.get("results_first_posted_date") and row["results_first_posted_date"] <= SNAPSHOT:
                category = "completed_results"
            if category:
                selected[row["nct_id"]] = {"category": category, "interventions": [], "conditions": []}
        for row in rows_from_zip(z, "interventions.txt"):
            if row["nct_id"] in selected:
                selected[row["nct_id"]]["interventions"].append({"type": row.get("intervention_type", "").upper(), "name": row.get("name", "")})
        for row in rows_from_zip(z, "conditions.txt"):
            if row["nct_id"] in selected:
                selected[row["nct_id"]]["conditions"].append(row.get("name", ""))

    counts = {"drug": {"terminated": Counter(), "completed_results": Counter()}, "disease": {"terminated": Counter(), "completed_results": Counter()}}
    coverage = Counter()
    for trial in selected.values():
        drug_nodes = set()
        disease_nodes = set()
        for item in trial["interventions"]:
            if item["type"] not in ACTIVE_TYPES or not item["name"] or NONACTIVE.search(item["name"]):
                continue
            node_id, _ = map_name(item["name"], "drug")
            if node_id:
                drug_nodes.add(node_id)
        for name in trial["conditions"]:
            node_id, _ = map_name(name, "disease")
            if node_id:
                disease_nodes.add(node_id)
        coverage[trial["category"]] += 1
        coverage[trial["category"] + "_with_mapped_drug"] += bool(drug_nodes)
        coverage[trial["category"] + "_with_mapped_disease"] += bool(disease_nodes)
        for node_id in drug_nodes:
            counts["drug"][trial["category"]][node_id] += 1
        for node_id in disease_nodes:
            counts["disease"][trial["category"]][node_id] += 1

    output_rows = []
    for node_type in ("drug", "disease"):
        identifiers = set(counts[node_type]["terminated"]) | set(counts[node_type]["completed_results"])
        for identifier in sorted(identifiers):
            terminated = counts[node_type]["terminated"][identifier]
            completed = counts[node_type]["completed_results"][identifier]
            node_name = node_names[node_type][identifier]
            output_rows.append({
                "node_type": node_type, "node_index": identifier, "node_name": node_name,
                "terminated_trial_count": terminated, "completed_with_results_trial_count": completed,
                "observed_trial_count": terminated + completed, "snapshot": SNAPSHOT,
            })
    output = project_path("data", "interim", "entity_observed_trial_counts_aact20231201.csv")
    with output.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    report = {
        "created_at": datetime.now(timezone.utc).isoformat(), "snapshot": SNAPSHOT,
        "coverage": dict(coverage), "entities_with_observed_trials": len(output_rows),
        "scope": "AACT interventional phase II/III terminated trials and completed trials with results posted by the snapshot",
        "mapping_note": "Frozen entity normalization was used without using current outcome labels",
    }
    report_path = project_path("reports", "e0_aact_observed_trial_counts.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

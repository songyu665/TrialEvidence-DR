import csv
import gzip
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from trialevidence.paths import project_path
from e0_build_adjudication_queue import ACTIVE_TYPES, NONACTIVE
from e0_map_chembl_primekg import normalize, read_nodes, unique_matches
from e0_map_negative_drug_synonyms import candidate_transformations

DATASETS = {
    "terminated": "terminated_phase23",
    "completed_results": "completed_results_phase23_with_comments",
}


def main() -> None:
    nodes = read_nodes()
    molecule_dir = project_path("data", "raw", "chembl", "molecule_dictionary")
    summary = json.loads((molecule_dir / "download_summary.json").read_text())
    molecules = []
    for batch in range(1, summary["next_batch"] + 1):
        with gzip.open(molecule_dir / f"batch_{batch:06d}.json.gz", "rt", encoding="utf-8") as stream:
            molecules.extend(json.load(stream)["molecules"])
    primekg_by_chembl = {}
    alias_to_chembl = defaultdict(set)
    for molecule in molecules:
        identifier = molecule["molecule_chembl_id"]
        names = {molecule.get("pref_name") or ""}
        for synonym in molecule.get("molecule_synonyms", []):
            for field in ("molecule_synonym", "synonyms"):
                if synonym.get(field):
                    names.add(synonym[field])
        names.discard("")
        for name in names:
            alias_to_chembl[normalize(name)].add(identifier)
        matches = unique_matches(names, nodes["drug"])
        if len(matches) == 1:
            primekg_by_chembl[identifier] = matches[0]

    drug_cache = {}
    disease_cache = {}

    def map_drug(raw: str):
        if raw in drug_cache:
            return drug_cache[raw]
        for candidate, rule in candidate_transformations(raw):
            direct = unique_matches({candidate}, nodes["drug"])
            chembl_ids = alias_to_chembl.get(normalize(candidate), set())
            indirect = [primekg_by_chembl[identifier] for identifier in chembl_ids if identifier in primekg_by_chembl]
            matches = {match["node_index"]: match for match in direct + indirect}
            if len(matches) == 1:
                match = next(iter(matches.values()))
                drug_cache[raw] = (match, rule)
                return drug_cache[raw]
        drug_cache[raw] = (None, "unmapped")
        return drug_cache[raw]

    def map_disease(raw: str):
        if raw not in disease_cache:
            matches = unique_matches({raw}, nodes["disease"])
            disease_cache[raw] = matches[0] if len(matches) == 1 else None
        return disease_cache[raw]

    counts = {"drug": {name: Counter() for name in DATASETS}, "disease": {name: Counter() for name in DATASETS}}
    coverage = {name: Counter() for name in DATASETS}
    for dataset, directory_name in DATASETS.items():
        directory = project_path("data", "raw", "clinicaltrials", directory_name)
        state = json.loads((directory / "download_checkpoint.json").read_text())
        if state["next_token"] is not None:
            raise RuntimeError(f"Dataset is incomplete: {dataset}")
        for page in range(1, state["pages_downloaded"] + 1):
            with gzip.open(directory / f"page_{page:06d}.json.gz", "rt", encoding="utf-8") as stream:
                studies = json.load(stream)["studies"]
            for study in studies:
                protocol = study.get("protocolSection", {})
                interventions = protocol.get("armsInterventionsModule", {}).get("interventions", [])
                conditions = protocol.get("conditionsModule", {}).get("conditions", [])
                drug_nodes = set()
                disease_nodes = set()
                for intervention in interventions:
                    name = (intervention.get("name") or "").strip()
                    if intervention.get("type") not in ACTIVE_TYPES or not name or NONACTIVE.search(name):
                        continue
                    drug, _ = map_drug(name)
                    if drug:
                        drug_nodes.add(drug["node_index"])
                for condition in conditions:
                    disease = map_disease(condition)
                    if disease:
                        disease_nodes.add(disease["node_index"])
                coverage[dataset]["trials"] += 1
                coverage[dataset]["with_mapped_drug"] += bool(drug_nodes)
                coverage[dataset]["with_mapped_disease"] += bool(disease_nodes)
                for node in drug_nodes:
                    counts["drug"][dataset][node] += 1
                for node in disease_nodes:
                    counts["disease"][dataset][node] += 1

    name_rows = []
    for raw, (node, rule) in sorted(drug_cache.items()):
        name_rows.append({"node_type": "drug", "raw_name": raw, "node_index": node["node_index"] if node else "", "node_name": node["node_name"] if node else "", "mapping_rule": rule})
    for raw, node in sorted(disease_cache.items()):
        name_rows.append({"node_type": "disease", "raw_name": raw, "node_index": node["node_index"] if node else "", "node_name": node["node_name"] if node else "", "mapping_rule": "exact_normalized_unique" if node else "unmapped"})
    name_path = project_path("data", "interim", "observed_trial_name_mapping.csv")
    with name_path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(name_rows[0]))
        writer.writeheader()
        writer.writerows(name_rows)

    output_rows = []
    for node_type in ("drug", "disease"):
        node_index = nodes[node_type]
        index_to_node = {row["node_index"]: row for rows in node_index.values() for row in rows}
        all_ids = set(counts[node_type]["terminated"]) | set(counts[node_type]["completed_results"])
        for identifier in sorted(all_ids):
            node = index_to_node[identifier]
            terminated = counts[node_type]["terminated"][identifier]
            completed = counts[node_type]["completed_results"][identifier]
            output_rows.append({
                "node_type": node_type, "node_index": identifier, "node_name": node["node_name"],
                "terminated_trial_count": terminated, "completed_with_results_trial_count": completed,
                "observed_trial_count": terminated + completed,
            })
    output = project_path("data", "interim", "entity_observed_trial_counts.csv")
    with output.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "coverage": {name: dict(values) for name, values in coverage.items()},
        "unique_raw_drug_names": len(drug_cache),
        "mapped_raw_drug_names": sum(node is not None for node, _ in drug_cache.values()),
        "unique_raw_disease_names": len(disease_cache),
        "mapped_raw_disease_names": sum(node is not None for node in disease_cache.values()),
        "entities_with_observed_trials": len(output_rows),
        "scope": "Counts cover only the terminated phase II/III and completed-with-results phase II/III query collections; they are not total ClinicalTrials.gov counts",
        "mapping_warning": "Unmapped combination regimens and ambiguous names reduce measured research activity; no fuzzy imputation is used",
    }
    report_path = project_path("reports", "e0_observed_trial_counts.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

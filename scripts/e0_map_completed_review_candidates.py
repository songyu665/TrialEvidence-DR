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
from e0_build_adjudication_queue import ACTIVE_TYPES, COMBINATION, NONACTIVE
from e0_map_chembl_primekg import normalize, read_nodes, unique_matches
from e0_map_negative_drug_synonyms import candidate_transformations

csv.field_size_limit(sys.maxsize)


def read_csv(name: str) -> list[dict]:
    with project_path("data", "interim", name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    analyses = [row for row in read_csv("completed_primary_analysis_review.csv") if row["review_priority"] in {"high_text_supported_N_review", "medium_statistical_N_review"} and row["tested_noninferiority"] != "True"]
    priorities = {}
    for row in analyses:
        priorities.setdefault(row["nct_id"], set()).add(row["review_priority"])
    trials = {row["nct_id"]: row for row in read_csv("completed_primary_outcome_review_queue.csv") if row["nct_id"] in priorities}

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

    def map_drug(raw: str):
        for candidate, rule in candidate_transformations(raw):
            direct = unique_matches({candidate}, nodes["drug"])
            chembl_ids = alias_to_chembl.get(normalize(candidate), set())
            indirect = [primekg_by_chembl[identifier] for identifier in chembl_ids if identifier in primekg_by_chembl]
            matches = {match["node_index"]: match for match in direct + indirect}
            if len(matches) == 1:
                match = next(iter(matches.values()))
                return match, rule
        return None, "unmapped"

    output_rows = []
    for nct_id, trial in sorted(trials.items()):
        interventions = json.loads(trial["interventions"])
        conditions = json.loads(trial["conditions"])
        for intervention in interventions:
            name = (intervention.get("name") or "").strip()
            if intervention.get("type") not in ACTIVE_TYPES or not name or NONACTIVE.search(name) or COMBINATION.search(name):
                continue
            drug, drug_rule = map_drug(name)
            for condition in conditions:
                disease_matches = unique_matches({condition}, nodes["disease"])
                disease = disease_matches[0] if len(disease_matches) == 1 else None
                output_rows.append({
                    "trial_id": nct_id,
                    "raw_intervention_name": name,
                    "raw_condition": condition,
                    "review_priorities": json.dumps(sorted(priorities[nct_id])),
                    "drug_node_index": drug["node_index"] if drug else "",
                    "drug_node_name": drug["node_name"] if drug else "",
                    "drug_mapping_rule": drug_rule,
                    "disease_node_index": disease["node_index"] if disease else "",
                    "disease_node_name": disease["node_name"] if disease else "",
                    "mapping_status": "both_mapped_unique" if drug and disease else "incomplete_or_ambiguous",
                    "label": "U_pending_review",
                    "evidence_url": trial["evidence_url"],
                })
    output = project_path("data", "interim", "completed_mapped_n_review_candidates.csv")
    temporary = output.with_suffix(".csv.tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    temporary.replace(output)

    mapped = [row for row in output_rows if row["mapping_status"] == "both_mapped_unique"]
    pairs = {(row["drug_node_index"], row["disease_node_index"]) for row in mapped}
    high_pairs = {(row["drug_node_index"], row["disease_node_index"]) for row in mapped if "high_text_supported_N_review" in json.loads(row["review_priorities"])}
    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "eligible_trials_with_statistical_or_text_review_priority": len(trials),
        "candidate_rows_after_nonactive_and_combination_exclusion": len(output_rows),
        "mapped_trial_pair_records": len(mapped),
        "mapped_unique_n_review_pairs": len(pairs),
        "mapped_unique_drugs": len({pair[0] for pair in pairs}),
        "mapped_unique_diseases": len({pair[1] for pair in pairs}),
        "high_priority_mapped_unique_pairs": len(high_pairs),
        "review_priority_trial_counts": dict(Counter(priority for values in priorities.values() for priority in values)),
        "warning": "These are N_review pairs for adjudication, not N_strict labels; statistical nonsignificance and text priority do not establish inefficacy",
    }
    report_path = project_path("reports", "e0_completed_n_review_mapping.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

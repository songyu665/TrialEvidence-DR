import csv
import gzip
import json
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from trialevidence.paths import project_path
from e0_map_chembl_primekg import normalize, read_nodes, unique_matches

SALT_SUFFIXES = ("hydrochloride", "dimeglumine", "phosphate", "mesylate", "citrate", "calcium", "sodium")
DOSE = re.compile(r"\b\d+(?:\.\d+)?\s*(?:mg|mcg|g|ml|milligrams?|micrograms?|grams?|milliliters?)\b", re.I)
FORMULATION = re.compile(r"\b(?:tablets?|capsules?|injection|infusion|inhalation suspension|extended release)\b", re.I)
ROUTE_PREFIX = re.compile(r"^(?:intravenous infusion of|inhaled|intravenous|oral)\s+", re.I)


def candidate_transformations(raw: str) -> list[tuple[str, str]]:
    candidates = [(raw, "original")]
    for match in re.findall(r"\(([^()]+)\)", raw):
        candidates.append((match, "parenthetical_alias"))
    cleaned = ROUTE_PREFIX.sub("", raw.strip())
    cleaned = DOSE.sub("", cleaned)
    cleaned = FORMULATION.sub("", cleaned)
    cleaned = re.sub(r"^drug\s+", "", cleaned, flags=re.I)
    cleaned = " ".join(cleaned.split()).strip(" -")
    candidates.append((cleaned, "dose_route_or_formulation_removed"))
    for suffix in SALT_SUFFIXES:
        if cleaned.casefold().endswith(" " + suffix):
            candidates.append((cleaned[: -len(suffix)].strip(), "salt_suffix_removed"))
    seen = set()
    return [(name, rule) for name, rule in candidates if name and not (normalize(name) in seen or seen.add(normalize(name)))]


def main() -> None:
    source = project_path("data", "interim", "mapped_pair_candidates.csv")
    with source.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames or []
        rows = list(reader)
    nodes = read_nodes()
    molecule_dir = project_path("data", "raw", "chembl", "molecule_dictionary")
    summary = json.loads((molecule_dir / "download_summary.json").read_text())
    if not summary["complete_query_download"]:
        raise RuntimeError("Molecule dictionary is incomplete")
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

    changed = 0
    rules = Counter()
    for row in rows:
        row["drug_mapping_transform"] = "original_exact" if row["drug_id"] else ""
        if row["active_intervention_candidate"] != "True" or row["drug_id"]:
            continue
        for candidate, rule in candidate_transformations(row["raw_intervention_name"]):
            direct_matches = unique_matches({candidate}, nodes["drug"])
            chembl_ids = alias_to_chembl.get(normalize(candidate), set())
            chembl_matches = [primekg_by_chembl[identifier] for identifier in chembl_ids if identifier in primekg_by_chembl]
            matches = {match["node_index"]: match for match in direct_matches + chembl_matches}
            if len(matches) != 1:
                continue
            match = next(iter(matches.values()))
            row["drug_id"] = match["node_index"]
            row["drug_node_name"] = match["node_name"]
            row["drug_mapping_status"] = "traceable_synonym_or_cleaned_unique"
            row["drug_mapping_transform"] = f"{rule}:{candidate}"
            row["mapping_status"] = "both_exact_normalized_unique" if row["disease_mapping_status"] == "exact_normalized_unique" else "incomplete_or_ambiguous"
            changed += 1
            rules[rule] += 1
            break

    output = project_path("data", "interim", "mapped_pair_candidates_v2.csv")
    temporary = output.with_suffix(".csv.tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields + ["drug_mapping_transform"])
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(output)
    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source_rows": len(rows),
        "newly_mapped_rows": changed,
        "transformation_counts": dict(rules),
        "method": "Unique exact alias after documented parenthetical, dose/route/formulation or salt transformation; ambiguous matches remain unmapped",
        "warning": "All transformed mappings remain review candidates and are not validated clinical entity mappings",
    }
    report_path = project_path("reports", "e0_negative_synonym_mapping.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

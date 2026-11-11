import csv
import gzip
import json
import re
import sys
import unicodedata
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path


def normalize(value: str) -> str:
    value = unicodedata.normalize("NFKC", value).casefold().replace("&", " and ")
    return " ".join(re.sub(r"[^\w]+", " ", value, flags=re.UNICODE).split())


def read_nodes() -> dict[str, defaultdict[str, list[dict]]]:
    path = project_path("data", "raw", "primekg", "nodes.tab")
    with path.open(encoding="utf-8", newline="") as stream:
        sample = stream.read(4096)
        stream.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        rows = list(csv.DictReader(stream, dialect=dialect))
    nodes = {"drug": defaultdict(list), "disease": defaultdict(list)}
    for row in rows:
        if row.get("node_type") in nodes and row.get("node_name"):
            nodes[row["node_type"]][normalize(row["node_name"])].append(row)
    return nodes


def unique_matches(names: set[str], index: dict[str, list[dict]]) -> list[dict]:
    matches = {}
    for name in names:
        for node in index.get(normalize(name), []):
            matches[node["node_index"]] = node
    return list(matches.values())


def main() -> None:
    molecule_dir = project_path("data", "raw", "chembl", "molecule_dictionary")
    summary = json.loads((molecule_dir / "download_summary.json").read_text(encoding="utf-8"))
    if not summary["complete_query_download"]:
        raise RuntimeError("Molecule dictionary download is incomplete")
    molecules = {}
    for batch in range(1, summary["next_batch"] + 1):
        with gzip.open(molecule_dir / f"batch_{batch:06d}.json.gz", "rt", encoding="utf-8") as stream:
            records = json.load(stream)["molecules"]
        for record in records:
            identifier = record["molecule_chembl_id"]
            if identifier in molecules:
                raise ValueError(f"Duplicate molecule record: {identifier}")
            molecules[identifier] = record

    nodes = read_nodes()
    molecule_mapping = {}
    for identifier, record in molecules.items():
        names = {record.get("pref_name") or ""}
        for synonym in record.get("molecule_synonyms", []):
            for field in ("molecule_synonym", "synonyms"):
                if synonym.get(field):
                    names.add(synonym[field])
        names.discard("")
        matches = unique_matches(names, nodes["drug"])
        status = "exact_or_synonym_unique" if len(matches) == 1 else "ambiguous" if matches else "unmapped"
        molecule_mapping[identifier] = {"status": status, "matches": matches}

    source = project_path("data", "interim", "chembl_approved_indication_candidates.csv")
    with source.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        source_fields = reader.fieldnames or []
        rows = list(reader)
    extra = ["drug_node_index", "drug_node_name", "drug_mapping_status", "disease_node_index", "disease_node_name", "disease_mapping_status"]
    output_rows = []
    for row in rows:
        molecule = molecule_mapping.get(row["parent_molecule_chembl_id"], {"status": "molecule_metadata_missing", "matches": []})
        drug_matches = molecule["matches"]
        disease_matches = unique_matches({row["mesh_heading"], row["efo_term"]}, nodes["disease"])
        row["drug_mapping_status"] = molecule["status"]
        row["disease_mapping_status"] = "exact_mesh_or_efo_unique" if len(disease_matches) == 1 else "ambiguous" if disease_matches else "unmapped"
        row["drug_node_index"] = drug_matches[0]["node_index"] if len(drug_matches) == 1 else ""
        row["drug_node_name"] = drug_matches[0]["node_name"] if len(drug_matches) == 1 else ""
        row["disease_node_index"] = disease_matches[0]["node_index"] if len(disease_matches) == 1 else ""
        row["disease_node_name"] = disease_matches[0]["node_name"] if len(disease_matches) == 1 else ""
        row["mapping_status"] = "both_mapped_unique" if len(drug_matches) == len(disease_matches) == 1 else "incomplete_or_ambiguous"
        output_rows.append(row)

    output = project_path("data", "interim", "chembl_primekg_positive_candidates.csv")
    temporary = output.with_suffix(".csv.tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=source_fields + extra)
        writer.writeheader()
        writer.writerows(output_rows)
    temporary.replace(output)

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "method": "Exact normalized matching using ChEMBL preferred names/synonyms and MeSH/EFO indication names; no fuzzy matching",
        "source_positive_pairs": len(rows),
        "molecule_metadata_records": len(molecules),
        "both_mapped_positive_pairs": sum(row["mapping_status"] == "both_mapped_unique" for row in rows),
        "mapped_unique_drugs": len({row["drug_node_index"] for row in rows if row["drug_node_index"]}),
        "mapped_unique_diseases": len({row["disease_node_index"] for row in rows if row["disease_node_index"]}),
        "drug_mapping_status": dict(Counter(row["drug_mapping_status"] for row in rows)),
        "disease_mapping_status": dict(Counter(row["disease_mapping_status"] for row in rows)),
        "warning": "Mapped positive candidates still require canonical-context review and conflict checks against negative evidence",
    }
    report_path = project_path("reports", "e0_positive_primekg_mapping.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

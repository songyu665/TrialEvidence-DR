import bz2
import csv
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path


def main() -> None:
    with open(project_path("data", "raw", "mondo", "mondo.json"), encoding="utf-8") as stream:
        mondo = json.load(stream)
    mondo_to_doid = defaultdict(set)
    for node in mondo["graphs"][0]["nodes"]:
        identifier = node["id"].rsplit("/", 1)[-1]
        if not identifier.startswith("MONDO_"):
            continue
        identifier = str(int(identifier.split("_", 1)[1]))
        for xref in node.get("meta", {}).get("xrefs", []):
            value = xref.get("val", "")
            if value.startswith("DOID:"):
                mondo_to_doid[identifier].add(value)
        for item in node.get("meta", {}).get("basicPropertyValues", []):
            value = item.get("val", "")
            if "/DOID_" in value:
                mondo_to_doid[identifier].add("DOID:" + value.rsplit("_", 1)[-1])

    with bz2.open(project_path("data", "raw", "hetionet", "hetionet-v1.0.json.bz2"), "rt", encoding="utf-8") as stream:
        graph = json.load(stream)
    compounds = {node["identifier"] for node in graph["nodes"] if node["kind"] == "Compound"}
    diseases = {node["identifier"] for node in graph["nodes"] if node["kind"] == "Disease"}

    with project_path("data", "raw", "primekg", "nodes.tab").open(encoding="utf-8", newline="") as stream:
        sample = stream.read(4096)
        stream.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        nodes = list(csv.DictReader(stream, dialect=dialect))
    primekg = {row["node_index"]: row for row in nodes}

    output_rows = []
    counts = Counter()
    with project_path("data", "interim", "selection_bias_pair_features.csv").open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            drug = primekg[row["drug_node_index"]]
            disease = primekg[row["disease_node_index"]]
            drug_id = drug["node_id"] if drug["node_id"] in compounds else ""
            mondo_id = str(int(disease["node_id"].removeprefix("MONDO_")))
            disease_candidates = mondo_to_doid.get(mondo_id, set()) & diseases
            disease_id = next(iter(disease_candidates)) if len(disease_candidates) == 1 else ""
            counts[(row["evidence_set"], "total")] += 1
            counts[(row["evidence_set"], "drug_id")] += bool(drug_id)
            counts[(row["evidence_set"], "disease_id")] += bool(disease_id)
            counts[(row["evidence_set"], "both")] += bool(drug_id and disease_id)
            if drug_id and disease_id:
                output_rows.append({
                    "evidence_set": row["evidence_set"], "drug_node_index": row["drug_node_index"],
                    "disease_node_index": row["disease_node_index"], "hetionet_compound_id": drug_id,
                    "hetionet_disease_id": disease_id,
                })
    output = project_path("data", "interim", "hetionet_id_mapped_pairs.csv")
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    summary = {}
    for name in ("P", "N_rule", "N_review", "U"):
        total = counts[(name, "total")]
        summary[name] = {
            "total": total,
            "drug_id_mapped": counts[(name, "drug_id")],
            "disease_id_mapped": counts[(name, "disease_id")],
            "both_mapped": counts[(name, "both")],
            "both_fraction": counts[(name, "both")] / total,
        }
    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "mapping": "PrimeKG DrugBank ID to Hetionet Compound plus MONDO to DOID exact cross-reference to Hetionet Disease",
        "coverage": summary,
        "mondo_nodes": len(mondo["graphs"][0]["nodes"]),
        "mondo_ids_with_doid_xref": len(mondo_to_doid),
    }
    with project_path("reports", "f3_hetionet_id_mapping.json").open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

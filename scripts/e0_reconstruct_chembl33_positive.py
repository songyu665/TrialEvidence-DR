import csv
import json
import sqlite3
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from trialevidence.paths import project_path
from e0_map_chembl_primekg import read_nodes, unique_matches


def read_csv(name: str) -> list[dict]:
    with project_path("data", "interim", name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    database = project_path("data", "raw", "chembl", "historical", "chembl_33.db")
    con = sqlite3.connect(f"file:{database}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    rows = con.execute(
        """
        SELECT d.drugind_id, d.molregno, d.mesh_id, d.mesh_heading, d.efo_id, d.efo_term,
               h.parent_molregno, r.ref_type, r.ref_id, r.ref_url
        FROM drug_indication d
        JOIN molecule_hierarchy h ON h.molregno=d.molregno
        JOIN indication_refs r ON r.drugind_id=d.drugind_id
        WHERE d.max_phase_for_ind=4 AND r.ref_type IN ('DailyMed','EMA','FDA')
        """
    ).fetchall()
    parent_ids = {row["parent_molregno"] or row["molregno"] for row in rows}
    dictionaries = {}
    for molregno in parent_ids:
        dictionaries[molregno] = con.execute(
            "SELECT molregno, chembl_id, pref_name, first_approval FROM molecule_dictionary WHERE molregno=?", (molregno,)
        ).fetchone()
    synonyms = defaultdict(set)
    for molregno in parent_ids:
        for row in con.execute("SELECT synonyms FROM molecule_synonyms WHERE molregno=?", (molregno,)):
            if row[0]:
                synonyms[molregno].add(row[0])
    con.close()

    nodes = read_nodes()
    pair_groups = {}
    mapping_counts = defaultdict(int)
    for row in rows:
        parent = row["parent_molregno"] or row["molregno"]
        molecule = dictionaries[parent]
        if molecule is None:
            mapping_counts["missing_parent_dictionary"] += 1
            continue
        names = {molecule["pref_name"] or ""} | synonyms[parent]
        names.discard("")
        drug_matches = unique_matches(names, nodes["drug"])
        disease_matches = unique_matches({row["mesh_heading"] or "", row["efo_term"] or ""}, nodes["disease"])
        if len(drug_matches) != 1 or len(disease_matches) != 1:
            mapping_counts["unmapped_or_ambiguous"] += 1
            continue
        drug = drug_matches[0]
        disease = disease_matches[0]
        key = (drug["node_index"], disease["node_index"])
        if key not in pair_groups:
            pair_groups[key] = {
                "drug_node_index": key[0], "disease_node_index": key[1],
                "drug_node_name": drug["node_name"], "disease_node_name": disease["node_name"],
                "chembl_parent_id": molecule["chembl_id"], "first_approval": molecule["first_approval"],
                "mesh_id": row["mesh_id"], "efo_id": row["efo_id"],
                "drugind_ids": set(), "reference_types": set(), "reference_ids": set(), "reference_urls": set(),
            }
        target = pair_groups[key]
        target["drugind_ids"].add(str(row["drugind_id"]))
        target["reference_types"].add(row["ref_type"])
        if row["ref_id"]:
            target["reference_ids"].add(row["ref_id"])
        if row["ref_url"]:
            target["reference_urls"].add(row["ref_url"])

    output_rows = []
    for row in pair_groups.values():
        for field in ("drugind_ids", "reference_types", "reference_ids", "reference_urls"):
            row[field] = json.dumps(sorted(row[field]), ensure_ascii=False)
        row["historical_version"] = "ChEMBL_33"
        row["candidate_label"] = "P_historical_regulatory_or_label_supported"
        output_rows.append(row)
    output_rows.sort(key=lambda row: (row["disease_node_name"], row["drug_node_name"]))
    output = project_path("data", "interim", "historical_chembl33_positive_pairs.csv")
    with output.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    current = {
        (row["drug_node_index"], row["disease_node_index"])
        for row in read_csv("chembl_primekg_positive_candidates.csv")
        if row["mapping_status"] == "both_mapped_unique"
    }
    historical = set(pair_groups)
    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "historical_version": "ChEMBL_33",
        "phase4_regulatory_reference_rows": len(rows),
        "historical_mapped_unique_positive_pairs": len(historical),
        "historical_unique_parent_molecules": len({row["chembl_parent_id"] for row in output_rows}),
        "historical_unique_diseases": len({row["disease_node_index"] for row in output_rows}),
        "current_mapped_positive_pairs": len(current),
        "current_positive_pairs_present_in_history": len(current & historical),
        "current_positive_recoverability_fraction": len(current & historical) / len(current),
        "historical_pairs_not_in_current_set": len(historical - current),
        "mapping_events_not_uniquely_mapped": dict(mapping_counts),
        "interpretation": "ChEMBL 33 provides a pre-cutoff historical version, but current-versus-historical differences can reflect mapping, curation and later additions",
        "time_experiment_allowed": False,
    }
    report_path = project_path("reports", "e0_chembl33_positive_reconstruction.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

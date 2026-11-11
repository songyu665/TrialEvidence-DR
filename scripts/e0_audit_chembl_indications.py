import csv
import gzip
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path

REGULATORY_REFS = {"DailyMed", "EMA", "FDA"}
FIELDS = [
    "pair_key", "parent_molecule_chembl_id", "efo_id", "efo_term", "mesh_id",
    "mesh_heading", "max_phase_for_ind", "drugind_ids", "reference_types",
    "reference_ids", "reference_urls", "candidate_label", "mapping_status",
]


def main() -> None:
    source_dir = project_path("data", "raw", "chembl", "drug_indication")
    with (source_dir / "download_checkpoint.json").open(encoding="utf-8") as stream:
        state = json.load(stream)
    records = []
    for page in range(1, state["pages_downloaded"] + 1):
        with gzip.open(source_dir / f"page_{page:06d}.json.gz", "rt", encoding="utf-8") as stream:
            records.extend(json.load(stream)["drug_indications"])

    phase_counts = Counter(str(record.get("max_phase_for_ind")) for record in records)
    reference_counts = Counter(
        ref.get("ref_type")
        for record in records
        if str(record.get("max_phase_for_ind")) == "4.0"
        for ref in record.get("indication_refs", [])
    )
    grouped = {}
    for record in records:
        if str(record.get("max_phase_for_ind")) != "4.0":
            continue
        refs = record.get("indication_refs", [])
        selected_refs = [ref for ref in refs if ref.get("ref_type") in REGULATORY_REFS]
        if not selected_refs:
            continue
        parent = record.get("parent_molecule_chembl_id") or record.get("molecule_chembl_id")
        mesh_id = record.get("mesh_id") or ""
        efo_id = record.get("efo_id") or ""
        if not parent or not (mesh_id or efo_id):
            continue
        pair_key = f"{parent}|{mesh_id or efo_id}"
        if pair_key not in grouped:
            grouped[pair_key] = {
                "pair_key": pair_key,
                "parent_molecule_chembl_id": parent,
                "efo_id": efo_id,
                "efo_term": record.get("efo_term") or "",
                "mesh_id": mesh_id,
                "mesh_heading": record.get("mesh_heading") or "",
                "max_phase_for_ind": "4.0",
                "drugind_ids": set(),
                "reference_types": set(),
                "reference_ids": set(),
                "reference_urls": set(),
                "candidate_label": "P_candidate_regulatory_or_label_supported",
                "mapping_status": "unmapped",
            }
        target = grouped[pair_key]
        target["drugind_ids"].add(str(record.get("drugind_id")))
        for ref in selected_refs:
            if ref.get("ref_type"):
                target["reference_types"].add(ref["ref_type"])
            if ref.get("ref_id"):
                target["reference_ids"].add(ref["ref_id"])
            if ref.get("ref_url"):
                target["reference_urls"].add(ref["ref_url"])

    rows = []
    for row in grouped.values():
        for field in ("drugind_ids", "reference_types", "reference_ids", "reference_urls"):
            row[field] = json.dumps(sorted(row[field]), ensure_ascii=False)
        rows.append(row)
    rows.sort(key=lambda row: row["pair_key"])

    output = project_path("data", "interim", "chembl_approved_indication_candidates.csv")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".csv.tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(output)

    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source_records": len(records),
        "source_total_reported_by_api": state["total_count"],
        "source_download_is_complete": state["next_url"] is None,
        "phase_counts": dict(phase_counts),
        "phase4_reference_counts": dict(reference_counts),
        "approved_candidate_pairs": len(rows),
        "unique_parent_molecules": len({row["parent_molecule_chembl_id"] for row in rows}),
        "unique_mesh_ids": len({row["mesh_id"] for row in rows if row["mesh_id"]}),
        "unique_efo_ids": len({row["efo_id"] for row in rows if row["efo_id"]}),
        "selection_rule": "max_phase_for_ind=4 and at least one DailyMed, EMA or FDA indication reference",
        "warning": "Candidates require entity/context mapping and do not by themselves establish universal clinical efficacy",
    }
    report = project_path("reports", "e0_chembl_approved_candidates.json")
    with report.open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

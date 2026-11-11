import csv
import json
import re
import sys
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from trialevidence.paths import project_path
from e0_audit_trial_sample import EXCLUSION_PATTERNS, STRICT_PATTERNS
from e0_build_adjudication_queue import ACTIVE_TYPES, COMBINATION, NONACTIVE
from e0_map_chembl_primekg import normalize, read_nodes, unique_matches

csv.field_size_limit(sys.maxsize)


def read_csv(name: str) -> list[dict]:
    with project_path("data", "interim", name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    known_drug = {}
    known_disease = {}
    for row in read_csv("observed_trial_name_mapping.csv"):
        if not row["node_index"]:
            continue
        if row["node_type"] == "drug":
            known_drug[row["raw_name"]] = (row["node_index"], row["node_name"])
        else:
            known_disease[row["raw_name"]] = (row["node_index"], row["node_name"])
    nodes = read_nodes()

    def map_name(raw: str, node_type: str):
        known = known_drug if node_type == "drug" else known_disease
        if raw in known:
            return known[raw]
        matches = unique_matches({raw}, nodes[node_type])
        return (matches[0]["node_index"], matches[0]["node_name"]) if len(matches) == 1 else (None, None)

    archive = project_path("data", "raw", "aact", "20231201_aact.zip")
    selected = {}
    with zipfile.ZipFile(archive) as z:
        with z.open("studies.txt") as stream:
            wrapper = (line.decode("utf-8-sig" if index == 0 else "utf-8") for index, line in enumerate(stream))
            reader = csv.DictReader(wrapper, delimiter="|")
            for row in reader:
                phase = row.get("phase", "")
                reason = row.get("why_stopped", "").strip()
                if row.get("study_type") != "Interventional" or row.get("overall_status") != "Terminated":
                    continue
                if not ("Phase 2" in phase or "Phase 3" in phase):
                    continue
                strict = any(pattern.search(reason) for pattern in STRICT_PATTERNS.values())
                excluded = any(pattern.search(reason) for pattern in EXCLUSION_PATTERNS.values())
                if strict and not excluded:
                    selected[row["nct_id"]] = {
                        "phase": phase, "why_stopped": reason,
                        "study_first_posted_date": row.get("study_first_posted_date"),
                        "results_first_posted_date": row.get("results_first_posted_date"),
                        "last_update_posted_date": row.get("last_update_posted_date"),
                        "interventions": [], "conditions": [],
                    }
        with z.open("interventions.txt") as stream:
            wrapper = (line.decode("utf-8-sig" if index == 0 else "utf-8") for index, line in enumerate(stream))
            reader = csv.DictReader(wrapper, delimiter="|")
            for row in reader:
                if row["nct_id"] in selected:
                    selected[row["nct_id"]]["interventions"].append({"type": row.get("intervention_type", "").upper(), "name": row.get("name", "")})
        with z.open("conditions.txt") as stream:
            wrapper = (line.decode("utf-8-sig" if index == 0 else "utf-8") for index, line in enumerate(stream))
            reader = csv.DictReader(wrapper, delimiter="|")
            for row in reader:
                if row["nct_id"] in selected:
                    selected[row["nct_id"]]["conditions"].append(row.get("name", ""))

    pair_groups = {}
    skipped = Counter()
    for nct_id, trial in selected.items():
        for intervention in trial["interventions"]:
            name = intervention["name"].strip()
            if intervention["type"] not in ACTIVE_TYPES or not name or NONACTIVE.search(name) or COMBINATION.search(name):
                continue
            drug_id, drug_name = map_name(name, "drug")
            if not drug_id:
                skipped["unmapped_drug"] += 1
                continue
            for condition in trial["conditions"]:
                disease_id, disease_name = map_name(condition, "disease")
                if not disease_id:
                    skipped["unmapped_disease"] += 1
                    continue
                key = (drug_id, disease_id)
                if key not in pair_groups:
                    pair_groups[key] = {
                        "drug_node_index": drug_id, "disease_node_index": disease_id,
                        "drug_node_name": drug_name, "disease_node_name": disease_name,
                        "trial_ids": set(), "raw_evidence": set(), "first_posted_dates": set(),
                    }
                pair_groups[key]["trial_ids"].add(nct_id)
                pair_groups[key]["raw_evidence"].add(trial["why_stopped"])
                pair_groups[key]["first_posted_dates"].add(trial["study_first_posted_date"])

    output_rows = []
    for row in pair_groups.values():
        for field in ("trial_ids", "raw_evidence", "first_posted_dates"):
            row[field] = json.dumps(sorted(row[field]), ensure_ascii=False)
        row["snapshot"] = "AACT_2023-12-01"
        row["label"] = "historical_N_strict_candidate_pending_same_conservative_review"
        output_rows.append(row)
    output_rows.sort(key=lambda row: (row["disease_node_name"], row["drug_node_name"]))
    output = project_path("data", "interim", "historical_20231201_n_strict_candidates.csv")
    with output.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    current = {
        (row["drug_node_index"], row["disease_node_index"])
        for row in read_csv("n_strict_pair_candidates_pending_review.csv")
    }
    historical = set(pair_groups)
    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "snapshot": "AACT_2023-12-01",
        "historical_text_candidate_trials": len(selected),
        "historical_mapped_unique_pair_candidates": len(historical),
        "current_pair_candidates": len(current),
        "current_pairs_present_in_historical_snapshot": len(current & historical),
        "current_pair_recoverability_fraction": len(current & historical) / len(current),
        "skipped_mapping_events": dict(skipped),
        "interpretation": "The snapshot proves availability by 2023-12-01, not the exact first-publication date of every field; historical pairs still require the same conservative review",
        "mapping_note": "Frozen traceable name mappings plus exact PrimeKG fallback were used for entity normalization; no current outcome labels were used",
        "time_experiment_allowed": False,
    }
    report_path = project_path("reports", "e0_aact_20231201_reconstruction.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

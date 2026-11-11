import csv
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


def read_primekg_nodes() -> dict[str, dict[str, list[dict]]]:
    path = project_path("data", "raw", "primekg", "nodes.tab")
    with path.open(encoding="utf-8", newline="") as stream:
        sample = stream.read(4096)
        stream.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        rows = list(csv.DictReader(stream, dialect=dialect))
    result = {"drug": defaultdict(list), "disease": defaultdict(list)}
    for row in rows:
        node_type = row.get("node_type")
        if node_type in result and row.get("node_name"):
            result[node_type][normalize(row["node_name"])].append(row)
    return result


def map_name(name: str, index: dict[str, list[dict]]) -> tuple[str, str]:
    matches = index.get(normalize(name), [])
    if len(matches) == 1:
        return matches[0]["node_index"], "exact_normalized_unique"
    if len(matches) > 1:
        return "", "ambiguous_multiple_primekg_nodes"
    return "", "unmapped"


def main() -> None:
    source = project_path("data", "interim", "label_adjudication_queue.csv")
    with source.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        input_fields = reader.fieldnames or []
        rows = list(reader)
    nodes = read_primekg_nodes()
    extra_fields = ["drug_mapping_status", "disease_mapping_status", "drug_node_name", "disease_node_name"]
    output_fields = input_fields + extra_fields
    for row in rows:
        if row["active_intervention_candidate"] != "True":
            row.update({"drug_mapping_status": "not_applicable", "disease_mapping_status": "not_applicable"})
            row["drug_node_name"] = ""
            row["disease_node_name"] = ""
            continue
        drug_index, drug_status = map_name(row["raw_intervention_name"], nodes["drug"])
        disease_index, disease_status = map_name(row["raw_condition"], nodes["disease"])
        row["drug_id"] = drug_index
        row["graph_node_id"] = disease_index
        row["drug_mapping_status"] = drug_status
        row["disease_mapping_status"] = disease_status
        row["drug_node_name"] = nodes["drug"][normalize(row["raw_intervention_name"])][0]["node_name"] if drug_index else ""
        row["disease_node_name"] = nodes["disease"][normalize(row["raw_condition"])][0]["node_name"] if disease_index else ""
        if drug_status == disease_status == "exact_normalized_unique":
            row["mapping_status"] = "both_exact_normalized_unique"
        else:
            row["mapping_status"] = "incomplete_or_ambiguous"

    output = project_path("data", "interim", "mapped_pair_candidates.csv")
    temporary = output.with_suffix(".csv.tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=output_fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(output)

    active = [row for row in rows if row["active_intervention_candidate"] == "True"]
    noncombination = [row for row in active if row["combination_regimen_flag"] != "True"]
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "method": "Unicode/punctuation normalization followed by exact unique PrimeKG node-name match; no fuzzy matching",
        "source_download_is_complete": json.loads(project_path("data", "raw", "clinicaltrials", "terminated_phase23", "download_summary.json").read_text())["complete_query_download"],
        "all_rows": len(rows),
        "active_rows": len(active),
        "both_mapped_active_rows": sum(row["mapping_status"] == "both_exact_normalized_unique" for row in active),
        "both_mapped_noncombination_rows": sum(row["mapping_status"] == "both_exact_normalized_unique" for row in noncombination),
        "drug_mapping_status_active_rows": dict(Counter(row["drug_mapping_status"] for row in active)),
        "disease_mapping_status_active_rows": dict(Counter(row["disease_mapping_status"] for row in active)),
        "warning": "Exact node-name mapping does not validate intervention attribution, canonical context or label semantics",
    }
    report = project_path("reports", "e0_primekg_exact_mapping.json")
    with report.open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

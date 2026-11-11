import csv
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.control import RunController
from trialevidence.paths import project_path

CUTOFF = "2023-12-31"


def read_csv(name: str) -> list[dict]:
    with project_path("data", "interim", name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def fetch(setid: str) -> bytes:
    url = "https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json?" + urllib.parse.urlencode({"setid": setid})
    request = urllib.request.Request(url, headers={"User-Agent": "TrialEvidence-DR-E0/0.1", "Accept": "application/json"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=45) as response:
                return response.read()
        except Exception:
            if attempt == 4:
                raise
            time.sleep(2**attempt)
    raise RuntimeError("Unreachable retry state")


def main() -> None:
    candidates = [row for row in read_csv("post_cutoff_positive_event_candidates.csv") if row["event_category"] == "new_indication_records_after_2023_version"]
    setid_to_pairs = defaultdict(set)
    pair_info = {}
    for row in candidates:
        pair_id = f"{row['drug_node_index']}|{row['disease_node_index']}"
        pair_info[pair_id] = row
        for url in json.loads(row["reference_urls"]):
            match = re.search(r"[?&]setid=([0-9a-f-]+)", url)
            if match:
                setid_to_pairs[match.group(1)].add(pair_id)

    controller = RunController("e0_dailymed_label_dates")
    controller.install_signal_handlers()
    output_dir = project_path("data", "raw", "dailymed", "spls")
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata_rows = []
    setids = sorted(setid_to_pairs)
    for index, setid in enumerate(setids):
        reason = controller.stop_reason
        if reason:
            controller.save_checkpoint({"next_index": index, "total": len(setids)})
            raise SystemExit(85)
        path = output_dir / f"{setid}.json"
        if path.exists():
            body = path.read_bytes()
        else:
            body = fetch(setid)
            temporary = path.with_suffix(".json.tmp")
            temporary.write_bytes(body)
            temporary.replace(path)
        payload = json.loads(body)
        records = payload.get("data", [])
        if records:
            for record in records:
                metadata_rows.append({
                    "setid": setid, "published_date": record.get("published_date"),
                    "spl_version": record.get("spl_version"), "title": record.get("title"),
                    "pair_ids": json.dumps(sorted(setid_to_pairs[setid])),
                    "source_url": payload.get("metadata", {}).get("current_url"),
                })
        if (index + 1) % 25 == 0:
            print(json.dumps({"processed": index + 1, "total": len(setids)}), flush=True)

    metadata_output = project_path("data", "interim", "dailymed_label_metadata.csv")
    with metadata_output.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(metadata_rows[0]))
        writer.writeheader()
        writer.writerows(metadata_rows)

    pair_dates = defaultdict(list)
    for row in metadata_rows:
        if not row["published_date"]:
            continue
        try:
            parsed = datetime.strptime(row["published_date"], "%B %d, %Y").date().isoformat()
        except ValueError:
            parsed = datetime.strptime(row["published_date"], "%b %d, %Y").date().isoformat()
        for pair_id in json.loads(row["pair_ids"]):
            pair_dates[pair_id].append(parsed)
    output_rows = []
    for pair_id in sorted(pair_info):
        dates = sorted(pair_dates.get(pair_id, []))
        if not dates:
            status = "unresolved_no_current_spl_metadata"
        elif any(date > CUTOFF for date in dates):
            status = "has_post_cutoff_label_publication_not_indication_approval"
        else:
            status = "all_available_labels_published_before_cutoff"
        source = pair_info[pair_id]
        output_rows.append({
            "pair_id": pair_id, "drug_node_name": source["drug_node_name"],
            "disease_node_name": source["disease_node_name"],
            "label_published_dates": json.dumps(dates), "date_audit_status": status,
            "indication_approval_date_confirmed": False,
        })
    output = project_path("data", "interim", "dailymed_label_date_audit.csv")
    with output.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "new_record_candidate_pairs": len(pair_info), "unique_dailymed_setids": len(setids),
        "setids_with_current_spl_metadata": len({row["setid"] for row in metadata_rows}),
        "pair_date_audit_counts": dict(Counter(row["date_audit_status"] for row in output_rows)),
        "date_meaning": "DailyMed SPL published_date, not indication-specific regulatory approval",
        "time_experiment_allowed": False,
    }
    report_path = project_path("reports", "e0_dailymed_label_date_audit.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    controller.mark_completed({"pairs": len(output_rows), "setids": len(setids)})
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

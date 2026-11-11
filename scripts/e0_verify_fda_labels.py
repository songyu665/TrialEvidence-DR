import csv
import hashlib
import json
import re
import subprocess
import sys
import time
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.control import RunController
from trialevidence.paths import project_path


def fetch(url: str) -> bytes:
    if url.startswith("http://"):
        url = "https://" + url.removeprefix("http://")
    request = urllib.request.Request(url, headers={"User-Agent": "TrialEvidence-DR-E0/0.1"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                return response.read()
        except Exception:
            if attempt == 4:
                raise
            time.sleep(2**attempt)
    raise RuntimeError("Unreachable retry state")


def main() -> None:
    with project_path("data", "interim", "post_cutoff_positive_event_candidates.csv").open(encoding="utf-8-sig", newline="") as stream:
        candidates = list(csv.DictReader(stream))
    url_to_pairs = defaultdict(list)
    pair_info = {}
    for row in candidates:
        pair_id = f"{row['drug_node_index']}|{row['disease_node_index']}"
        pair_info[pair_id] = row
        for url in json.loads(row["reference_urls"]):
            if "fda.gov" in url.lower() and ".pdf" in url.lower():
                url_to_pairs[url].append(pair_id)

    controller = RunController("e0_fda_label_metadata")
    controller.install_signal_handlers()
    raw_dir = project_path("data", "raw", "fda", "labels")
    text_dir = project_path("data", "raw", "fda", "label_text")
    raw_dir.mkdir(parents=True, exist_ok=True)
    text_dir.mkdir(parents=True, exist_ok=True)
    metadata = {}
    urls = sorted(url_to_pairs)
    for index, url in enumerate(urls):
        if controller.stop_reason:
            controller.save_checkpoint({"next_index": index, "total": len(urls)})
            raise SystemExit(85)
        key = hashlib.sha256(url.encode()).hexdigest()[:16]
        path = raw_dir / f"{key}.pdf"
        body = path.read_bytes() if path.exists() else fetch(url)
        if not body.startswith(b"%PDF"):
            raise ValueError(f"Not a PDF: {url}")
        if not path.exists():
            temporary = path.with_suffix(".pdf.tmp")
            temporary.write_bytes(body)
            temporary.replace(path)
        process = subprocess.run(["pdftotext", "-f", "1", "-l", "3", str(path), "-"], capture_output=True, check=True)
        text = process.stdout.decode("utf-8", "replace")
        (text_dir / f"{key}.txt").write_text(text, encoding="utf-8")
        initial = re.search(r"Initial U\.S\. Approval[:\s]+(?:\(.*?\)\s*)?(\d{4})", text, re.I)
        revised = re.search(r"\bRevised[:\s]+((?:\d{1,2}/)?(?:\d{4})|(?:[A-Z][a-z]+\s+\d{4}))", text, re.I)
        url_year = re.search(r"/label/(\d{4})/", url)
        metadata[url] = {
            "initial_us_approval_year": initial.group(1) if initial else "",
            "revised_date_text": revised.group(1) if revised else "",
            "label_url_year": url_year.group(1) if url_year else "",
            "pdf_sha256": hashlib.sha256(body).hexdigest(),
            "text_file": str((text_dir / f"{key}.txt").relative_to(project_path())),
        }
        print(json.dumps({"processed": index + 1, "total": len(urls), **metadata[url]}), flush=True)

    output_rows = []
    for url, pairs in sorted(url_to_pairs.items()):
        for pair_id in pairs:
            row = pair_info[pair_id]
            output_rows.append({
                "pair_id": pair_id, "drug_node_name": row["drug_node_name"],
                "disease_node_name": row["disease_node_name"], "fda_label_url": url,
                **metadata[url], "indication_event_date_confirmed": False,
                "date_interpretation": "Label metadata only; not indication-specific event confirmation",
            })
    output = project_path("data", "interim", "fda_label_metadata.csv")
    with output.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "candidate_pairs_with_fda_labels": len({row["pair_id"] for row in output_rows}),
        "unique_fda_label_urls": len(urls),
        "with_initial_us_approval_year": sum(bool(item["initial_us_approval_year"]) for item in metadata.values()),
        "with_revised_date": sum(bool(item["revised_date_text"]) for item in metadata.values()),
        "label_url_year_counts": dict(Counter(item["label_url_year"] for item in metadata.values())),
        "indication_event_date_confirmed": False, "time_experiment_allowed": False,
    }
    with project_path("reports", "e0_fda_label_metadata.json").open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    controller.mark_completed({"pairs": report["candidate_pairs_with_fda_labels"], "urls": len(urls)})
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

import csv
import json
import re
import sys
import time
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.control import RunController
from trialevidence.paths import project_path


class TextParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.parts = []

    def handle_data(self, data):
        text = " ".join(data.split())
        if text:
            self.parts.append(text)


def fetch(url: str) -> bytes:
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


def after(parts: list[str], label: str, pattern: str) -> str:
    for index, part in enumerate(parts):
        if label in part:
            for value in parts[index + 1:index + 5]:
                match = re.search(pattern, value)
                if match:
                    return match.group()
    return ""


def main() -> None:
    with project_path("data", "interim", "post_cutoff_positive_event_candidates.csv").open(encoding="utf-8-sig", newline="") as stream:
        candidates = list(csv.DictReader(stream))
    url_to_pairs = defaultdict(list)
    pair_info = {}
    for row in candidates:
        pair_id = f"{row['drug_node_index']}|{row['disease_node_index']}"
        pair_info[pair_id] = row
        for url in json.loads(row["reference_urls"]):
            if "ema.europa.eu" in url and "/EPAR/" in url:
                url_to_pairs[url].append(pair_id)

    controller = RunController("e0_ema_event_metadata")
    controller.install_signal_handlers()
    raw_dir = project_path("data", "raw", "ema", "epar")
    raw_dir.mkdir(parents=True, exist_ok=True)
    metadata = {}
    urls = sorted(url_to_pairs)
    for index, url in enumerate(urls):
        if controller.stop_reason:
            controller.save_checkpoint({"next_index": index, "total": len(urls)})
            raise SystemExit(85)
        key = re.sub(r"[^A-Za-z0-9_.-]+", "_", url.rsplit("/", 1)[-1])
        path = raw_dir / f"{key}.html"
        body = path.read_bytes() if path.exists() else fetch(url)
        if not path.exists():
            temporary = path.with_suffix(".html.tmp")
            temporary.write_bytes(body)
            temporary.replace(path)
        parser = TextParser()
        parser.feed(body.decode("utf-8", "replace"))
        metadata[url] = {
            "marketing_authorisation_issued": after(parser.parts, "Marketing authorisation issued", r"\b\d{2}/\d{2}/\d{4}\b"),
            "page_last_updated": after(parser.parts, "This page was last updated on", r"\b\d{2}/\d{2}/\d{4}\b"),
            "raw_html_file": str(path.relative_to(project_path())),
        }
        if (index + 1) % 10 == 0:
            print(json.dumps({"processed": index + 1, "total": len(urls)}), flush=True)

    output_rows = []
    for url, pairs in sorted(url_to_pairs.items()):
        for pair_id in pairs:
            row = pair_info[pair_id]
            output_rows.append({
                "pair_id": pair_id, "drug_node_name": row["drug_node_name"],
                "disease_node_name": row["disease_node_name"], "epar_url": url,
                **metadata[url],
                "indication_event_date_confirmed": False,
                "date_interpretation": "Initial marketing authorisation and page-update dates, not indication-specific event dates",
            })
    output = project_path("data", "interim", "ema_event_metadata.csv")
    with output.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "candidate_pairs_with_ema_metadata": len({row["pair_id"] for row in output_rows}),
        "unique_epar_urls": len(urls),
        "with_initial_authorisation_date": sum(bool(item["marketing_authorisation_issued"]) for item in metadata.values()),
        "with_page_update_date": sum(bool(item["page_last_updated"]) for item in metadata.values()),
        "indication_event_date_confirmed": False,
        "time_experiment_allowed": False,
    }
    with project_path("reports", "e0_ema_event_metadata.json").open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    controller.mark_completed({"pairs": report["candidate_pairs_with_ema_metadata"], "urls": len(urls)})
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

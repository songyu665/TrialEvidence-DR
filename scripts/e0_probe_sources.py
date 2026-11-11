import json
import ssl
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path


def probe(source: dict) -> dict:
    request = urllib.request.Request(
        source["url"],
        method=source.get("method", "GET"),
        headers={"User-Agent": "TrialEvidence-DR-E0/0.1", "Accept": "application/json,text/html,*/*"},
    )
    started = time.monotonic()
    result = {"name": source["name"], "role": source["role"], "url": source["url"]}
    try:
        with urllib.request.urlopen(request, timeout=45, context=ssl.create_default_context()) as response:
            parse_json = source["name"] in {"ClinicalTrials.gov_API_v2", "PrimeKG_Dataverse_metadata"}
            if parse_json:
                body = response.read()
            else:
                body = response.read(2048) if request.get_method() == "GET" else b""
            result.update(
                {
                    "reachable": True,
                    "http_status": response.status,
                    "content_type": response.headers.get("Content-Type"),
                    "content_length": response.headers.get("Content-Length"),
                    "bytes_read": len(body),
                }
            )
            if source["name"] == "ClinicalTrials.gov_API_v2":
                payload = json.loads(body.decode("utf-8"))
                result["reported_total_studies"] = payload.get("totalCount")
            elif source["name"] == "PrimeKG_Dataverse_metadata":
                payload = json.loads(body.decode("utf-8"))
                files = payload.get("data", {}).get("latestVersion", {}).get("files", [])
                result["metadata_file_count"] = len(files)
                result["data_file_names"] = [item.get("dataFile", {}).get("filename") for item in files]
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError) as error:
        result.update({"reachable": False, "error": f"{type(error).__name__}: {error}"})
    result["elapsed_seconds"] = round(time.monotonic() - started, 3)
    return result


def main() -> None:
    config_path = project_path("configs", "data_sources.json")
    with config_path.open(encoding="utf-8") as stream:
        sources = json.load(stream)["sources"]
    results = [probe(source) for source in sources]
    payload = {
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "scope": "connectivity_and_metadata_only",
        "results": results,
    }
    output = project_path("reports", "e0_source_probe.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".json.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    temporary.replace(output)
    for result in results:
        print(json.dumps(result, ensure_ascii=False))
    if not all(result.get("reachable") for result in results):
        raise SystemExit(1)


if __name__ == "__main__":
    main()

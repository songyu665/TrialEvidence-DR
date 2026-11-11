import csv
import gzip
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path

STRICT_PATTERNS = {
    "explicit_futility": re.compile(r"\bfutility\b|\bfutile\b", re.I),
    "explicit_lack_of_efficacy": re.compile(r"lack of efficacy|lack of effectiveness|insufficient efficacy|no efficacy|inefficacy|ineffective|not effective|no (?:clinical |treatment )?benefit", re.I),
}
REVIEW_PATTERN = re.compile(r"primary (?:efficacy )?endpoint|co-primary|failed to meet|did not meet|not met", re.I)
EXCLUSION_PATTERNS = {
    "safety": re.compile(r"safety|toxicity|adverse|tolerability", re.I),
    "recruitment": re.compile(r"recruit|enrol|accrual|participant availability", re.I),
    "administrative_or_business": re.compile(r"business|commercial|funding|financial|sponsor decision|administrative|strategic", re.I),
}


def main() -> None:
    input_dir = project_path("data", "raw", "clinicaltrials", "terminated_phase23")
    with (input_dir / "download_summary.json").open(encoding="utf-8") as stream:
        download = json.load(stream)
    records = []
    for page in range(1, download["pages_downloaded"] + 1):
        with gzip.open(input_dir / f"page_{page:06d}.json.gz", "rt", encoding="utf-8") as stream:
            records.extend(json.load(stream)["studies"])

    candidates = []
    reason_counts = Counter()
    phase_counts = Counter()
    condition_counts = Counter()
    missing = Counter()
    results_available = 0
    for study in records:
        protocol = study.get("protocolSection", {})
        identification = protocol.get("identificationModule", {})
        status = protocol.get("statusModule", {})
        design = protocol.get("designModule", {})
        conditions_module = protocol.get("conditionsModule", {})
        interventions_module = protocol.get("armsInterventionsModule", {})
        nct_id = identification.get("nctId", "")
        reason = status.get("whyStopped", "").strip()
        if not reason:
            missing["whyStopped"] += 1
        conditions = conditions_module.get("conditions", [])
        if not conditions:
            missing["conditions"] += 1
        interventions = interventions_module.get("interventions", [])
        if not interventions:
            missing["interventions"] += 1
        phases = design.get("phases", [])
        if not phases:
            missing["phases"] += 1
        if "resultsSection" in study:
            results_available += 1
        for phase in phases:
            phase_counts[phase] += 1
        for condition in conditions:
            condition_counts[condition.strip().casefold()] += 1

        strict_matches = [name for name, pattern in STRICT_PATTERNS.items() if pattern.search(reason)]
        exclusion_matches = [name for name, pattern in EXCLUSION_PATTERNS.items() if pattern.search(reason)]
        review_match = bool(REVIEW_PATTERN.search(reason))
        if strict_matches and not exclusion_matches:
            bucket = "N_strict_candidate_requires_three_layer_review"
        elif strict_matches:
            bucket = "conflicting_reason_candidates_exclude_pending_review"
        elif review_match:
            bucket = "N_review_candidate_not_automatic_negative"
        elif exclusion_matches:
            bucket = "N_excl_candidate"
        elif reason:
            bucket = "other_or_unclear"
        else:
            bucket = "missing_reason"
        reason_counts[bucket] += 1
        candidates.append(
            {
                "nct_id": nct_id,
                "phases": json.dumps(phases, ensure_ascii=False),
                "conditions": json.dumps(conditions, ensure_ascii=False),
                "interventions": json.dumps(
                    [{"type": item.get("type"), "name": item.get("name")} for item in interventions],
                    ensure_ascii=False,
                ),
                "overall_status": status.get("overallStatus"),
                "why_stopped": reason,
                "results_available": "resultsSection" in study,
                "start_date": status.get("startDateStruct", {}).get("date"),
                "completion_date": status.get("completionDateStruct", {}).get("date"),
                "results_first_posted": status.get("resultsFirstPostDateStruct", {}).get("date"),
                "candidate_bucket": bucket,
                "strict_text_signals": json.dumps(strict_matches),
                "exclusion_text_signals": json.dumps(exclusion_matches),
                "requires_manual_review": True,
            }
        )

    ids = [row["nct_id"] for row in candidates]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate NCT IDs in sample")
    output_csv = project_path("data", "interim", "trial_label_candidates.csv")
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(candidates[0]))
        writer.writeheader()
        writer.writerows(candidates)

    complete = download["complete_query_download"]
    limits = [
        "Text candidates are not adjudicated N_strict labels",
        "Intervention attribution and indication-context aggregation remain unverified",
    ]
    if complete:
        sample_design = "Complete census of the explicit terminated phase II/III API query"
        limits.append("Results do not cover completed trials with negative findings or trials outside this query")
    else:
        sample_design = "First consecutive API pages; not a random or representative sample"
        limits.append("Consecutive-page frequencies must not be extrapolated to the full query population")
    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "sample_design": sample_design,
        "query_total_count": download["total_count"],
        "records_audited": len(records),
        "results_section_available": results_available,
        "missing_fields": dict(missing),
        "phase_counts": dict(phase_counts),
        "distinct_raw_conditions_casefolded": len(condition_counts),
        "candidate_buckets": dict(reason_counts),
        "interpretation_limits": limits,
    }
    output_json = project_path("reports", "e0_trial_sample_audit.json")
    with output_json.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

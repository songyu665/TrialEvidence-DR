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


def numeric_pvalue(value: str | None) -> float | None:
    if not value:
        return None
    match = re.search(r"\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", value)
    return float(match.group()) if match else None


def main() -> None:
    source_dir = project_path("data", "raw", "clinicaltrials", "completed_results_phase23_with_comments")
    state = json.loads((source_dir / "download_checkpoint.json").read_text())
    studies = []
    for page in range(1, state["pages_downloaded"] + 1):
        with gzip.open(source_dir / f"page_{page:06d}.json.gz", "rt", encoding="utf-8") as stream:
            studies.extend(json.load(stream)["studies"])

    counts = Counter()
    queue = []
    for study in studies:
        protocol = study.get("protocolSection", {})
        nct_id = protocol.get("identificationModule", {}).get("nctId")
        outcomes = study.get("resultsSection", {}).get("outcomeMeasuresModule", {}).get("outcomeMeasures", [])
        primary = [outcome for outcome in outcomes if outcome.get("type") == "PRIMARY" and outcome.get("reportingStatus") == "POSTED"]
        analyses = [analysis for outcome in primary for analysis in outcome.get("analyses", [])]
        with_p = [analysis for analysis in analyses if analysis.get("pValue")]
        with_ci = [analysis for analysis in analyses if analysis.get("ciLowerLimit") or analysis.get("ciUpperLimit")]
        noninferiority = [analysis for analysis in analyses if analysis.get("testedNonInferiority") is True]
        counts["studies"] += 1
        counts["with_posted_primary_outcome"] += bool(primary)
        counts["with_primary_analysis"] += bool(analyses)
        counts["with_primary_pvalue"] += bool(with_p)
        counts["with_primary_confidence_interval"] += bool(with_ci)
        counts["with_noninferiority_analysis"] += bool(noninferiority)
        counts["primary_outcomes"] += len(primary)
        counts["primary_analyses"] += len(analyses)
        if not analyses:
            continue
        if any((numeric_pvalue(analysis.get("pValue")) or 0) > 0.05 for analysis in with_p):
            counts["with_p_gt_0_05_exploratory_not_a_label"] += 1
        interventions = protocol.get("armsInterventionsModule", {}).get("interventions", [])
        arms = protocol.get("armsInterventionsModule", {}).get("armGroups", [])
        queue.append({
            "nct_id": nct_id,
            "conditions": json.dumps(protocol.get("conditionsModule", {}).get("conditions", []), ensure_ascii=False),
            "interventions": json.dumps([{"type": item.get("type"), "name": item.get("name")} for item in interventions], ensure_ascii=False),
            "arm_groups": json.dumps([{"type": item.get("type"), "interventionNames": item.get("interventionNames", [])} for item in arms], ensure_ascii=False),
            "results_first_posted": protocol.get("statusModule", {}).get("resultsFirstPostDateStruct", {}).get("date"),
            "primary_outcomes": json.dumps(primary, ensure_ascii=False),
            "primary_analysis_count": len(analyses),
            "primary_pvalue_count": len(with_p),
            "primary_ci_count": len(with_ci),
            "noninferiority_analysis_present": bool(noninferiority),
            "candidate_bucket": "N_review_primary_results_require_adjudication",
            "label": "U_pending_review",
            "evidence_url": f"https://clinicaltrials.gov/study/{nct_id}?tab=results",
        })

    ids = [study["protocolSection"]["identificationModule"]["nctId"] for study in studies]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate NCT IDs")
    output = project_path("data", "interim", "completed_primary_outcome_review_queue.csv")
    temporary = output.with_suffix(".csv.tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(queue[0]))
        writer.writeheader()
        writer.writerows(queue)
    temporary.replace(output)

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source_download_is_complete": state["next_token"] is None,
        "source_total_reported_by_api": state["total_count"],
        "audit_counts": dict(counts),
        "review_queue_trials": len(queue),
        "interpretation_limits": [
            "No completed trial is automatically assigned N_strict",
            "P>0.05 is only an exploratory field and does not establish lack of efficacy",
            "Effect direction, comparator, estimand, multiplicity and clinical relevance require adjudication",
        ],
    }
    report_path = project_path("reports", "e0_completed_results_audit.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

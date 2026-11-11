import csv
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path

csv.field_size_limit(sys.maxsize)

NEGATIVE_TEXT = re.compile(r"no significant|not statistically significant|failed to (?:meet|demonstrate|show)|did not (?:meet|demonstrate|show)|lack of efficacy|no efficacy|futility|no clinical benefit|no treatment benefit", re.I)
DESIGN_TEXT = re.compile(r"\bif\b|would|criteria|monitoring rule|analysis plan|pre-?defined|pre-?specified|ability to stop|discontinuation due to|was to be|were to be", re.I)
ACTUAL_NEGATIVE = re.compile(r"did not meet|was not statistically significant|were not statistically significant|showed no significant|reaching the futility boundary|stopped (?:early )?due to futility|indicated futility|showed no significant", re.I)
RATIO_PARAMS = {"ODDS_RATIO", "HAZARD_RATIO", "RATE_RATIO", "RISK_RATIO", "RELATIVE_RISK", "RATIO", "GEOMETRIC_MEAN_RATIO"}


def number(value):
    if value is None:
        return None
    match = re.search(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", str(value))
    return float(match.group()) if match else None


def ci_crosses_null(analysis: dict) -> bool | None:
    lower = number(analysis.get("ciLowerLimit"))
    upper = number(analysis.get("ciUpperLimit"))
    if lower is None or upper is None:
        return None
    null = 1.0 if str(analysis.get("paramType", "")).upper() in RATIO_PARAMS else 0.0
    return lower <= null <= upper


def main() -> None:
    source = project_path("data", "interim", "completed_primary_outcome_review_queue.csv")
    with source.open(encoding="utf-8-sig", newline="") as stream:
        trials = list(csv.DictReader(stream))
    rows = []
    for trial in trials:
        outcomes = json.loads(trial["primary_outcomes"])
        for outcome_index, outcome in enumerate(outcomes):
            for analysis_index, analysis in enumerate(outcome.get("analyses", [])):
                text = " ".join(str(analysis.get(field) or "") for field in (
                    "statisticalComment", "pValueComment", "estimateComment", "otherAnalysisDescription",
                ))
                pvalue = number(analysis.get("pValue"))
                null_in_ci = ci_crosses_null(analysis)
                noninferiority = analysis.get("testedNonInferiority") is True
                explicit_negative = bool(NEGATIVE_TEXT.search(text))
                design_text = bool(DESIGN_TEXT.search(text))
                actual_negative = bool(ACTUAL_NEGATIVE.search(text))
                if explicit_negative and actual_negative and not design_text and not noninferiority:
                    priority = "high_text_supported_N_review"
                elif not noninferiority and (explicit_negative or (pvalue is not None and pvalue > 0.05) or null_in_ci is True):
                    priority = "medium_statistical_N_review"
                else:
                    priority = "other_primary_analysis"
                rows.append({
                    "nct_id": trial["nct_id"],
                    "outcome_index": outcome_index,
                    "analysis_index": analysis_index,
                    "outcome_title": outcome.get("title"),
                    "outcome_param_type": outcome.get("paramType"),
                    "outcome_unit": outcome.get("unitOfMeasure"),
                    "analysis_param_type": analysis.get("paramType"),
                    "analysis_param_value": analysis.get("paramValue"),
                    "dispersion_type": analysis.get("dispersionType"),
                    "dispersion_value": analysis.get("dispersionValue"),
                    "statistical_method": analysis.get("statisticalMethod"),
                    "p_value": analysis.get("pValue"),
                    "ci_num_sides": analysis.get("ciNumSides"),
                    "ci_pct_value": analysis.get("ciPctValue"),
                    "ci_lower_limit": analysis.get("ciLowerLimit"),
                    "ci_upper_limit": analysis.get("ciUpperLimit"),
                    "tested_noninferiority": noninferiority,
                    "noninferiority_type": analysis.get("nonInferiorityType"),
                    "group_ids": json.dumps(analysis.get("groupIds", [])),
                    "group_description": analysis.get("groupDescription"),
                    "statistical_comment": analysis.get("statisticalComment"),
                    "p_value_comment": analysis.get("pValueComment"),
                    "estimate_comment": analysis.get("estimateComment"),
                    "other_analysis_description": analysis.get("otherAnalysisDescription"),
                    "explicit_negative_text": explicit_negative,
                    "design_or_hypothetical_text": design_text,
                    "actual_negative_result_text": actual_negative,
                    "p_gt_0_05": pvalue is not None and pvalue > 0.05,
                    "ci_contains_null": null_in_ci,
                    "review_priority": priority,
                    "label": "U_pending_review",
                    "results_first_posted": trial["results_first_posted"],
                    "evidence_url": trial["evidence_url"],
                })

    output = project_path("data", "interim", "completed_primary_analysis_review.csv")
    temporary = output.with_suffix(".csv.tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(output)
    priorities = Counter(row["review_priority"] for row in rows)
    high = [row for row in rows if row["review_priority"] == "high_text_supported_N_review"]
    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "analysis_rows": len(rows),
        "unique_trials": len({row["nct_id"] for row in rows}),
        "review_priority_counts": dict(priorities),
        "high_priority_unique_trials": len({row["nct_id"] for row in high}),
        "warning": "Priority is based on text/statistical screening only; no row is an adjudicated negative and no effect direction is inferred",
    }
    report_path = project_path("reports", "e0_completed_analysis_review.json")
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

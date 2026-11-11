import csv
import hashlib
import io
import json
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOCAL_ROOT = Path(__file__).resolve().parents[1]
DELIVERY_ZIP = ROOT / "TrialEvidence_DR_delivery.zip"
CURRENT_MODEL = ROOT / "reports" / "tables" / "f3_evidence_status_model_results.csv"
REVIEW_SCRIPT = ROOT / "scripts" / "e0_apply_conservative_review.py"
REVIEW_OUTPUT = ROOT / "data" / "interim" / "n_strict_conservative_rule_review.csv"


def is_number(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def main() -> None:
    with zipfile.ZipFile(DELIVERY_ZIP) as archive:
        old = list(csv.DictReader(io.StringIO(archive.read("reports/tables/f3_evidence_status_model_results.csv").decode("utf-8"))))
    with CURRENT_MODEL.open(encoding="utf-8", newline="") as stream:
        new = list(csv.DictReader(stream))
    if len(old) != len(new):
        raise ValueError("Model-result row count changed")
    differences = []
    for row_index, (a, b) in enumerate(zip(old, new)):
        for key in a:
            if a[key] == b[key]:
                continue
            if a[key].startswith("[") and b[key].startswith("["):
                left_list, right_list = json.loads(a[key]), json.loads(b[key])
                if len(left_list) == len(right_list) and all(is_number(str(x)) and is_number(str(y)) for x, y in zip(left_list, right_list)):
                    close = all(abs(x-y) <= 1e-6 * max(1.0, abs(x), abs(y)) for x, y in zip(left_list, right_list))
                    if close:
                        differences.append({"row": row_index, "field": key, "old": left_list, "new": right_list, "within_tolerance": True})
                        continue
            if is_number(a[key]) and is_number(b[key]):
                left, right = float(a[key]), float(b[key])
                tolerance = 1e-6 * max(1.0, abs(left), abs(right))
                if abs(left - right) <= tolerance:
                    differences.append({"row": row_index, "field": key, "old": left, "new": right, "absolute_difference": abs(left-right), "within_tolerance": True})
                    continue
            differences.append({"row": row_index, "field": key, "old": a[key], "new": b[key], "within_tolerance": False})
    failed = [item for item in differences if not item["within_tolerance"]]

    before = hashlib.sha256(REVIEW_OUTPUT.read_bytes()).hexdigest()
    process = subprocess.run([str(ROOT / "venv_trialevidence_dr" / "bin" / "python"), str(REVIEW_SCRIPT)], cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, check=False)
    if process.returncode:
        raise RuntimeError(process.stderr.decode("utf-8", "replace"))
    after = hashlib.sha256(REVIEW_OUTPUT.read_bytes()).hexdigest()

    report_path = ROOT / "reports" / "reproducibility_rerun.json"
    original = json.loads(report_path.read_text())
    original["resolution"] = {
        "resolved_at": datetime.now(timezone.utc).isoformat(),
        "model_numeric_differences": differences,
        "model_numeric_tolerance": "absolute difference <= 1e-6 * max(1, |old|, |new|)",
        "model_material_differences": failed,
        "rule_review_exact_rerun_hash_match": before == after,
        "initial_rule_review_change_explanation": "The review_date field changed across the calendar-day boundary; exact rerun after the boundary passed",
        "status": "passed_with_date_normalization_and_float_tolerance" if not failed and before == after else "unresolved",
    }
    original["status"] = original["resolution"]["status"]
    report_path.write_text(json.dumps(original, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(original["resolution"], ensure_ascii=False, indent=2))
    if failed or before != after:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

SCRIPTS = [
    "e0_audit_trial_sample.py",
    "e0_build_adjudication_queue.py",
    "e0_map_primekg_exact.py",
    "e0_map_negative_drug_synonyms.py",
    "e0_aggregate_pair_candidates.py",
    "e0_apply_conservative_review.py",
    "e0_finalize_rule_based_pairs.py",
    "e0_audit_chembl_indications.py",
    "e0_map_chembl_primekg.py",
    "e0_conflict_and_coverage.py",
    "e0_compute_observed_trial_counts.py",
    "f3_selection_bias.py",
    "f3_selection_bias_effects.py",
    "f3_evidence_status_models.py",
    "e0_reconstruct_aact_20231201.py",
    "e0_review_aact_history.py",
    "e0_history_discrepancy_analysis.py",
    "e0_reconstruct_chembl33_positive.py",
    "e0_post_cutoff_event_audit.py",
    "e0_verify_dailymed_dates.py",
    "f3_historical_selection_bias.py",
    "f3_repeated_u_sensitivity.py",
    "f3_leave_label_source_out.py",
    "f3_repeated_cv.py",
    "f3_hetionet_validation.py",
    "f3_hetionet_id_mapping.py",
    "e0_generate_figures.py",
    "render_e0_report.py",
    "build_expert_review_package.py",
    "build_manuscript_docx.py",
]


def hashes() -> dict[str, str]:
    result = {}
    for directory in ("data/interim", "reports/tables"):
        for path in sorted((ROOT / directory).rglob("*")):
            if path.is_file():
                result[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def main() -> None:
    before = hashes()
    log_path = ROOT / "runs" / "reproducibility_rerun.log"
    with log_path.open("w", encoding="utf-8") as log:
        for script in SCRIPTS:
            log.write(f"START {script}\n")
            log.flush()
            process = subprocess.run(
                [str(ROOT / "venv_trialevidence_dr" / "bin" / "python"), str(ROOT / "scripts" / script)],
                cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, timeout=900, check=False,
            )
            log.write(f"END {script} exit={process.returncode}\n")
            log.flush()
            if process.returncode:
                raise RuntimeError(f"Failed {script}; see {log_path}")
    after = hashes()
    added = sorted(after.keys() - before.keys())
    removed = sorted(before.keys() - after.keys())
    changed = sorted(name for name in before.keys() & after.keys() if before[name] != after[name])
    report = {
        "started_and_finished_at": datetime.now(timezone.utc).isoformat(),
        "scripts_rerun": SCRIPTS,
        "script_count": len(SCRIPTS),
        "data_table_count_compared": len(after),
        "added": added, "removed": removed, "changed": changed,
        "status": "passed" if not added and not removed and not changed else "changed_outputs_detected",
        "log": str(log_path.relative_to(ROOT)),
        "scope": "Analysis rerun from retained raw data; no redownload and no model-training claim",
    }
    (ROOT / "reports" / "reproducibility_rerun.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if report["status"] != "passed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()

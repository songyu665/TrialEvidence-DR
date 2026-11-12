# TrialEvidence-DR

Reproducible code for auditing clinical-trial evidence labels, drug–disease mapping, selection bias, historical recoverability, and later-evidence ranking.

Public repository: https://github.com/songyu665/TrialEvidence-DR

Current source release: https://github.com/songyu665/TrialEvidence-DR/releases/tag/v0.1.1

## Requirements

- Python 3.12
- Linux recommended for exact reproduction of the recorded environment

## Setup

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r configs/requirements-freeze.txt
export TRIALEVIDENCE_ROOT=/absolute/path/to/project
```

Use an active virtual environment for all Python execution. Keep downloaded data and generated outputs inside the project root. The code does not require server credentials or any Kimi-specific tooling.

## Data layout

- `data/raw/`: source downloads and original checksums
- `data/interim/`: mapped and adjudicated tables
- `data/processed/`: frozen analysis tables
- `reports/`: audit reports, tables and figures
- `runs/`: logs and pause/resume checkpoints

## Reproduction

After restoring the documented raw data layout:

```bash
.venv/bin/python scripts/reproduce_analysis.py
```

The script reruns the main analyses without redownloading and compares table hashes. It writes `reports/reproducibility_rerun.json` and `runs/reproducibility_rerun.log`; these generated outputs are not included in this code-only repository.

## Evidence boundaries

Rule-based labels are not clinical gold standards. Unknown pairs are not true negatives. Evidence-status models do not predict therapeutic efficacy. Historical-version and resampling analyses are not independent external validation.

AI-assisted semantic screening is part of the documented research workflow where applicable; code cleaning does not alter that methodological disclosure.

import os
from pathlib import Path

PROJECT_ROOT = Path(os.environ.get("TRIALEVIDENCE_ROOT", Path(__file__).resolve().parents[2])).resolve()
ENV_PREFIX = PROJECT_ROOT / "venv_trialevidence_dr"


def project_path(*parts: str) -> Path:
    path = PROJECT_ROOT.joinpath(*parts).resolve()
    if path != PROJECT_ROOT and PROJECT_ROOT not in path.parents:
        raise ValueError(f"Path escapes the project root: {path}")
    return path


def require_project_path(path: str | Path) -> Path:
    resolved = Path(path).resolve()
    if resolved != PROJECT_ROOT and PROJECT_ROOT not in resolved.parents:
        raise ValueError(f"Path escapes the project root: {resolved}")
    return resolved

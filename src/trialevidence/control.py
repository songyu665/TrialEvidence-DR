import json
import os
import re
import signal
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .paths import project_path

RUN_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as stream:
            json.dump(payload, stream, ensure_ascii=False, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        if os.path.exists(temporary):
            os.unlink(temporary)
        raise


class RunController:
    def __init__(self, run_name: str) -> None:
        if not RUN_NAME.fullmatch(run_name):
            raise ValueError("Invalid run name")
        self.run_name = run_name
        self.run_dir = project_path("runs", run_name)
        self.run_dir.mkdir(parents=True, exist_ok=True)
        self.pause_file = self.run_dir / "PAUSE_REQUESTED"
        self.stop_file = self.run_dir / "STOP_REQUESTED"
        self.checkpoint_file = self.run_dir / "checkpoint.json"
        self._signal_received: str | None = None

    def install_signal_handlers(self) -> None:
        def handler(signum: int, _frame: Any) -> None:
            self._signal_received = signal.Signals(signum).name

        signal.signal(signal.SIGINT, handler)
        signal.signal(signal.SIGTERM, handler)

    @property
    def stop_reason(self) -> str | None:
        if self._signal_received:
            return self._signal_received
        if self.stop_file.exists():
            return "STOP_REQUESTED"
        if self.pause_file.exists():
            return "PAUSE_REQUESTED"
        return None

    def request_pause(self) -> None:
        _atomic_json(self.pause_file, {"requested_at": self._now()})

    def request_stop(self) -> None:
        _atomic_json(self.stop_file, {"requested_at": self._now()})

    def clear_requests(self) -> None:
        for path in (self.pause_file, self.stop_file):
            if path.exists():
                path.unlink()

    def save_checkpoint(self, state: dict[str, Any]) -> None:
        _atomic_json(
            self.checkpoint_file,
            {
                "format_version": 1,
                "run_name": self.run_name,
                "saved_at": self._now(),
                "state": state,
            },
        )

    def load_checkpoint(self) -> dict[str, Any]:
        if not self.checkpoint_file.exists():
            raise FileNotFoundError(f"No checkpoint for run {self.run_name}")
        with self.checkpoint_file.open(encoding="utf-8") as stream:
            payload = json.load(stream)
        if payload.get("format_version") != 1 or payload.get("run_name") != self.run_name:
            raise ValueError("Checkpoint is incompatible with this run")
        state = payload.get("state")
        if not isinstance(state, dict):
            raise ValueError("Checkpoint state must be an object")
        return state

    def mark_completed(self, summary: dict[str, Any]) -> None:
        _atomic_json(
            self.run_dir / "COMPLETED.json",
            {"completed_at": self._now(), "run_name": self.run_name, "summary": summary},
        )
        self.clear_requests()

    @staticmethod
    def _now() -> str:
        return datetime.now(timezone.utc).isoformat()

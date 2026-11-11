import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.control import RunController
from trialevidence.paths import ENV_PREFIX, PROJECT_ROOT, project_path


def main() -> None:
    assert PROJECT_ROOT == ROOT
    assert Path(sys.prefix).resolve() != Path(sys.base_prefix).resolve()
    assert project_path("runs", "probe").is_relative_to(PROJECT_ROOT)
    try:
        project_path("..", "outside")
    except ValueError:
        pass
    else:
        raise AssertionError("Project path guard failed")

    controller = RunController("runtime_control_smoke")
    controller.clear_requests()
    controller.save_checkpoint({"epoch": 3, "status": "before_pause"})
    controller.request_pause()
    assert controller.stop_reason == "PAUSE_REQUESTED"
    state = controller.load_checkpoint()
    assert state["epoch"] == 3
    controller.clear_requests()
    state["status"] = "resumed"
    controller.save_checkpoint(state)
    assert controller.stop_reason is None
    assert controller.load_checkpoint()["status"] == "resumed"
    controller.mark_completed({"verification": "passed"})
    print(json.dumps({"environment": sys.prefix, "project_root": str(PROJECT_ROOT), "pause_resume": "passed"}, indent=2))


if __name__ == "__main__":
    main()

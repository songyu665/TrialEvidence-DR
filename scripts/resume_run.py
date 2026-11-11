import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from trialevidence.control import RunController


def main() -> None:
    parser = argparse.ArgumentParser(description="Clear pause requests and verify the latest checkpoint")
    parser.add_argument("run_name")
    args = parser.parse_args()
    controller = RunController(args.run_name)
    state = controller.load_checkpoint()
    controller.clear_requests()
    print(json.dumps({"run_name": args.run_name, "checkpoint_state": state}, ensure_ascii=False, indent=2))
    print("Restart the original experiment command with --resume.")


if __name__ == "__main__":
    main()

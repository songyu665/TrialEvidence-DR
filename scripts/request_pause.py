import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from trialevidence.control import RunController


def main() -> None:
    parser = argparse.ArgumentParser(description="Request a checkpointed pause for a run")
    parser.add_argument("run_name")
    parser.add_argument("--stop", action="store_true", help="Request termination instead of a resumable pause")
    args = parser.parse_args()
    controller = RunController(args.run_name)
    if args.stop:
        controller.request_stop()
    else:
        controller.request_pause()


if __name__ == "__main__":
    main()

"""Run each complete ACCESS member/scenario as its downloaded files become ready.

Launch fetch_access_wpsh_winds.py separately first. Default timeout is two hours.
Completed data are revalidated by the scientific pipeline; no partial years run.
"""
import argparse
from pathlib import Path
import subprocess
import sys
import time

from wpsh_streamfunction import FUTURE_PERIOD, HISTORICAL_PERIOD

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/wpsh_winds/access"


def ready(member, scenario):
    for variable in ("ua", "va"):
        for experiment, period in (("historical", HISTORICAL_PERIOD), (scenario, FUTURE_PERIOD)):
            paths = [p for p in DATA.glob(f"{variable}850_{experiment}_{member}_*.nc") if ".part." not in p.name]
            years = []
            for path in paths:
                first, last = map(int, path.stem.split("_")[-2:])
                years.extend(range(first, last+1))
            if sorted(years) != list(range(period[0], period[1]+1)):
                return False
    return True


def main(timeout):
    pending = [(s, m) for s in ("ssp245", "ssp585") for m in ("r1", "r2", "r5")]
    deadline = time.monotonic()+timeout
    while pending:
        ran = False
        for scenario, member in pending.copy():
            if ready(member, scenario):
                for method in ("harmonic", "nearest"):
                    print(f"RUN {scenario} {member} {method}", flush=True)
                    subprocess.run([sys.executable, "-u", "scripts/task7_wpsh_eof.py", "--source", "access",
                                    "--members", member, "--scenarios", scenario, "--wind-fill", method], cwd=ROOT, check=True)
                pending.remove((scenario, member))
                ran = True
        if pending and not ran:
            if time.monotonic() > deadline:
                raise TimeoutError(f"Incomplete downloads for: {pending}")
            print(f"Waiting for complete winds: {pending}", flush=True)
            time.sleep(30)
    subprocess.run([sys.executable, "scripts/summarize_access_wpsh.py"], cwd=ROOT, check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout", type=float, default=7200)
    main(parser.parse_args().timeout)

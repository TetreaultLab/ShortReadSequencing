import argparse
from datetime import datetime
import subprocess
import os
from pathlib import Path

parser = argparse.ArgumentParser(
    prog="Filter snpEff output.",
    description="Filter snpEff output.",
)

parser.add_argument("--toml_config", type=str, required=True, help="Toml config file.")
parser.add_argument("--start", type=str, required=True, help="Start time.")
parser.add_argument("--sample", type=str, required=True, help="Sample name.")

args = parser.parse_args()

toml_config = args.toml_config
start = args.start
sample = args.sample


def get_time():
    now = datetime.now()
    return now


# Move results to projects directory
print("\n>>> Transferring from scratch to results")
current_directory = os.getcwd()
results = Path(str(current_directory).replace("/work", "/results/"))
output = Path(toml_config["general"]["output"] + "/" + sample)
results.parent.mkdir(parents=True, exist_ok=True)
cmd = ["rsync", "-avxH", "--no-g", "--no-p", "--partial", str(output), str(results)]
subprocess.run(cmd, check=True)
final = f"{str(results)}/{sample}"
print(f"\n\n>>> Results found in {final}")

end = get_time()
total_time = end - start
end_str = ">>> {}-seq pipeline for {} completed in {}.".format(
    toml_config["general"]["sequencing"], sample, total_time
)
print("=" * len(end_str) + "\n" + end_str + "\n" + "=" * len(end_str))

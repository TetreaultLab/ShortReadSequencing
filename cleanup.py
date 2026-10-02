import argparse
from datetime import datetime
import subprocess
import os
from pathlib import Path

parser = argparse.ArgumentParser(
    prog="Clean up SRS output.",
    description="TRansfering files from scratch to results in projects.",
)

parser.add_argument("--output", type=str, required=True, help="Output path.")
parser.add_argument("--sample", type=str, required=True, help="Sample name.")

args = parser.parse_args()

output = args.output
sample = args.sample

# Move results to projects directory
print("\n>>> Transferring from scratch to results")
current_directory = os.getcwd()
results = Path(str(current_directory).replace("/work", "/results/"))
output = Path(output + "/" + sample)
results.parent.mkdir(parents=True, exist_ok=True)
cmd = ["rsync", "-avxH", "--no-g", "--no-p", "--partial", str(output), str(results)]
subprocess.run(cmd, check=True)
final = f"{str(results)}/{sample}"
print(f"\n\n>>> Results found in {final}")

end_str = ">>> Pipeline for {} completed.".format(sample)
print("=" * len(end_str) + "\n" + end_str + "\n" + "=" * len(end_str))

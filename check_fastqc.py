import argparse
from pathlib import Path
import zipfile
import re

parser = argparse.ArgumentParser(
    prog="GenerateJob",
    description="Generate SLURM job for samples based on a TOML project config file.",
)

parser.add_argument(
    "--path", type=str, required=True, help="Path for fastqc report zip."
)
parser.add_argument("--sample", type=str, required=True, help="Sample name.")

args = parser.parse_args()

output = args.path
sample = args.sample

fastqc_reports = list(Path(output + "/").glob(sample + "*_fastqc.zip"))


for zip_path in fastqc_reports:
    with zipfile.ZipFile(zip_path, "r") as z:
        file_list = z.namelist()
        summary_file = [f for f in file_list if f.endswith("summary.txt")][0]
        try:
            with z.open(summary_file) as f:
                content = f.read().decode("utf-8")
                print(content.rstrip())
        except KeyError:
            print(f"No summary.txt found in {zip_path}")

        # Search for WARN or FAIL
        if re.search(r"\b(FAIL)\b", content):
            print(
                ">>> FastQC failures detected.\n>>> Please check FastQC report (and adjust for trimming if necessary) before resubmitting (with --redo).\n"
            )
            # sys.exit(1)
        elif re.search(r"\b(WARN)\b", content):
            print(">>> FastQC warnings detected\n")
        else:
            print(">>> FastQC passed with no warnings or failures.\n")

import argparse
import toml
import os
import subprocess

parser = argparse.ArgumentParser(
    prog="GenerateJob",
    description="Generate SLURM job for samples based on a TOML project config file.",
)

parser.add_argument("--sample", type=str, required=True, help="Sample name.")
parser.add_argument(
    "--config", type=str, required=True, help="Project config file, including path."
)
parser.add_argument(
    "--launch",
    action="store_true",
    required=False,
    help="Decide if it launches the job for sample or not. True will automatically launch the slurm file. False will only create the file",
)

args = parser.parse_args()

work_dir = os.getcwd()

# SAMPLE NAME
sample_name = args.sample

# LOAD CONFIG FILE
toml_file = args.config

with open(toml_file, "r") as f:
    toml_config = toml.load(f)

email = toml_config["general"]["email"]
time = toml_config["general"]["time"]
project = toml_config["general"]["project"]

f.close()

path_config = work_dir + "/" + toml_file

slurm = """#!/bin/sh
#SBATCH -N 1 #Number of nodes
#SBATCH --cpus-per-task 8 # number of cores
#SBATCH --mem 64G # memory pool for all cores
#SBATCH -t {1} # time (DD-HH:MM)
#SBATCH -o {2}_{3}.%N.%j.log
#SBATCH -e {2}_{3}.%N.%j.log
#SBATCH --mail-type=FAIL
#SBATCH --account=rrg-tetreaum
#SBATCH --mail-user={4}
#
### Load environnment
#
module load StdEnv/2023
module load python/3.11.5
module load fastqc/0.12.1 bbmap/39.06 star/2.7.11a bwa-mem2/2.2.1 gcc/12.3 openmpi/4.1.5 salmon/1.10.2 samtools/1.18 gatk/4.4.0.0 subread/2.0.6 bcftools/1.18 freebayes/1.3.7
source /lustre03/project/6019267/shared/tools/PIPELINES/ShortReadSequencing/bin/activate
#
### Launch script
#
# newgrp rrg-tetreaum
#
python -u /lustre03/project/6019267/shared/tools/PIPELINES/ShortReadSequencing/ShortReadSequencing/pipeline_short.py --sample {2} --config {5}
#
""".format(
    time, sample_name, project, email, path_config
)

print(slurm, file=open(work_dir + "/" + sample_name + "_" + project + ".slurm", "w"))

if args.launch is True:
    subprocess.run(["sbatch", work_dir + "/" + sample_name + "_" + project + ".slurm"])

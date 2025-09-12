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
cores = toml_config["general"]["threads"]
memory = toml_config["general"]["memory"]
time = toml_config["general"]["time"]
project = toml_config["general"]["project"]

f.close()

path_config = work_dir + "/" + toml_file

slurm = """#!/bin/sh
#SBATCH -N 1 #Number of nodes
#SBATCH --cpus-per-task {0} # number of cores
#SBATCH --mem {1}G # memory pool for all cores
#SBATCH -t {2} # time (DD-HH:MM)
#SBATCH -o {3}_{4}.%N.%j.log
#SBATCH -e {3}_{4}.%N.%j.log
#SBATCH --mail-type=FAIL
#SBATCH --account=rrg-tetreaum
#SBATCH --mail-user={5}
#
### Load environnment
#
module load StdEnv/2023
module load python/3.10.13
module load fastqc/0.12.1 bbmap/39.06 star/2.7.11a bwa-mem2/2.2.1 gcc/12.3 openmpi/4.1.5 salmon/1.10.2 samtools/1.22.1 gatk/4.6.1.0 subread/2.0.6 htslib/1.22.1 bcftools/1.22 freebayes/1.3.7
source /lustre09/project/6019267/shared/tools/main_pipelines/long-read/launch_pipeline_env/bin/activate
#
### Launch script
#
# newgrp rrg-tetreaum
#
python -u /lustre09/project/6019267/shared/tools/main_pipelines/short-read/ShortReadSequencing/pipeline_short.py --sample {3} --config {6}
#
""".format(
    cores, memory, time, sample_name, project, email, path_config
)

print(slurm, file=open(work_dir + "/" + sample_name + "_" + project + ".slurm", "w"))

if args.launch is True:
    subprocess.run(["sbatch", work_dir + "/" + sample_name + "_" + project + ".slurm"])

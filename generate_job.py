import argparse
import toml
import os

parser = argparse.ArgumentParser(
    prog="GenerateJob",
    description="Generate SLURM job for samples based on a TOML project config file.",
)

parser.add_argument("--sample", type=str, required=True, help="Sample name.")
parser.add_argument(
    "--config", type=str, required=True, help="Project config file, including path."
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

# time
if toml_config["general"]["sequencing"] == "RNA":
    time = "00-23:59"
elif toml_config["general"]["sequencing"] == "Exome":
    time = "00-15:59"
elif toml_config["general"]["sequencing"] == "Genome":
    time = "01-11:59"

f.close()

path_config = work_dir + "/" + toml_file

# slurm = """#!/bin/sh
# #SBATCH -N 1 #Number of nodes
# #SBATCH -n {0} # number of cores
# #SBATCH --mem {1}G # memory pool for all cores
# #SBATCH -t {2} # time (DD-HH:MM)
# #SBATCH -o {3}.%N.%j.log
# #SBATCH -e {4}.%N.%j.log
# #SBATCH --mail-type=FAIL
# #SBATCH --account=rrg-tetreaum
# #SBATCH --mail-user={5}
# #
# ### Load environnment
# #
# module load python/3.10.2 StdEnv/2020 fastqc/0.11.9 bbmap/38.86 star/2.7.9a bwa-mem2/2.2.1 gcc/9.3.0 openmpi/4.0.3 salmon/1.4.0 subread/2.0.3 samtools/1.17
# source /lustre03/project/6019267/shared/tools/PIPELINES/ShortReadSequencing/bin/activate
# #
# ### Variables
# #
# sample='{6}'
# config='{7}'
# #
# ### Launch script
# #
# python /lustre03/project/6019267/shared/tools/PIPELINES/shortReads/pipeline_short_v202402.py --sample ${sample} --config ${config}
# #
# """.format(  # noqa: F524
#     cores, memory, time, sample_name, sample_name, email, sample_name, path_config
# )

slurm = """#!/bin/sh
#
### Load environnment
#
module load python/3.10.2 StdEnv/2020 fastqc/0.11.9 bbmap/38.86 star/2.7.9a bwa-mem2/2.2.1 gcc/9.3.0 openmpi/4.0.3 salmon/1.4.0 subread/2.0.3 samtools/1.17
source /lustre03/project/6019267/shared/tools/PIPELINES/ShortReadSequencing/bin/activate
#
### Launch script
#
python /lustre04/scratch/mlab/pipeline2024/ShortReadSequencing/pipeline_short_v202402.py --sample {0} --config {1}
#
""".format(sample_name, path_config)  # noqa: F524

print(slurm)
print(slurm, file=open(work_dir + "/" + sample_name + ".slurm", "w"))

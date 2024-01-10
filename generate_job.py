import argparse
import toml
import os

parser = argparse.ArgumentParser(
                    prog='GenerateJob',
                    description='Generate SLURM job for samples based on a TOML project config file.')

parser.add_argument('--sample', type=str, required=True, help='Sample name.')
parser.add_argument('--config', type=str, required=True, help='Project config file, including path.')

args = parser.parse_args()

work_dir = os.getcwd()

# SAMPLE NAME
sample_name = args.sample

# LOAD CONFIG FILE
toml_file = args.config

with open(toml_file, 'r') as f:
    toml_config = toml.load(f)

email = toml_config['general']['email']
cores = toml_config['general']['threads']
memory = toml_config['general']['memory']

# time
if toml_config["general"]["sequencing"] == "RNA":
    time = "00-23:59"
elif toml_config["general"]["sequencing"] == "Exome":
    time = "00-15:59"
elif toml_config["general"]["sequencing"] == "Genome":
    time = "01-11:59"

f.close()

path_config = work_dir + "/" + toml_file

slurm = '''#!/bin/sh
#SBATCH -N 1 #Number of nodes
#SBATCH -n {0} # number of cores
#SBATCH --mem {1}G # memory pool for all cores
#SBATCH -t {2} # time (DD-HH:MM)
#SBATCH -o {3}.%N.%j.log
#SBATCH -e {3}.%N.%j.log
#SBATCH --mail-type=FAIL
#SBATCH --account=rrg-tetreaum
#SBATCH --mail-user={4}
#
### Load environnment
#
module load python/3.10.2 StdEnv/2020 fastqc/0.11.9 bbmap/38.86 star/2.7.9a bwa/0.7.17 gcc/9.3.0 openmpi/4.0.3 salmon/1.4.0 samtools/1.17
source /lustre03/project/6019267/shared/tools/PIPELINES/ShortReadSequencing/bin/activate
#
### Variables
#
sample='{3}'
config='{5}'
#
### Launch script
#
python /lustre03/project/6019267/shared/tools/PIPELINES/shortReads/pipeline_short_v202402.py --sample ${sample} --config ${config}
#
'''.format(cores, memory, time, sample_name, email, path_config)

print(slurm)
print(slurm,  file=open(work_dir+"/"+sample_name+".slurm", 'w'))
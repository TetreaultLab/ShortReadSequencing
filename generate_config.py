import argparse
import os

# Arguments
parser = argparse.ArgumentParser(
    prog="GenerateConfig",
    description="Generate TOML config file based on tools selected.",
)

parser.add_argument("--project", type=str, required=True, help="Project name")
parser.add_argument(
    "--trimming",
    type=str,
    required=True,
    help='choose trimming tool. Options: "none", "bbduk".',
)
parser.add_argument(
    "--alignment",
    type=str,
    required=True,
    help='choose alignment tool. Options: "none", "star", "bwa".',
)
parser.add_argument(
    "--pseudo",
    type=str,
    required=True,
    help='choose pseudoalignment tool. Options: "none", "salmon".',
)
parser.add_argument(
    "--quantification",
    type=str,
    required=True,
    help='choose quantification tool. Options: "none", "featurecounts".',
)

args = parser.parse_args()

work_dir = os.getcwd()
print(work_dir)

# Project NAME
project_name = args.project

# TRIMMING
if args.trimming == "none":
    t = "# No trimming\n"
elif args.trimming == "bbduk":
    t = """# Trimming
[bbduk]
    ordered = "f"   # Set to true to output reads in same order as input.
    kmers = 27  # Kmer length used for finding contaminants.  Contaminants shorter than k will not be found.  k must be at least 1.
    qtrim = "r"    # Trim read ends to remove bases with quality below trimq. Performed AFTER looking for kmers.  Values: rl (trim both ends), f (neither end), r (right end only), l (left end only), w (sliding window).
    trimq = 6 # Regions with average quality BELOW this will be trimmed, if qtrim is set to something other than f.
    minlength = 10 # Reads shorter than this after trimming will be discarded.  Pairs will be discarded if both are shorter.
    mlf = 0 # (minlengthfraction) Reads shorter than this fraction of original length after trimming will be discarded.
    minavgquality = 0   # (maq) Reads with average quality (after trimming) below this will be discarded.
"""

else:
    print(
        f"'{args.trimming}' is not a valid trimming option. Check generate_config .py --help for options."
    )
    exit()

# ALIGNMENT
if args.alignment == "none":
    a = "# No alignment\n"
elif args.alignment == "star":
    a = """# Alignment
[star]
    outSAMtype = "BAM SortedByCoordinate"   # output file type
    quantMode  = "GeneCounts"   # how it counts Gene vs Transcript

"""
elif args.alignment == "bwa":
    a = """# Alignment
[bwa-mem]

"""
else:
    print(
        f"'{args.alignment}' is not a valid alignment option. Check generate_config .py --help for options."
    )
    exit()

# PSEUDOALIGNMENT
if args.pseudo == "none":
    p = "# No pseudoalignment\n"
elif args.pseudo == "salmon":
    p = """# Pseudoalignment
[salmon]

"""
else:
    print(
        f"'{args.pseudo}' is not a valid pseudoalignment option. Check generate_config .py --help for options."
    )
    exit()

# QUANTIFICATION
if args.quantification == "none":
    q = "# No quantification\n"
elif args.quantification == "featurecounts":
    q = """# Quantification
[featurecounts]

"""
else:
    print(
        f"'{args.quantification}' is not a valid quantification option. Check generate_config .py --help for options."
    )
    exit()


# Create config file
general = """# TOML config file for {}

# Required information
[general]
    project     = "{}"   # Choose a project name. Should not start with a number.
    fastq      = "/path/to/raw/data"   # Path to raw fastq.gz files. Ex: /lustre03/project/6019267/shared/data/<project>.
    output     = "/path/to/output" # Preferably use your scratch. The directory will be created. Exemple: /lustre04/scratch/<user>/<project>/output/.
    temporary  = "/path/to/output/tmp" # Preferably use your scratch. The directory will be created. Exemple: /lustre04/scratch/<user>/<project>/tmp/. 
    sequencing = "" # Type of sequencing. Short read RNA or DNA. Possible values: ["RNA", "Exome", "Genome"].
    reads      = "" # Type of reads sequencing. Either single-end or paired-end. Possible values: ["SE", "PE"].
    reference  = "" # Possible values: Human: ["hg19", "hg38"]. Mouse: ["mm39"]. Worm: ["ce11"]. Zebrafish: ["danRer11"].
    threads    = 8  # Number of threads.
    memory     = 32 # Total memory needed.
    email      = "" # You e-mail adress to receive notification

# Quality control
[fastqc]
    kmers = 7   # Specifies the length of Kmer to look for in the Kmer content module. Specified Kmer length must be between 2 and 10.

[multiqc]

# Sorting and indexing
[samtools]

# MarkDuplicates
[markduplicates]

""".format(project_name, project_name)

toml_config = general + "\n" + t + "\n" + a + "\n" + p + "\n" + q
print(toml_config)

print(toml_config, file=open(work_dir + "/" + project_name + ".config.toml", "w"))
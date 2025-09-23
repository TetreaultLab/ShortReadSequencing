import argparse
import os

# Arguments
parser = argparse.ArgumentParser(
    prog="GenerateConfig",
    description="Generate TOML config file based on tools selected.",
)

parser.add_argument("--project", type=str, required=True, help="Project name")

# Choose only rna, exome or genome (will help determine the time needed)
group = parser.add_mutually_exclusive_group(required=True)

group.add_argument("--rna", action="store_true", help="Run RNA-seq analysis")
group.add_argument("--exome", action="store_true", help="Run Exome-seq analysis")
group.add_argument("--genome", action="store_true", help="Run WGS analysis")

parser.add_argument(
    "--trimming",
    action='store_true',
    help='Enable trimming with bbduk.',
)
parser.add_argument(
    "--pseudo",
    action='store_true',
    help='Enable pseudoalignment with salmon.',
)
parser.add_argument(
    "--quantification",
    action='store_true',
    help='Enable gene quantification with featurecounts.',
)
parser.add_argument(
    "--variants",
    action='store_true',
    help='Enable variant calling.',
)

args = parser.parse_args()

# if nothing is passed, defaults to trim, quantification and variant calling
if not any([args.trimming, args.pseudo, args.quantification, args.variants]):
    args.trim = True
    args.quantification = True
    args.variants = True

work_dir = os.getcwd()

# Project NAME
project_name = args.project

options = ""

# TRIMMING
if args.trimming:
    options += """# Trimming
[bbduk]
    ordered = "f"   # Set to true to output reads in same order as input.
    kmers = 8  # Kmer length used for finding contaminants.  Contaminants shorter than k will not be found.  k must be at least 1.
    mink = 7  # look for shorter kmers at read tips to this min 
    ktrim = "r"  # trim bases that match adapters, trim to the right
    qtrim = "rl"    # Trim read ends to remove bases with quality below trimq. Performed AFTER looking for kmers.  Values: rl (trim both ends), f (neither end), r (right end only), l (left end only), w (sliding window).
    trimq = 10 # Regions with average quality BELOW this will be trimmed, if qtrim is set to something other than f.
    minlength = 10 # Reads shorter than this after trimming will be discarded.  Pairs will be discarded if both are shorter.
    mlf = 0 # (minlengthfraction) Reads shorter than this fraction of original length after trimming will be discarded.
    minavgquality = 0   # (maq) Reads with average quality (after trimming) below this will be discarded.
    \n
"""
else:
    options += "# No trimming\n"

# ALIGNMENT
if args.rna:
    sequencing = "rna"
    alignment = "star"
    options += """# Alignment
[star]
    outSAMtype1 = "BAM"   # type of SAM/BAM output
    outSAMtype2 = "SortedByCoordinate"
    twopassMode = "Basic"    # 2-pass mapping mode. None or Basic
    outSJtype = "Standard"  # type of splice junction output
    quantMode  = "GeneCounts"   # types of quantification requested. -, TranscriptomeSAM and/or GeneCounts
    \n
"""
else: 
    if args.exome:
        sequencing = "exome"
    elif args.genome:
        sequencing = "genome"
    
    alignment = "bwa"
    options += """# Alignment
[bwa-mem]
\n
"""

# PSEUDOALIGNMENT
if args.pseudo:
    options += """# Pseudoalignment
[salmon]
    minScoreFraction = 0.65	# The fraction of the optimal possible alignment score that a mapping must achieve in order to be considered "valid" --- should be in [0,1]. Salmon Default 0.65 and Alevin Default 0.87.
    \n
"""
else:
    options += "# No pseudoalignment\n"
# QUANTIFICATION
if args.quantification:
    options += """# Quantification
[featurecounts]
    features = "gene"   # Specify feature type(s) in a GTF annotation. If multiple types are provided, they should be separated by ',' with no space in between. 'exon' by default. Rows in the annotation with a matched feature will be extracted and used for read mapping.
    attribute = "gene_id"   # Specify attribute type in GTF annotation. 'gene_id' by default. Meta-features used for read counting will be extracted from annotation using the provided value.
    overlap = 1 # Minimum number of overlapping bases in a read that is required for read assignment. 1 by default. Number of overlapping bases is counted from both reads if paired end.
    \n
"""
else:
    options += "# No quantification\n"

# VARIANT CALLING
if args.variants: 
    options += """# Variant calling
[snpeff]
    \n
"""
    
else: options += "# No variant calling\n"
    

# Create config file
general = """# TOML config file for {0}

# Required information
[general]
    threads = 8  # Number of threads.
    memory = 64 # Total memory needed.
    time = "02-23:59" # DD-HH:MM. Default: 2 days 23h59
    project = "{0}"   # Choose a project name. Should not start with a number.
    output = "/lustre10/scratch/$USER/{0}/output" # Path to your output, preferably use your scratch. The directory will be created. Exemple: /lustre10/scratch/<user>/<project>/output/.
    temporary = "/lustre10/scratch/$USER/{0}/tmp" # The directory will be created. Exemple: /lustre10/scratch/<user>/<project>/tmp/. 
    sequencing = "{1}" # Type of sequencing. Short read RNA or DNA.
    trimming = "{2}"
    alignment = "{3}"
    pseudo = "{4}"
    quantification = "{5}"
    variants = "{6}"
    # To fill
    fastq = "" # Path to raw fastq.gz files. Ex: /lustre09/project/6019267/shared/data/<project>.
    reads = "" # Type of reads sequencing. Either single-end or paired-end. Possible values: ["SE", "PE"].
    reference = "" # Possible values: Human: ["grch37", "grch38"]. Mouse: ["grcm39"]. Worm: ["wbcel235"]. Zebrafish: ["grcz11"].
    email = "" # You e-mail address to receive notification
    

# Quality control
[fastqc]
    kmers = 7   # Specifies the length of Kmer to look for in the Kmer content module. Specified Kmer length must be between 2 and 10.

[multiqc]

# Sorting and indexing
[samtools]

# MarkDuplicates
[markduplicates]
    index = "true"    # Whether to create an index when writing VCF or coordinate sorted BAM output.
    strategy = "SUM_OF_BASE_QUALITIES"   # The scoring strategy for choosing the non-duplicate among candidates.  Default value: SUM_OF_BASE_QUALITIES. Possible values: [SUM_OF_BASE_QUALITIES, TOTAL_MAPPED_REFERENCE_LENGTH, RANDOM].
    remove = "false"
\n
""".format(
    project_name,
    sequencing,
    args.trimming,
    alignment,
    args.pseudo,
    args.quantification,
    args.variants,
)  # noqa: F524

toml_config = general + options

print(toml_config, file=open(work_dir + "/" + project_name + ".config.toml", "w"))

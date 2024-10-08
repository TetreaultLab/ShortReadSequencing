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
parser.add_argument(
    "--variant",
    type=str,
    required=True,
    help='choose variant tool. Options: "none", "yes".',
)


args = parser.parse_args()

work_dir = os.getcwd()

# Project NAME
project_name = args.project

options = ""
# TRIMMING
if args.trimming == "none":
    options += "# No trimming\n"
elif args.trimming == "bbduk":
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
    print(
        f"'{args.trimming}' is not a valid trimming option. Check generate_config.py --help for options."
    )
    exit()

# ALIGNMENT
if args.alignment == "none":
    options += "# No alignment\n"
elif args.alignment == "star":
    options += """# Alignment
[star]
    outSAMtype1 = "BAM"   # type of SAM/BAM output
    outSAMtype2 = "SortedByCoordinate"
    twopassMode = "Basic"    # 2-pass mapping mode. None or Basic
    outSJtype = "Standard"  # type of splice junction output
    quantMode  = "GeneCounts"   # types of quantification requested. -, TranscriptomeSAM and/or GeneCounts
    \n
"""
elif args.alignment == "bwa":
    options += """# Alignment
[bwa-mem]
\n
"""
else:
    print(
        f"'{args.alignment}' is not a valid alignment option. Check generate_config.py --help for options."
    )
    exit()

# PSEUDOALIGNMENT
if args.pseudo == "none":
    options += "# No pseudoalignment\n"
elif args.pseudo == "salmon":
    options += """# Pseudoalignment
[salmon]
    minScoreFraction = 0.65	# The fraction of the optimal possible alignment score that a mapping must achieve in order to be considered "valid" --- should be in [0,1]. Salmon Default 0.65 and Alevin Default 0.87.
    \n
"""
else:
    print(
        f"'{args.pseudo}' is not a valid pseudoalignment option. Check generate_config.py --help for options."
    )
    exit()

# QUANTIFICATION
if args.quantification == "none":
    options += "# No quantification\n"
elif args.quantification == "featurecounts":
    options += """# Quantification
[featurecounts]
    features = "gene"   # Specify feature type(s) in a GTF annotation. If multiple types are provided, they should be separated by ',' with no space in between. 'exon' by default. Rows in the annotation with a matched feature will be extracted and used for read mapping.
    attribute = "gene_id"   # Specify attribute type in GTF annotation. 'gene_id' by default. Meta-features used for read counting will be extracted from annotation using the provided value.
    overlap = 1 # Minimum number of overlapping bases in a read that is required for read assignment. 1 by default. Number of overlapping bases is counted from both reads if paired end.
    \n
"""
else:
    print(
        f"'{args.quantification}' is not a valid quantification option. Check generate_config.py --help for options."
    )
    exit()

# VARIANT CALLING
if args.variant == "none":
    options += "# No variant calling\n"
elif args.variant == "yes":
    options += """# Variant calling
[snpeff]
    \n
"""

# Create config file
general = """# TOML config file for {0}

# Required information
[general]
    project = "{0}"   # Choose a project name. Should not start with a number.
    fastq = "" # Path to raw fastq.gz files. Ex: /lustre03/project/6019267/shared/data/<project>.
    output = "/lustre04/scratch/<USER>/{0}/output" # # Path to your output, preferably use your scratch. The directory will be created. Exemple: /lustre04/scratch/<user>/<project>/output/.
    temporary = "/lustre04/scratch/<USER>/{0}/tmp" # The directory will be created. Exemple: /lustre04/scratch/<user>/<project>/tmp/. 
    sequencing = "" # Type of sequencing. Short read RNA or DNA. Possible values: ["RNA", "Exome", "Genome"].
    reads = "" # Type of reads sequencing. Either single-end or paired-end. Possible values: ["SE", "PE"].
    reference = "" # Possible values: Human: ["grch37", "grch38"]. Mouse: ["grcm39"]. Worm: ["wbcel235"]. Zebrafish: ["grcz11"].
    trimming = "{1}"
    alignment = "{2}"
    pseudo = "{3}"
    quantification = "{4}"
    variant = "{5}"
    threads = 8  # Number of threads.
    memory = 64 # Total memory needed.
    time = "02-23:59" # DD-HH:MM. Default: 2 days 23h59
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
\n
""".format(
    project_name,
    args.trimming,
    args.alignment,
    args.pseudo,
    args.quantification,
    args.variant,
)  # noqa: F524

toml_config = general + options

print(toml_config, file=open(work_dir + "/" + project_name + ".config.toml", "w"))

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
    help='choose variant tool. Options: "none", "snpeff".',
)


args = parser.parse_args()

work_dir = os.getcwd()

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
    outSAMtype1 = "BAM"   # type of SAM/BAM output
    outSAMtype2 = "SortedByCoordinate"
    twopassMode = "None"    # 2-pass mapping mode. None or Basic
    outSJtype = "Standard"  # type of splice junction output
    quantMode  = "GeneCounts"   # types of quantification requested. -, TranscriptomeSAM and/or GeneCounts
"""
elif args.alignment == "bwa":
    a = """# Alignment
[bwa-mem]

"""
else:
    print(
        f"'{args.alignment}' is not a valid alignment option. Check generate_config.py --help for options."
    )
    exit()

# PSEUDOALIGNMENT
if args.pseudo == "none":
    p = "# No pseudoalignment\n"
elif args.pseudo == "salmon":
    p = """# Pseudoalignment
[salmon]
    minScoreFraction = 0.65	# The fraction of the optimal possible alignment score that a mapping must achieve in order to be considered "valid" --- should be in [0,1]. Salmon Default 0.65 and Alevin Default 0.87.
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
    features = "gene"   # Specify feature type(s) in a GTF annotation. If multiple types are provided, they should be separated by ',' with no space in between. 'exon' by default. Rows in the annotation with a matched feature will be extracted and used for read mapping.
    attribute = "gene_name"   # Specify attribute type in GTF annotation. 'gene_id' by  default. Meta-features used for read counting will be  extracted from annotation using the provided value.
    overlap = 1 # Minimum number of overlapping bases in a read that is  required for read assignment. 1 by default. Number of overlapping bases is counted from both reads if paired end.
"""
else:
    print(
        f"'{args.quantification}' is not a valid quantification option. Check generate_config .py --help for options."
    )
    exit()

# VARIANT CALLING
if args.variant == "none":
    q = "# No variant calling\n"
elif args.variant == "snpeff":
    q = """# Variant calling
[snpeff]
    
"""


# Create config file
general = """# TOML config file for {0}

# Required information
[general]
    project = "{0}"   # Choose a project name. Should not start with a number.
    fastq = "/lustre03/project/6019267/shared/data/testShortReadSeq" # "/path/to/raw/data"   # Path to raw fastq.gz files. Ex: /lustre03/project/6019267/shared/data/<project>.
    output = "/lustre04/scratch/mlab/pipeline2024/{0}/output" # "/path/to/output" # Preferably use your scratch. The directory will be created. Exemple: /lustre04/scratch/<user>/<project>/output/.
    temporary = "/lustre04/scratch/mlab/pipeline2024/{0}/tmp" # Preferably use your scratch. The directory will be created. Exemple: /lustre04/scratch/<user>/<project>/tmp/. 
    sequencing = "RNA" # Type of sequencing. Short read RNA or DNA. Possible values: ["RNA", "Exome", "Genome"].
    reads = "PE" # Type of reads sequencing. Either single-end or paired-end. Possible values: ["SE", "PE"].
    reference = "hg38" # Possible values: Human: ["hg19", "hg38"]. Mouse: ["mm39"]. Worm: ["ce11"]. Zebrafish: ["danRer11"].
    trimming = "{1}"
    alignment = "{2}"
    pseudo = "{3}"
    quantification = "{4}"
    variant = "{5}"
    threads = 4  # Number of threads.
    memory = 64 # Total memory needed.
    email = "marjorie.labrecque@umontreal.ca" # You e-mail adress to receive notification
    

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

""".format(
    project_name,
    args.trimming,
    args.alignment,
    args.pseudo,
    args.quantification,
    args.variant,
)  # noqa: F524

toml_config = general + "\n" + t + "\n" + a + "\n" + p + "\n" + q

print(toml_config, file=open(work_dir + "/" + project_name + ".config.toml", "w"))

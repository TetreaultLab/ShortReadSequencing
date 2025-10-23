import argparse
from pathlib import Path
import re
import sys
import zipfile
import toml
import os
import subprocess


def main():
    parser = argparse.ArgumentParser(
        prog="GenerateJob",
        description="Generate SLURM job for samples based on a TOML project config file.",
    )

    parser.add_argument("--sample", type=str, required=True, help="Sample name.")
    parser.add_argument(
        "--config", type=str, required=True, help="Project config file, including path."
    )
    parser.add_argument(
        "--test",
        action="store_true",
        required=False,
        help="Will not launch, only create the script.",
    )
    parser.add_argument(
        "--redo",
        action="store_true",
        required=False,
        help="Removes steps_done.txt file and start from scratch",
    )

    args = parser.parse_args()

    # SAMPLE NAME
    sample = args.sample

    # LOAD CONFIG FILE
    toml_file = args.config

    with open(toml_file, "r") as f:
        toml_config = toml.load(f)
    f.close()

    # Start of pipeline
    print(
        ">>> {}-seq pipeline starting for {}".format(
            toml_config["general"]["sequencing"], sample
        ),
        file=sys.stdout,
    )

    # Creating output and tmp directories for sample
    output = toml_config["general"]["output"] + "/" + sample
    tmp = toml_config["general"]["temporary"] + "/" + sample
    subprocess.run(["mkdir", "-p", output])
    subprocess.run(["mkdir", "-p", tmp])
    print(f"\n>>> Output saved to {output}\n")

    # Open file for steps done
    if args.redo:
        subprocess.run(["rm", output + "/steps_done.txt"])
    steps = open(output + "/steps_done.txt", "a")
    steps.close()

    done = []
    with open(output + "/steps_done.txt", "r") as f:
        for line in f:
            done.append(line.strip())

    # Get tools and versions
    print(">>> Parameters:")

    genome = get_reference(toml_config["general"]["reference"], "")["fasta"]
    print(f"\t>>> Reference genome version: {genome}")

    # Calling main and downstream
    try:
        steps_main(sample, toml_config, done)
        steps_downstream(sample, toml_config, done)
    except Exception as e:
        print(f"Error: {e}")
        exit(1)

    if args.test:
        print("\n\nTESTING MODE!\nThe pipeline will not be launched\n\n")
    else:
        # Call main.sh (Launch the pipeline)
        subprocess.run(["bash", output + "/scripts/main.sh"])


# Functions needed by others
def title(message):
    print(f"\n>>> Running {message} [{get_time()}]")


def get_file_trimmed(toml_config, output, sample):
    files = {}
    if toml_config["general"]["trimming"] == "bbduk":
        files = {
            "I1_toAlign": toml_config["general"]["output"]
            + "/"
            + sample
            + "/Trimmed/"
            + sample
            + "_trimmed_R1.fastq.gz",
            "I2_toAlign": toml_config["general"]["output"]
            + "/"
            + sample
            + "/Trimmed/"
            + sample
            + "_trimmed_R2.fastq.gz",
            "I_toAlign": toml_config["general"]["output"]
            + "/"
            + sample
            + "/Trimmed/"
            + sample
            + "_trimmed.fastq.gz",
            "O_aligned": output + "/" + sample + "_",
        }

    else:
        files = {
            "I1_toAlign": toml_config["general"]["fastq"]
            + "/"
            + sample
            + "_R1.fastq.gz",
            "I2_toAlign": toml_config["general"]["fastq"]
            + "/"
            + sample
            + "_R2.fastq.gz",
            "I_toAlign": toml_config["general"]["fastq"] + "/" + sample + ".fastq.gz",
            "O_aligned": output + "/" + sample + "_",
        }

    return files


def get_reference(ref, tool):
    path = "/lustre09/project/6019267/shared/tools/references/ensembl_release114/"
    reference: {}  # type: ignore
    match ref:
        case "grch37":
            reference = {
                "fasta": path + "GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa",
                "index": path + "index_" + tool + "/GRCh37",
                "gtf": path + "GRCh37/Homo_sapiens.GRCh37.87.gtf",
                "gff3": path + "GRCh37/Homo_sapiens.GRCh37.87.gff3.gz",
            }
        case "grch38":
            reference = {
                "fasta": path + "GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
                "index": path + "index_" + tool + "/GRCh38",
                "gtf": path + "GRCh38/Homo_sapiens.GRCh38.114.gtf",
                "gff3": path + "GRCh38/Homo_sapiens.GRCh38.114.gff3",
            }
        case "grcm39":
            reference = {
                "fasta": path + "GRCm39/Mus_musculus.GRCm39.dna.primary_assembly.fa",
                "index": path + "index_" + tool + "/GRCm39",
                "gtf": path + "GRCm39/Mus_musculus.GRCm39.114.gtf",
                "gff3": path + "GRCm39/Mus_musculus.GRCm39.114.gff3",
            }
        case "grcz11":
            reference = {
                "fasta": path + "GRCz11/Danio_rerio.GRCz11.dna.primary_assembly.fa",
                "index": path + "index_" + tool + "/GRCz11",
                "gtf": path + "GRCz11/Danio_rerio.GRCz11.114.gtf",
                "gff3": path + "GRCz11/Danio_rerio.GRCz11.114.gff3",
            }
        case "wbcel235":
            reference = {
                "fasta": path
                + "WBcel235/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa",
                "index": path + "index_" + tool + "/WBcel235",
                "gtf": path + "WBcel235/Caenorhabditis_elegans.WBcel235.114.gtf",
                "gff3": path + "WBcel235/Caenorhabditis_elegans.WBcel235.114.gff3",
            }

    return reference


def create_script(cores, memory, time, sample_name, step, email, command, output):
    steps_done = output + "/steps_done.txt"
    work_dir = os.getcwd()
    job = work_dir + "/" + sample_name + "_" + step + ".slurm"

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
    module load python/3.11.5 apptainer
    module load fastqc/0.12.1 bbmap/39.06 star/2.7.11a bwa-mem2/2.2.1 gcc/12.3 openmpi/4.1.5 salmon/1.10.2 samtools/1.22.1 gatk/4.6.1.0 subread/2.0.6 htslib/1.22.1 bcftools/1.22 freebayes/1.3.7
    source /lustre09/project/6019267/shared/tools/main_pipelines/long-read/launch_pipeline_env/bin/activate
    export JAVA_TOOL_OPTIONS="-Xmx{1}g -XX:ParallelGCThreads={0}"
    #
    ### Launch script
    #
    newgrp rrg-tetreaum
    #
    {6}
    #
    """.format(cores, memory, time, sample_name, step, email, command)

    slurm += f'\nif [ $? -eq 0 ]; then echo "{step}" >> "{steps_done}"; fi\n\n'

    with open(job, "w") as o:
        o.write(slurm)

    return job


## Tool functions
def fastqc(sample, toml_config):
    title("FastQC")

    output = toml_config["general"]["output"] + "/" + sample + "/QC/fastQC"
    subprocess.run(["mkdir", "-p", output])

    temporary = toml_config["general"]["temporary"] + "/" + sample
    subprocess.run(["mkdir", "-p", temporary])

    if toml_config["general"]["reads"] == "PE":
        Read1 = toml_config["general"]["fastq"] + "/" + sample + "_R1.fastq.gz"
        Read2 = toml_config["general"]["fastq"] + "/" + sample + "_R2.fastq.gz"

        command = [
            "fastqc",
            "-o",
            output,
            "--noextract",
            "--threads",
            "6",
            "--dir",
            temporary,
            "--kmers",
            str(toml_config["fastqc"]["kmers"]),
            Read1,
            Read2,
        ]
        command_str = " ".join(command)

    else:
        Read = toml_config["general"]["fastq"] + "/" + sample + ".fastq.gz"

        command = [
            "fastqc",
            "-o",
            output,
            "--noextract",
            "--threads",
            "6",
            "--dir",
            temporary,
            "--kmers",
            str(toml_config["fastqc"]["kmers"]),
            Read,
        ]

        command_str = " ".join(command)

    steps_done = toml_config["general"]["output"] + "/" + sample + "/steps_done.txt"
    command_str += f'\nif [ $? -eq 0 ]; then echo "FastQC" >> "{steps_done}"; fi\n\n'

    return command_str


def bbduk(sample, toml_config):
    title("BBDuk")

    output = toml_config["general"]["output"] + "/" + sample + "/Trimmed"
    subprocess.run(["mkdir", "-p", output])

    if toml_config["general"]["reads"] == "PE":
        I1duk = toml_config["general"]["fastq"] + "/" + sample + "_R1.fastq.gz"
        I2duk = toml_config["general"]["fastq"] + "/" + sample + "_R2.fastq.gz"
        O1duk = output + "/" + sample + "_trimmed_R1.fastq.gz"
        O2duk = output + "/" + sample + "_trimmed_R2.fastq.gz"

        command = [
            "bbduk.sh",
            "in=" + I1duk,
            "in2=" + I2duk,
            "out=" + O1duk,
            "out2=" + O2duk,
            "stats=" + output + "/contaminants_stats.txt",
            "threads=" + str(toml_config["general"]["threads"]),
            "ordered=" + toml_config["bbduk"]["ordered"],
            "ref=/lustre09/project/6019267/shared/tools/main_pipelines/short-read/adapters.fa",
            "k=" + str(toml_config["bbduk"]["kmers"]),
            "mink=" + str(toml_config["bbduk"]["mink"]),
            "ktrim=" + str(toml_config["bbduk"]["ktrim"]),
            "trimq=" + str(toml_config["bbduk"]["trimq"]),
            "qtrim=" + toml_config["bbduk"]["qtrim"],
            "minlength=" + str(toml_config["bbduk"]["minlength"]),
            "mlf=" + str(toml_config["bbduk"]["mlf"]),
            "minavgquality=" + str(toml_config["bbduk"]["minavgquality"]),
            "tbo=t",
            "tpe=t",
        ]
    else:
        Iduk = toml_config["general"]["fastq"] + "/" + sample + ".fastq.gz"
        Oduk = output + "/" + sample + "_trimmed.fastq.gz"
        command = [
            "bbduk.sh",
            "in=" + Iduk,
            "out=" + Oduk,
            "stats=" + output + "/contaminants_stats.txt",
            "threads=" + str(toml_config["general"]["threads"]),
            "ordered=" + toml_config["bbduk"]["ordered"],
            "ref=/lustre09/project/6019267/shared/tools/main_pipelines/short-read/adapters.fa",
            "k=" + str(toml_config["bbduk"]["kmers"]),
            "mink=" + str(toml_config["bbduk"]["mink"]),
            "ktrim=" + str(toml_config["bbduk"]["ktrim"]),
            "trimq=" + str(toml_config["bbduk"]["trimq"]),
            "qtrim=" + toml_config["bbduk"]["qtrim"],
            "minlength=" + str(toml_config["bbduk"]["minlength"]),
            "mlf=" + str(toml_config["bbduk"]["mlf"]),
            "minavgquality=" + str(toml_config["bbduk"]["minavgquality"]),
            "tbo=t",
            "tpe=t",
        ]

    command_str = " ".join(command)

    steps_done = toml_config["general"]["output"] + "/" + sample + "/steps_done.txt"
    command_str += f'\nif [ $? -eq 0 ]; then echo "BBDuk" >> "{steps_done}"; fi\n\n'

    return command_str


def star(sample, toml_config):
    title("STAR")
    output = toml_config["general"]["output"] + "/" + sample + "/Aligned"
    subprocess.run(["mkdir", "-p", output])

    temporary = toml_config["general"]["temporary"] + "/" + sample + "/star_tmp"
    subprocess.run(["rm", "-r", temporary])

    ref = get_reference(toml_config["general"]["reference"], "star")["index"]

    files = get_file_trimmed(toml_config, output, sample)
    I1_toAlign = files["I1_toAlign"]
    I2_toAlign = files["I2_toAlign"]
    I_toAlign = files["I_toAlign"]
    O_aligned = files["O_aligned"]

    if toml_config["general"]["reads"] == "PE":
        command = [
            "STAR",
            "--runMode",
            "alignReads",
            "--runThreadN",
            str(toml_config["general"]["threads"]),
            "--limitBAMsortRAM",
            "60000000000",
            "--genomeDir",
            ref,
            "--outFileNamePrefix",
            O_aligned,
            "--outSAMattributes",
            "Standard",
            "--outTmpDir",
            temporary,
            "--outReadsUnmapped",
            "Fastx",
            "--sjdbOverhang",
            "99",
            "--readFilesCommand",
            "zcat",
            "--readFilesIn",
            I1_toAlign,
            I2_toAlign,
            "--outSAMtype",
            toml_config["star"]["outSAMtype1"],
            toml_config["star"]["outSAMtype2"],
            "--twopassMode",
            toml_config["star"]["twopassMode"],
            "--outSJtype",
            toml_config["star"]["outSJtype"],
            "--quantMode",
            toml_config["star"]["quantMode"],
        ]
    else:
        command = [
            "STAR",
            "--runMode",
            "alignReads",
            "--runThreadN",
            str(toml_config["general"]["threads"]),
            "--limitBAMsortRAM",
            "60000000000",
            "--genomeDir",
            ref,
            "--outSAMattributes",
            "Standard",
            "--outFileNamePrefix",
            O_aligned,
            "--outTmpDir",
            temporary,
            "--readFilesCommand",
            "zcat",
            "--readFilesIn",
            I_toAlign,
            "--outSAMtype",
            toml_config["star"]["outSAMtype1"],
            toml_config["star"]["outSAMtype2"],
            "--twopassMode",
            toml_config["star"]["twopassMode"],
            "--outSJtype",
            toml_config["star"]["outSJtype"],
            "--quantMode",
            toml_config["star"]["quantMode"],
        ]
    command_str = " ".join(command)

    command_str += (
        f"\nmv {output}/{sample}_Aligned.sortedByCoord.out.bam {output}/{sample}.bam"
    )
    command_str += f"\nmv {output}/{sample}_Log.final.out {output}/{sample}_summary_mapping_stats.out"
    command_str += (
        f"\nmv {output}/{sample}_Log.out {output}/{sample}_run_information.out"
    )
    command_str += f"\nrm {output}/{sample}_Log.progress.out"
    command_str += f"\nrm {output}/{sample}_Unmapped.out.mate1"
    command_str += f"\nrm {output}/{sample}_Unmapped.out.mate2"
    command_str += f"\nrm -r {output}/{sample}__STARpass1"
    command_str += f"\nmv {output}/{sample}__STARgenome/sjdbList.out.tab {output}/{sample}_sjdbList.out.tab"

    steps_done = toml_config["general"]["output"] + "/" + sample + "/steps_done.txt"
    command_str += f'\nif [ $? -eq 0 ]; then echo "STAR" >> "{steps_done}"; fi\n\n'
    return command_str


def bwa(sample, toml_config):
    title("BWA-MEM2")

    output = toml_config["general"]["output"] + "/" + sample + "/Aligned"
    subprocess.run(["mkdir", "-p", output])

    ref = get_reference(toml_config["general"]["reference"], "bwa-mem2")["index"]

    files = get_file_trimmed(toml_config, output, sample)
    I1_toAlign = files["I1_toAlign"]
    I2_toAlign = files["I2_toAlign"]
    I_toAlign = files["I_toAlign"]

    if toml_config["general"]["reads"] == "PE":
        command = [
            "bwa-mem2",
            "mem",
            "-o",
            output + "/" + sample + ".sam",
            "-t",
            str(toml_config["general"]["threads"]),
            "-v",
            "1",
            ref,
            I1_toAlign,
            I2_toAlign,
        ]
    else:
        command = [
            "bwa-mem2",
            "mem",
            "-o",
            output + "/" + sample + ".sam",
            "-t",
            str(toml_config["general"]["threads"]),
            "-v",
            "1",
            ref,
            I_toAlign,
        ]

    command_str = " ".join(command)

    # Change SAM to BAM format and remove SAM
    command_str += (
        f"\nsamtools view -S -b {output}/{sample}.sam -o {output}/{sample}.bam",
    )
    command_str += f"\nrm {output}/{sample}.sam"

    steps_done = toml_config["general"]["output"] + "/" + sample + "/steps_done.txt"
    command_str += f'\nif [ $? -eq 0 ]; then echo "BWA-MEM2" >> "{steps_done}"; fi\n\n'
    return command_str


def salmon(sample, toml_config):
    title("Salmon")

    output = toml_config["general"]["output"] + "/" + sample + "/Salmon"
    subprocess.run(["mkdir", "-p", output])

    ref = get_reference(toml_config["general"]["reference"], "salmon")["index"]

    files = get_file_trimmed(toml_config, output, sample)
    I1_toAlign = files["I1_toAlign"]
    I2_toAlign = files["I2_toAlign"]
    I_toAlign = files["I_toAlign"]

    if toml_config["general"]["reads"] == "PE":
        command = [
            "salmon",
            "quant",
            "--libType",
            "A",
            "--validateMappings",
            "--index",
            ref,
            "--threads",
            str(toml_config["general"]["threads"]),
            "--auxDir",
            "salmon_tmp",
            "--minScoreFraction",
            str(toml_config["salmon"]["minScoreFraction"]),
            "--mates1",
            I1_toAlign,
            "--mates2",
            I2_toAlign,
            "--output",
            output,
        ]
    else:
        command = [
            "salmon",
            "quant",
            "--libType",
            "A",
            "--index",
            ref,
            "--threads",
            str(toml_config["general"]["threads"]),
            "--auxDir",
            "salmon_tmp",
            "--minScoreFraction",
            str(toml_config["salmon"]["minScoreFraction"]),
            "--unmatedReads",
            I_toAlign,
            "--output",
            output,
        ]

    command_str = " ".join(command)

    command_str += f"\nmv {output}/logs/salmon_quant.log {output}/{sample}_log.out"
    command_str += f"\nmv  {output}/quant.sf {output}/{sample}_transcript_quant.sf"
    command_str += f"\nrm {output}/cmd_info.json"
    command_str += f"\nrm -r {output}/libParams"
    command_str += f"\nrm -r {output}/logs"
    command_str += f"\nrm -r {output}/salmon_tmp"

    steps_done = toml_config["general"]["output"] + "/" + sample + "/steps_done.txt"
    command_str += f'\nif [ $? -eq 0 ]; then echo "Salmon" >> "{steps_done}"; fi\n\n'
    return command_str


def samtools(sample, toml_config):
    title("Samtools")
    in_out = toml_config["general"]["output"] + "/" + sample + "/Aligned"

    inBAM = in_out + "/" + sample + ".bam"
    bamCoord = in_out + "/" + sample + "_sortedCoordinate.bam"
    stats = in_out + "/" + sample + "_stats.txt"

    # Sort by coordinate
    command_str = f"\nsamtools sort {inBAM} -o {bamCoord}"

    # Index bam sorted by coordinates
    command_str = f"\nsamtools", "index", "-b", {bamCoord}, "-o", {bamCoord} + ".bai"

    # alignment stats
    command_str = f"\nsamtools stats {bamCoord} | grep ^SN | cut -f, 2- > {stats}"
    command_str = f"\nrm {inBAM}"

    steps_done = toml_config["general"]["output"] + "/" + sample + "/steps_done.txt"
    command_str += f'\nif [ $? -eq 0 ]; then echo "Samtools" >> "{steps_done}"; fi\n\n'
    return command_str


def bamqc(sample, toml_config):
    title("FastQC for bam")

    output = toml_config["general"]["output"] + "/" + sample + "/QC/fastQC"
    subprocess.run(["mkdir", "-p", output])

    temporary = toml_config["general"]["temporary"] + "/" + sample
    subprocess.run(["mkdir", "-p", temporary])

    input = (
        toml_config["general"]["output"]
        + "/"
        + sample
        + "/Aligned/"
        + sample
        + "_sortedCoordinate.bam"
    )

    command = [
        "fastqc",
        "-o",
        output,
        "--noextract",
        "--threads",
        "2",
        "--dir",
        temporary,
        "--kmers",
        str(toml_config["fastqc"]["kmers"]),
        input,
    ]

    command_str = " ".join(command)
    steps_done = output + "/steps_done.txt"
    command_str += (
        f'\nif [ $? -eq 0 ]; then echo "FastQC for bam" >> "{steps_done}"; fi\n\n'
    )

    command_str += f"python -u /lustre09/project/6019267/shared/tools/main_pipelines/short-read/ShortReadSequencing/check_fasqc.py --path {output} --sample {sample}"

    return command_str


def markduplicates(sample, toml_config):
    title("MarkDuplicates")

    output = toml_config["general"]["output"] + "/" + sample + "/MarkDuplicates/"
    subprocess.run(["mkdir", "-p", output])

    temporary = toml_config["general"]["temporary"] + "/" + sample + "/md_tmp"
    subprocess.run(["mkdir", "-p", temporary])

    input = toml_config["general"]["output"] + "/" + sample + "/Aligned"
    bamCoord = input + "/" + sample + "_sortedCoordinate.bam"
    metrics = output + sample + "_duplicates_metrics.txt"
    records = output + sample + "_markDuplicates.bam"

    command = [
        "gatk",
        "MarkDuplicates",
        "--INPUT",
        bamCoord,
        "--METRICS_FILE",
        metrics,
        "--REMOVE_DUPLICATES",
        toml_config["markduplicates"]["remove"],
        "--TMP_DIR",
        temporary,
        "--OUTPUT",
        records,
        "--CREATE_INDEX",
        toml_config["markduplicates"]["index"],
        "--DUPLICATE_SCORING_STRATEGY",
        toml_config["markduplicates"]["strategy"],
    ]

    command_str = " ".join(command)
    print(f">>> {command_str}\n")
    subprocess.run(command, check=True)

    subprocess.run("rm", "-r", temporary)

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("MarkDuplicates\n")


def featurecounts(sample, toml_config):
    title("FeatureCounts")

    temporary = toml_config["general"]["temporary"] + "/" + sample
    output = toml_config["general"]["output"] + "/" + sample + "/FeatureCounts"
    subprocess.run(["mkdir", "-p", output])

    input = (
        toml_config["general"]["output"]
        + "/"
        + sample
        + "/Aligned/"
        + sample
        + "_sortedCoordinate.bam"
    )

    gtf = get_reference(toml_config["general"]["reference"], "")["gtf"]

    if toml_config["general"]["reads"] == "PE":
        command = [
            "featureCounts",
            "-t",
            toml_config["featurecounts"]["features"],
            "-O",
            "--countReadPairs",
            "-F",
            "GTF",
            "-g",
            toml_config["featurecounts"]["attribute"],
            "--minOverlap",
            str(toml_config["featurecounts"]["overlap"]),
            "-p",
            "-T",
            str(toml_config["general"]["threads"]),
            "--tmpDir",
            temporary,
            "-a",
            gtf,
            "-o",
            output + "/" + sample + "_geneID.txt",
            input,
        ]
    else:
        command = [
            "featureCounts",
            "-t",
            toml_config["featurecounts"]["features"],
            "-O",
            "-F",
            "GTF",
            "-g",
            toml_config["featurecounts"]["attribute"],
            "--minOverlap",
            str(toml_config["featurecounts"]["overlap"]),
            "-T",
            str(toml_config["general"]["threads"]),
            "--tmpDir",
            temporary,
            "-a",
            gtf,
            "-o",
            output + "/" + sample + "_geneID.txt",
            input,
        ]

    command_str = " ".join(command)
    print(f">>> {command_str}\n")
    subprocess.run(command, check=True)

    # Keep Geneid and counts
    with open(output + "/" + sample + "_0.counts", "w") as outfile:
        subprocess.run(
            ["tail", "-n", "+2", output + "/" + sample + "_geneID.txt"], stdout=outfile
        )

    with open(output + "/" + sample + "_counts.txt", "w") as outfile:
        subprocess.run(
            ["cut", "-f1,7", output + "/" + sample + "_0.counts"], stdout=outfile
        )

    subprocess.run(
        [
            "sed",
            "-i",
            "1s|.*|gene_id\t" + sample + "|",
            output + "/" + sample + "_counts.txt",
        ]
    )
    subprocess.run(["rm", output + "/" + sample + "_0.counts"])

    # Steps done
    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("FeatureCounts\n")


def multiqc(sample, toml_config):
    title("MultiQC")
    input = toml_config["general"]["output"] + "/" + sample + "/"
    output = toml_config["general"]["output"] + "/" + sample + "/QC/multiQC/"
    subprocess.run(["mkdir", "-p", output])

    with open(output + "/my_file_list.txt", "w") as f:
        f.write(input + "QC/fastQC/\n")
        if toml_config["general"]["trimming"] != "False":
            f.write(input + "Trimmed/\n")
        if toml_config["general"]["pseudo"] != "False":
            f.write(input + "Salmon/\n")
        if toml_config["general"]["quantification"] != "False":
            f.write(input + "FeatureCounts/\n")
        f.write(input + "Aligned/\n")
        f.write(input + "Samtools/\n")
        f.write(input + "MarkDuplicates/\n")
    f.close()

    subprocess.run(
        [
            "apptainer",
            "run",
            "/lustre09/project/6019267/shared/tools/others/multiqc/multiqc.sif",
            "multiqc",
            "--file-list",
            output + "my_file_list.txt",
            "--force",
            "--filename",
            sample + "_multiqc_report",
            "--outdir",
            output,
        ]
    )

    subprocess.run(["rm", output + "my_file_list.txt"])
    subprocess.run(["rm", "-r", output + "/" + sample + "_multiqc_report_data"])

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("MultiQC\n")


def bcftools(sample, toml_config):
    title("BCFtools")

    input = (
        toml_config["general"]["output"]
        + "/"
        + sample
        + "/Aligned/"
        + sample
        + "_sortedCoordinate.bam"
    )
    output = toml_config["general"]["output"] + "/" + sample + "/Variants/"
    subprocess.run(["mkdir", "-p", output])

    ref = get_reference(toml_config["general"]["reference"], "")["fasta"]

    with open(
        toml_config["general"]["output"] + "/" + sample + "/sample.txt", "w"
    ) as sample_file:
        subprocess.run(["echo", sample], stdout=sample_file)

    mpileup = [
        "bcftools",
        "mpileup",
        "--threads",
        str(toml_config["general"]["threads"]),
        "-d",
        "1000000",
        "--max-idepth",
        "1000000",
        "--count-orphans",
        "-o",
        output + sample + "_mpileup.vcf",
        "-f",
        ref,
        input,
    ]

    call = [
        "bcftools",
        "call",
        "--threads",
        str(toml_config["general"]["threads"]),
        "--multiallelic-caller",
        "--variants-only",
        "-o",
        output + sample + "_bcftools.vcf",
        output + sample + "_mpileup.vcf",
    ]

    command_1 = " ".join(mpileup)
    print(f">>> {command_1}\n")
    subprocess.run(mpileup, check=True)

    command_2 = " ".join(call)
    print(f">>> {command_2}\n")
    subprocess.run(call, check=True)

    subprocess.run(["rm", output + sample + "_mpileup.vcf"])
    subprocess.run(
        [
            "bcftools",
            "reheader",
            "--samples",
            toml_config["general"]["output"] + "/" + sample + "/sample.txt",
            "-o",
            output + sample + "_bcftools_header.vcf",
            output + sample + "_bcftools.vcf",
        ],
        check=True,
    )

    command_norm = [
        "bcftools",
        "norm",
        "-m-any",
        "-o",
        output + sample + "_bcftools_norm.vcf",
        output + sample + "_bcftools_header.vcf",
    ]

    command_norm_str = " ".join(command_norm)
    print(f">>> {command_norm_str}\n")
    subprocess.run(command_norm, check=True)

    subprocess.run(["bgzip", "-f", output + sample + "_bcftools_norm.vcf"], check=True)
    subprocess.run(
        ["tabix", "-p", "vcf", output + sample + "_bcftools_norm.vcf.gz"], check=True
    )

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("BCFtools\n")


def freebayes(sample, toml_config):
    title("FreeBayes")

    input = (
        toml_config["general"]["output"]
        + "/"
        + sample
        + "/Aligned/"
        + sample
        + "_sortedCoordinate.bam"
    )
    output = toml_config["general"]["output"] + "/" + sample + "/Variants/"
    subprocess.run(["mkdir", "-p", output])

    ref = get_reference(toml_config["general"]["reference"], "")["fasta"]

    # freebayes with legacy genotype calling and naive variant calling
    command = [
        "freebayes",
        "-f",
        ref,
        # "--legacy-gls",
        # "--min-alternate-fraction",
        # "0.01",
        # "--ploidy",
        # "50",
        # "--min-repeat-entropy",
        # "1",
        # "--no-partial-observations",
        input,
        "--vcf",
        output + sample + "_freebayes.vcf",
    ]

    command_str = " ".join(command)
    print(f">>> {command_str}\n")
    subprocess.run(command, check=True)

    subprocess.run(
        [
            "bcftools",
            "reheader",
            "--samples",
            toml_config["general"]["output"] + "/" + sample + "/sample.txt",
            "-o",
            output + sample + "_freebayes_header.vcf",
            output + sample + "_freebayes.vcf",
        ],
        check=True,
    )

    command_norm = [
        "bcftools",
        "norm",
        "-m-any",
        "-o",
        output + sample + "_freebayes_norm.vcf",
        output + sample + "_freebayes_header.vcf",
    ]

    command_norm_str = " ".join(command_norm)
    print(f">>> {command_norm_str}\n")
    subprocess.run(command_norm, check=True)

    subprocess.run(["bgzip", "-f", output + sample + "_freebayes_norm.vcf"], check=True)
    subprocess.run(
        ["tabix", "-p", "vcf", output + sample + "_freebayes_norm.vcf.gz"], check=True
    )
    subprocess.run(
        ["rm", toml_config["general"]["output"] + "/" + sample + "/sample.txt"]
    )

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("FreeBayes\n")


def bcftools_filter(sample, toml_config):
    title("BCFtools filters")

    output = toml_config["general"]["output"] + "/" + sample + "/Variants/"
    subprocess.run(["mkdir", "-p", output])

    command_concat = [
        "bcftools",
        "concat",
        "--threads",
        str(toml_config["general"]["threads"]),
        "--allow-overlaps",
        "--rm-dups",
        "exact",
        output + sample + "_bcftools_norm.vcf.gz",
        output + sample + "_freebayes_norm.vcf.gz",
        "-o",
        output + sample + "_merged.vcf.gz",
    ]

    command_concat_str = " ".join(command_concat)
    print(f">>> {command_concat_str}\n")
    subprocess.run(command_concat, check=True)

    command_filter = [
        "bcftools",
        "filter",
        "-i",
        "QUAL>1 && INFO/DP>0 && AC>0",
        "-o",
        output + sample + "_unannotated.vcf",
        output + sample + "_merged.vcf.gz",
    ]

    command_filter_str = " ".join(command_filter)
    print(f">>> {command_filter_str}\n")
    subprocess.run(command_filter, check=True)

    # remove intermediate vcf files
    # subprocess.run(["rm", output + sample + "_unannotated.vcf"])
    # subprocess.run(["rm", output + sample + "_merged.vcf.gz"])
    # subprocess.run(["rm", output + sample + "_freebayes.vcf"])
    # subprocess.run(["rm", output + sample + "_bcftools.vcf"])

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("BCFtools filters\n")


def snpeff(sample, toml_config):
    title("SnpEff")

    snpeff = "/lustre09/project/6019267/shared/tools/main_pipelines/short-read/snpEff"

    genome = toml_config["general"]["reference"]

    if genome == "grch37":
        ref = "GRCh37.87"
    if genome == "grch38":
        ref = "GRCh38.105"
    if genome == "grcz11":
        ref = "GRCz11.105"
    if genome == "grcm39":
        ref = "GRCm39.105"
    if genome == "wbcel235":
        ref = "WBcel235.105"

    path = toml_config["general"]["output"] + "/" + sample + "/Variants"

    cmd_snpeff = [
        "java",
        "-jar",
        snpeff + "/snpEff.jar",
        "-noLog",
        ref,
        "-c",
        snpeff + "/snpEff.config",
        "-stats",
        path + "/" + sample + "_summary.html",
        "-csvStats",
        path + "/" + sample + "_summary.csv",
        path + "/" + sample + "_unannotated.vcf",
    ]

    command_str1 = " ".join(cmd_snpeff)
    print(f">>> {command_str1}\n")
    with open(path + "/" + sample + "_snpeff.vcf", "w") as outfile:
        subprocess.run(
            cmd_snpeff,
            stdout=outfile,
            check=True,
        )

    cmd_vartype = [
        "java",
        "-jar",
        snpeff + "/SnpSift.jar",
        "varType",
        "-noLog",
        path + "/" + sample + "_snpeff.vcf",
    ]

    command_str2 = " ".join(cmd_vartype)
    print(f">>> {command_str2}\n")
    with open(path + "/" + sample + "_annotated.vcf", "w") as outfile:
        subprocess.run(
            cmd_vartype,
            stdout=outfile,
            check=True,
        )

    cmd_extract = [
        "java",
        "-jar",
        snpeff + "/SnpSift.jar",
        "extractFields",
        "-noLog",
        "-s",
        "|",
        "-e",
        ".",
        path + "/" + sample + "_annotated.vcf",
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "QUAL",
        "HOM",
        "DP",
        # "TYPE",
        "ANN[*].GENE",
        "ANN[*].GENEID",
        "ANN[*].FEATUREID",
        "ANN[*].EFFECT",
        "ANN[*].IMPACT",
        "ANN[*].BIOTYPE",
        "ANN[*].HGVS_C",
        "ANN[*].HGVS_P",
    ]

    command_str3 = " ".join(cmd_extract)
    print(f">>> {command_str3}\n")
    with open(path + "/" + sample + "_annotated.txt", "w") as outfile:
        subprocess.run(
            cmd_extract,
            stdout=outfile,
            check=True,
        )

    # subprocess.run(["rm", path + "/" + sample + "_summary.csv"])
    # subprocess.run(["rm", path + "/" + sample + "_snpeff.vcf"])

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("SnpEff\n")


def dbNSFP(sample, toml_config):
    title("Add dbNSFP to snpEff output")

    genome = toml_config["general"]["reference"]
    path = toml_config["general"]["output"] + "/" + sample + "/Variants"

    chromosomes = [
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "10",
        "11",
        "12",
        "13",
        "14",
        "15",
        "16",
        "17",
        "18",
        "19",
        "20",
        "21",
        "22",
        "X",
        "Y",
        "MT",
    ]

    var = pd.read_csv(
        path + "/" + sample + "_annotated.txt", header=0, sep="\t", low_memory=False
    )
    var = var.rename(columns={"CHROM": "CHROM_" + genome, "POS": "POS_" + genome})

    var = var.astype({"CHROM_" + genome: str, "POS_" + genome: str})

    appended_data = []
    ref = "/lustre09/project/6019267/shared/tools/main_pipelines/short-read/dbNSFP"

    for chromosome in chromosomes:
        print(chromosome)
        db = pd.read_csv(
            ref + "/dbNSFP4.7a_variant.chr" + chromosome + "_small.txt",
            header=0,
            sep="\t",
            low_memory=False,
        )

        db = db.astype({"CHROM_" + genome: str, "POS_" + genome: str})

        var_chr = var[var["CHROM_" + genome] == chromosome]

        m = pd.merge(
            var_chr,
            db,
            how="left",
            on=["CHROM_" + genome, "POS_" + genome, "REF", "ALT"],
        )
        appended_data.append(m)

    final = pd.concat(appended_data)
    final.to_csv(path + "/" + sample + "_annotated_dbNSFP.txt", sep="\t", index=False)

    n_rows = len(final.index)
    rows_ann = final["rsID"].count()
    percentage = round(rows_ann / n_rows * 100, 2)

    print(
        "\nFinal\nannotated variants: ",
        rows_ann,
        "\ntotal variants:",
        n_rows,
        "\npercentage dbNSFP annotation: ",
        percentage,
        "%",
    )

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("dbNSFP\n")


def formatting(sample, toml_config):
    title("Format variants output")

    genome = toml_config["general"]["reference"]
    path = toml_config["general"]["output"] + "/" + sample + "/Variants"

    if genome == "grch37" or genome == "grch38":
        final = pd.read_csv(
            path + "/" + sample + "_annotated_dbNSFP.txt",
            sep="\t",
            header=0,
            low_memory=False,
        )

        final = final.rename(
            columns={
                "REF": "Ref",
                "ALT": "Alt",
                "QUAL": "Quality",
                # "TYPE": "Variation",
                "HOM": "Zygosity",
                "DP": "Read_depth",
                "ANN[*].GENE": "Gene_name",
                "ANN[*].GENEID": "Gene_id",
                "ANN[*].FEATUREID": "Transcript",
                "ANN[*].EFFECT": "Effect",
                "ANN[*].IMPACT": "Impact",
                "ANN[*].BIOTYPE": "Biotype",
                "ANN[*].HGVS_C": "Codon_change",
                "ANN[*].HGVS_P": "Protein_change",
            }
        )

        final["Position_" + genome] = (
            final["CHROM_" + genome].astype(str)
            + ":"
            + final["POS_" + genome].astype(str)
        )
        final["Quality"] = round(final["Quality"], 2)
        final = final.replace({"Zygosity": {True: "Hom", False: "Het"}})

        columns = [
            "Gene_name",
            "Gene_id",
            "Transcript",
            "Effect",
            "Impact",
            "Biotype",
            "Codon_change",
            "Protein_change",
        ]

        for index, row in final.iterrows():
            infos = []

            for c in columns:
                cell_value = str(row[c])
                str_split = [sub.replace('"', "") for sub in cell_value.split("|")]

                # Append to infos, or take the first split value if there's only one
                infos.append(str_split if len(str_split) > 1 else [cell_value])

                # Update the final DataFrame with the first part of split value
                final.at[index, c] = str_split[0]

                # Transpose and concatenate corresponding elements from each list
                if len(infos) > 1:
                    info_concat = [
                        "|".join([info[i] for info in infos])
                        for i in range(len(infos[0]))
                    ]
                    final.at[index, "Infos"] = "; ".join(info_concat)
                else:
                    final.at[index, "Infos"] = "; ".join(infos[0])

        final = final[
            [
                "Position_" + genome,
                "Ref",
                "Alt",
                "Quality",
                "Read_depth",
                "Zygosity",
                "rsID",
                # "Variation",
                "Gene_name",
                "Gene_id",
                "Transcript",
                "Effect",
                "Impact",
                "Biotype",
                "Codon_change",
                "Protein_change",
                "Infos",
                "genename",
                "Ensembl_geneid",
                "Ensembl_transcriptid",
                "Ensembl_proteinid",
                "codon_change",
                "protein_change",
                "SIFT_score",
                "PolyPhen2_HDIV_score",
                "GERP_score",
                "phyloP_score",
                "phastCons_score",
                "MutationTaster_score",
                "FATHMM_score",
                "REVEL_score",
                "AlphaMissense_score",
                "CADD_phred_score",
                "1000G_AF",
                "ExAC_AF",
                "gnomAD_exomes_AF",
                "gnomAD_exomes_NFE_AF",
                "gnomAD_genomes_AF",
                "gnomAD_genomes_NFE_AF",
                "clinvar_id",
                "clinvar_trait",
                "clinvar_clnsig",
                "OMIM",
            ]
        ]

        # count "Het" for each gene
        final["het_count"] = final.groupby("Gene_name")["Zygosity"].transform(
            lambda x: (x == "Het").sum()
        )

        # Remplace "Het" for "Multiple-het" if het_count > 1
        final.loc[
            (final["het_count"] > 1) & (final["Zygosity"] == "Het"), "Zygosity"
        ] = "Multiple-het"

        # Remove het_count
        final = final.drop(columns=["het_count"])

        final = final.replace(".", np.nan)
        final.to_csv(path + "/" + sample + "_variants_all.txt", sep="\t", index=False)

        ## Filtering
        # https://www.htslib.org/workflow/filter.html
        # https://jp.support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/QUAL_QD_GQ_Formulation_fDG.htm

        df_filtered = final[(final["Quality"] > 20) & (final["Read_depth"] > 5)]

        df_filtered.to_csv(
            path + "/" + sample + "_variants_filtered.txt", sep="\t", index=False
        )

        print(">>> Filters:")
        print("\t>>> Quality Phred > 20")
        print("\t>>> Number of total reads > 5")

        all_var = len(final.index)
        filtered_var = len(df_filtered.index)
        percentage = round(filtered_var / all_var * 100, 2)

        print(
            "\nNumber of variants in all: ",
            all_var,
            "\nNumber of variants in filtered: ",
            filtered_var,
            "\npercentage filtered variants: ",
            percentage,
            "%",
        )

    else:
        final = pd.read_csv(
            path + "/" + sample + "_annotated.txt", header=0, sep="\t", low_memory=False
        )

        final = final.rename(
            columns={
                "REF": "Ref",
                "ALT": "Alt",
                "QUAL": "Quality",
                "TYPE": "Variation",
                "HOM": "Zygosity",
                "DP": "Read_depth",
                "ANN[*].GENE": "Gene_name",
                "ANN[*].GENEID": "Gene_id",
                "ANN[*].FEATUREID": "Transcript",
                "ANN[*].EFFECT": "Effect",
                "ANN[*].IMPACT": "Impact",
                "ANN[*].BIOTYPE": "Biotype",
                "ANN[*].HGVS_C": "Codon_change",
                "ANN[*].HGVS_P": "Protein_change",
            }
        )

        final["Position"] = final["CHROM"].astype(str) + ":" + final["POS"].astype(str)
        final["Quality"] = round(final["Quality"], 2)
        final = final.replace({"Zygosity": {True: "Hom", False: "Het"}})

        columns = [
            "Gene_name",
            "Gene_id",
            "Transcript",
            "Effect",
            "Impact",
            "Biotype",
            "Codon_change",
            "Protein_change",
        ]

        for index, row in final.iterrows():
            infos = []

            for c in columns:
                cell_value = str(row[c])
                str_split = [sub.replace('"', "") for sub in cell_value.split("|")]

                # Append to infos, or take the first split value if there's only one
                infos.append(str_split if len(str_split) > 1 else [cell_value])

                # Update the final DataFrame with the first part of split value
                final.at[index, c] = str_split[0]

                # Transpose and concatenate corresponding elements from each list
                if len(infos) > 1:
                    info_concat = [
                        "|".join([info[i] for info in infos])
                        for i in range(len(infos[0]))
                    ]
                    final.at[index, "Infos"] = "; ".join(info_concat)
                else:
                    final.at[index, "Infos"] = "; ".join(infos[0])

        final = final[
            [
                "Position",
                "Ref",
                "Alt",
                "Quality",
                "Read_depth",
                "Zygosity",
                "Variation",
                "Gene_name",
                "Gene_id",
                "Transcript",
                "Effect",
                "Impact",
                "Biotype",
                "Codon_change",
                "Protein_change",
                "Infos",
            ]
        ]

        # count "Het" for each gene
        final["het_count"] = final.groupby("Gene_name")["Zygosity"].transform(
            lambda x: (x == "Het").sum()
        )

        # Remplace "Het" for "Multiple-het" if het_count > 1
        final.loc[
            (final["het_count"] > 1) & (final["Zygosity"] == "Het"), "Zygosity"
        ] = "Multiple-het"

        # Remove het_count
        final = final.drop(columns=["het_count"])

        final = final.replace(".", np.nan)
        print(final)

        final.to_csv(path + "/" + sample + "_variants_all.txt", sep="\t", index=False)

        df_filtered = final[(final["Quality"] > 20) & (final["Read_depth"] > 5)]
        df_filtered.to_csv(
            path + "/" + sample + "_variants_filtered.txt", sep="\t", index=False
        )

        print(">>> Filters:")
        print("\t>>> Quality Phred > 20")
        print("\t>>> Number of total reads > 5")

        all_var = len(final.index)
        filtered_var = len(df_filtered.index)
        percentage = round(filtered_var / all_var * 100, 2)

        print(
            "\nNumber of variants in all: ",
            all_var,
            "\nNumber of variants in filtered: ",
            filtered_var,
            "\nPercentage filtered variants: ",
            percentage,
            "%",
        )

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("formatting\n")


def steps_main(sample, toml_config, done):
    step = "main"

    command = ""

    # Quality control
    if "FastQC" not in done:
        print("\t>>> Quality control: FastQC (v0.12.1)")
        command += fastqc(sample, toml_config)

    # Trimming
    if toml_config["general"]["trimming"] == "True":
        print("\t>>> Trimming: BBDuk (v39.06)")
        if "BBDuk" not in done:
            command += bbduk(sample, toml_config)
    else:
        print("\t>>> Trimming: none")

    # Alignment
    if toml_config["general"]["alignment"] == "star":
        print("\t>>> Alignment: STAR (v2.7.11a)")
        if "STAR" not in done:
            command += star(sample, toml_config)
    elif toml_config["general"]["alignment"] == "bwa":
        print("\t>>> Alignment: BWA-MEM2 (v2.2.1)")
        if "BWA-MEM2" not in done:
            command += bwa(sample, toml_config)

    # Sorting and indexing
    print("\t>>> Sorting/Indexing: Samtools (v1.18)")
    if "Samtools" not in done:
        command += samtools(sample, toml_config)

    # Alignment QC
    print("\t>>> Quality control: FastQC (v0.12.1)")
    if "FastQC for bam" not in done:
        command += bamqc(sample, toml_config)

    print(command)

    cores = toml_config["general"]["threads"]
    memory = toml_config["general"]["memory"]
    time = toml_config["general"]["time"]
    email = toml_config["general"]["email"]
    output = toml_config["general"]["output"] + "/" + sample

    job = create_script(cores, memory, time, sample, step, email, command, output)

    with open(output + "/main.sh", "a") as f:
        f.write("\n# Main steps")
        f.write(f"\nmain=$(sbatch --parsable {job})\n")


def steps_downstream(sample, toml_config, done):
    step = "down"

    command = ""

    # Pseudo alignedment
    if toml_config["general"]["pseudo"] == "True":
        print("\t>>> Pseudoalignment: Salmon (v1.10.2)")
        if "Salmon" not in done:
            command += salmon(sample, toml_config)
    else:
        print("\t>>> Pseudoalignment: none")

    # MarkDuplicates
    # print("\t>>> MarkDuplicates: GATK (4.6.1.0) & Picard (v3.0.0)")
    # if "MarkDuplicates" not in done:
    #     command += markduplicates(sample, toml_config)

    # Quantification
    if toml_config["general"]["quantification"] == "True":
        print("\t>>> Quantification: featureCounts (v2.0.6)")
        if "FeatureCounts" not in done:
            command += featurecounts(sample, toml_config)
    else:
        print("\t>>> Quantification: none")

    # MultiQC
    print("\t>>> Quality control report: MultiQC (v1.31)")
    if "MultiQC" not in done:
        command += multiqc(sample, toml_config)

    # Variant Calling
    if toml_config["general"]["variants"] == "True":
        print("\t>>> Variant Calling: BCFtools (v1.22) & FreeBayes (v1.37)")
        if "BCFtools" not in done:
            command += bcftools(sample, toml_config)

        if "FreeBayes" not in done:
            command += freebayes(sample, toml_config)

        # Variant filtering
        print("\t>>> Variant Filtering: BCFtools (v1.18)")
        if "BCFtools filters" not in done:
            command += bcftools_filter(sample, toml_config)

        # Variant Annotation : SnpEff
        print("\t>>> Variant Annotation: SnpEff + SnpSift (v5.2a)")
        if "SnpEff" not in done:
            command += snpeff(sample, toml_config)

        # Variant formatting
        if (
            toml_config["general"]["reference"] == "grch37"
            or toml_config["general"]["reference"] == "grch38"
        ):
            if "dbNSFP" not in done:
                command += dbNSFP(sample, toml_config)
            if "formatting" not in done:
                command += formatting(sample, toml_config)
        else:
            if "formatting" not in done:
                command += formatting(sample, toml_config)
    else:
        print("\t>>> Variant Calling: none")

    print(command)

    cores = "1"
    memory = "4"
    time = toml_config["general"]["time"]
    email = toml_config["general"]["email"]
    output = toml_config["general"]["output"] + "/" + sample

    job = create_script(cores, memory, time, sample, step, email, command, output)

    with open(output + "/main.sh", "a") as f:
        f.write("\n# Main steps")
        if "main" not in done:
            f.write(f"\nsbatch --dependency=afterok:$main {job}\n")
        else:
            f.write(f"\nsbatch {job}\n")


if __name__ == "__main__":
    main()

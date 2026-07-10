import argparse
import re
import zipfile
import toml
import time
from datetime import datetime
import subprocess
import sys
import pandas as pd
import numpy as np
import os
from pathlib import Path

TOOL_PATH = "/lustre09/project/6019267/shared/tools/"
work_dir = os.getcwd()
username = os.environ.get("USER")


def main():
    start = get_time()

    parser = argparse.ArgumentParser(
        prog="PipelineShort",
        description="ShortReadSequencing pipeline for RNA or DNA data.",
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

    if "/" in toml_file:
        path_config = toml_file
    else:
        path_config = work_dir + "/" + toml_file

    with open(path_config, "r") as f:
        toml_config_initial = toml.load(f)

    if path_config != work_dir + "/config_final.toml":
        # Create final TOML config
        toml_config = create_config_final(path_config)

        if not args.test:
            print(
                "\n\n\n!!! WARNING !!!\nIf you want to change the parameters: Press CTRL+C now!\nModify config_final.toml and launch_pipeline with that config file.\n\nOtherwise it will run with default parameters.\n\n"
            )
            time.sleep(10)
    else:
        toml_config = toml_config_initial

    # Creating output and tmp directories for sample
    output = toml_config["general"]["output"] + "/" + sample
    tmp = toml_config["general"]["temporary"] + "/" + sample
    subprocess.run(["mkdir", "-p", output])
    subprocess.run(["mkdir", "-p", tmp])

    # Open file for steps done
    info = open(output + "/infos.txt", "a")
    start_str = ">>> {}-seq pipeline starting for {} at {}.".format(
        toml_config["general"]["sequencing"], sample, start
    )
    print(
        "=" * len(start_str) + "\n" + start_str + "\n" + "=" * len(start_str),
        file=sys.stdout,
    )
    print(f"\n>>> Output saved to {output}\n")

    info.close()

    done = []
    with open(output + "/infos.txt", "r") as f:
        for line in f:
            done.append(line.strip())

    with open(output + "/infos.txt", "a") as f:
        # Get tools and versions
        function_queue = []
        f.write(">>> Parameters:")

        genome = get_reference(toml_config["general"]["reference"], "")["fasta"]
        f.write(f"\t>>> Reference genome version: {genome}")

        # Quality control
        f.write("\t>>> Quality control: FastQC (v0.12.1)")
        if "FastQC" not in done:
            function_queue.append(fastqc)

        # Trimming
        if toml_config["general"]["trimming"] == "True":
            f.write("\t>>> Trimming: fastp (v1.0.1)")
            if "fastp" not in done:
                function_queue.append(fastp)
        else:
            f.write("\t>>> Trimming: none")

        # Alignment
        if toml_config["general"]["alignment"] == "star":
            f.write("\t>>> Alignment: STAR (v2.7.11a)")
            if "STAR" not in done:
                function_queue.append(star)
        elif toml_config["general"]["alignment"] == "bwa":
            f.write("\t>>> Alignment: BWA-MEM2 (v2.2.1)")
            if "BWA-MEM2" not in done:
                function_queue.append(bwa)
        else:
            f.write("\t>>> Alignment: none")

        # Pseudo alignedment
        if toml_config["general"]["pseudo"] == "True":
            f.write("\t>>> Pseudoalignment: Salmon (v1.10.2)")
            if "Salmon" not in done:
                function_queue.append(salmon)
        else:
            f.write("\t>>> Pseudoalignment: none")

        # Sorting and indexing
        f.write("\t>>> Sorting/Indexing: Samtools (v1.18)")
        if "Samtools" not in done:
            function_queue.append(samtools)

        # Alignment QC
        f.write("\t>>> Quality control: FastQC (v0.12.1)")
        if "FastQC for bam" not in done:
            function_queue.append(bamqc)

        # MarkDuplicates
        f.write("\t>>> MarkDuplicates: GATK (4.6.1.0) & Picard (v3.0.0)")
        if "MarkDuplicates" not in done:
            function_queue.append(markduplicates)

        # Quantification
        if toml_config["general"]["quantification"] == "True":
            f.write("\t>>> Quantification: featureCounts (v2.0.6)")
            if "FeatureCounts" not in done:
                function_queue.append(featurecounts)
        else:
            f.write("\t>>> Quantification: none")

        # MultiQC
        f.write("\t>>> Quality control report: MultiQC (v1.31)")
        if "MultiQC" not in done:
            function_queue.append(multiqc)

        # Variant Calling
        if toml_config["general"]["variants"] == "True":
            f.write("\t>>> Variant Calling: BCFtools (v1.22) & FreeBayes (v1.37)")
            if "BCFtools" not in done:
                function_queue.append(bcftools)

            if "FreeBayes" not in done:
                function_queue.append(freebayes)

            # Variant filtering
            f.write("\t>>> Variant Filtering: BCFtools (v1.18)")
            if "BCFtools filters" not in done:
                function_queue.append(bcftools_filter)

            # Variant Annotation
            if (
                toml_config["general"]["reference"] == "grch37"
                or toml_config["general"]["reference"] == "grch38"
            ):
                # VEP
                f.write("\t>>> Variant Annotation: VEP (v115.2)")
                if "vep" not in done:
                    function_queue.append(vep)

                # openCravat
                f.write("\t>>> Variant Annotation: openCravat (v2.17.0)")
                if "openCravat" not in done:
                    function_queue.append(openCravat)
            else:
                # SnpEff
                f.write("\t>>> Variant Annotation: SnpEff + SnpSift (v5.2a)")
                if "SnpEff" not in done:
                    function_queue.append(snpeff)
                if "formatting" not in done:
                    function_queue.append(formatting)
        else:
            f.write("\t>>> Variant Calling: none")

        # Create main.sh
        with open(f"{work_dir}/scripts/{sample}.sh", "w") as f:
            f.write("#!/bin/sh\n")
            f.write("\nDEPS=()\n")  # Add downstream dependencies

        # Calling each steps
        for func in function_queue:
            try:
                func(sample, toml_config, done)
            except Exception as e:
                print(f"Error: {e}")
                exit(1)

        # Check if in testing mode
        if args.test:
            print("\n\nTESTING MODE!\nThe pipeline will not be launched\n\n")
        else:
            # Call main.sh (Launch the pipeline)
            subprocess.run(["bash", f"{work_dir}/scripts/{sample}.sh"])


def create_config_final(path_config):
    with open(path_config, "r") as f:
        toml_config = toml.load(f)

    project = toml_config["general"]["project"]
    sequencing = toml_config["general"]["sequencing"]

    # set output in scratch
    toml_config["general"]["output"] = f"/lustre10/scratch/{username}/{project}/output"

    # set tmp in scratch
    toml_config["general"]["tmp"] = f"/lustre10/scratch/{username}/{project}/tmp"

    # set trimming, alignment, pseudoalignment, quantification and variants setting
    toml_config["general"]["trimming"] = "True"

    if sequencing == "rna":
        toml_config["general"]["alignment"] = "star"
        toml_config["star"]["outSAMtype1"] = "BAM"
        toml_config["star"]["outSAMtype2"] = "SortedByCoordinate"
        toml_config["star"]["twopassMode"] = "Basic"
        toml_config["star"]["outSJtype"] = "Standard"
        toml_config["star"]["quantMode"] = "GeneCounts"
        toml_config["fastp"]["phred"] = 15
        toml_config["fastp"]["length"] = 50
        toml_config["general"]["quantification"] = "True"
        toml_config["featurecounts"]["features"] = "gene"
        toml_config["featurecounts"]["attribute"] = "gene_id"
        toml_config["featurecounts"]["overlap"] = 1

    else:
        toml_config["general"]["alignment"] = "bwa"
        toml_config["general"]["quantification"] = "False"
        toml_config["fastp"]["phred"] = 20
        toml_config["fastp"]["length"] = 36
        toml_config["bwa-mem"]

    toml_config["general"]["pseudo"] = "False"
    toml_config["salmon"]["minScoreFraction"] = 0.65

    toml_config["general"]["variants"] = "True"

    with open(work_dir + "/config_final.toml", "w") as f:
        toml.dump(toml_config, f)

    return toml_config


def fill_template(tool, toml_config, sample, cpu, mem, time_allocated, env, command):
    email = toml_config["general"]["email"]
    job = f"{work_dir}/scripts/{tool}_{sample}.slurm"
    infos = toml_config["general"]["output"] + "/" + sample + "/infos.txt"

    with open(
        f"{TOOL_PATH}main_pipelines/short-read/ShortReadSequencing/template.txt",
        "r",
    ) as f:
        slurm = f.read()
        slurm_filled = slurm.format(
            cpu, mem, time_allocated, tool, sample, email, env, command, infos
        )

        with open(job, "w") as o:
            o.write(slurm_filled)

    return job


def get_time():
    now = datetime.now()
    return now


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


def get_file_trimmed(toml_config, output, sample):
    files = {}
    if toml_config["general"]["trimming"] == "True":
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


def fastqc(sample, toml_config, done):
    tool = "FastQC"

    cpu = ""
    mem = ""
    time_allocated = ""
    env = ""

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
            str(toml_config["general"]["threads"]),
            "--dir",
            temporary,
            "--kmers",
            str(toml_config["fastqc"]["kmers"]),
            Read1,
            Read2,
        ]
        command_str = " ".join(command)
        print(f">>> {command_str}\n")
        subprocess.run(command, check=True)

    else:
        Read = toml_config["general"]["fastq"] + "/" + sample + ".fastq.gz"

        command = [
            "fastqc",
            "-o",
            output,
            "--noextract",
            "--threads",
            str(toml_config["general"]["threads"]),
            "--dir",
            temporary,
            "--kmers",
            str(toml_config["fastqc"]["kmers"]),
            Read,
        ]

        command_str = " ".join(command)

    job = fill_template(
        tool, toml_config, sample, cpu, mem, time_allocated, env, command_str
    )

    if tool not in done:
        print(f"To-Do: {tool}")
        with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
            f.write(f"\n# Running {tool} for {sample}")
            f.write(f"\nDEPS+=($(sbatch --parsable {job}))\n")
    else:
        print(f"Done: {tool}")

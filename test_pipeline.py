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
        "--config", type=str, required=True, help="Project config file."
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
    tmp = toml_config["general"]["tmp"] + "/" + sample
    subprocess.run(["mkdir", "-p", output])
    subprocess.run(["mkdir", "-p", tmp])

    # Open file for steps done
    with open(output + "/infos.txt", "a") as f:
        start_str = ">>> {}-seq pipeline starting for {} at {}.".format(
            toml_config["general"]["sequencing"], sample, start
        )
        f.write(
            "\n" + "=" * len(start_str) + "\n" + start_str + "\n" + "=" * len(start_str)
        )
        f.write(f"\n>>> Output saved to {output}\n")

    done = []
    with open(output + "/infos.txt", "r") as f:
        for line in f:
            done.append(line.strip())

    with open(output + "/infos.txt", "a") as f:
        # Get tools and versions
        function_queue = []
        f.write(">>> Parameters:")

        genome = get_reference(toml_config["general"]["reference"], "")["fasta"]
        f.write(f"\n>>> Reference genome version: {genome}")

        # Quality control
        f.write("\n>>> Quality control: FastQC (v0.12.1)")
        if "FastQC" not in done:
            function_queue.append(fastqc)

        # Trimming
        if toml_config["general"]["trimming"] == "True":
            f.write("\n>>> Trimming: fastp (v1.0.1)")
            if "fastp" not in done:
                function_queue.append(fastp)
        else:
            f.write("\n>>> Trimming: none")

        # Alignment
        if toml_config["general"]["alignment"] == "star":
            f.write("\n>>> Alignment: STAR (v2.7.11a)")
            if "STAR" not in done:
                function_queue.append(star)
        elif toml_config["general"]["alignment"] == "bwa":
            f.write("\n>>> Alignment: BWA-MEM2 (v2.2.1)")
            if "BWA-MEM2" not in done:
                function_queue.append(bwa)
        else:
            f.write("\n>>> Alignment: none")

        # Pseudo alignedment
        if toml_config["general"]["pseudo"] == "True":
            f.write("\n>>> Pseudoalignment: Salmon (v1.10.2)")
            if "Salmon" not in done:
                function_queue.append(salmon)
        else:
            f.write("\n>>> Pseudoalignment: none")

        # Sorting and indexing
        f.write("\n>>> Sorting/Indexing: Samtools (v1.18)")
        if "Samtools" not in done:
            function_queue.append(samtools)

        # Alignment QC
        f.write("\n>>> Quality control: FastQC (v0.12.1)")
        if "FastQC_bam" not in done:
            function_queue.append(bamqc)

        # MarkDuplicates
        f.write("\n>>> MarkDuplicates: GATK (4.6.1.0) & Picard (v3.0.0)")
        if "MarkDuplicates" not in done:
            function_queue.append(markduplicates)

        # Quantification
        if toml_config["general"]["quantification"] == "True":
            f.write("\n>>> Quantification: featureCounts (v2.0.6)")
            if "FeatureCounts" not in done:
                function_queue.append(featurecounts)
        else:
            f.write("\n>>> Quantification: none")

        # MultiQC
        f.write("\n>>> Quality control report: MultiQC (v1.31)")
        if "MultiQC" not in done:
            function_queue.append(multiqc)

        # Variant Calling
        if toml_config["general"]["variants"] == "True":
            f.write("\n>>> Variant Calling: BCFtools (v1.22) & FreeBayes (v1.37)")
            if "BCFtools" not in done:
                function_queue.append(bcftools)

            if "FreeBayes" not in done:
                function_queue.append(freebayes)

            # Variant filtering
            f.write("\n>>> Variant Filtering: BCFtools (v1.18)")
            if "BCFtools_filters" not in done:
                function_queue.append(bcftools_filter)

            # Variant Annotation
            if (
                toml_config["general"]["reference"] == "grch37"
                or toml_config["general"]["reference"] == "grch38"
            ):
                # VEP
                f.write("\n>>> Variant Annotation: VEP (v115.2)")
                if "vep" not in done:
                    function_queue.append(vep)

                # openCravat
                f.write("\n>>> Variant Annotation: openCravat (v2.17.0)")
                if "openCravat" not in done:
                    function_queue.append(openCravat)
            else:
                # SnpEff
                f.write("\n>>> Variant Annotation: SnpEff + SnpSift (v5.2a)")
                if "SnpEff" not in done:
                    function_queue.append(snpeff)
        else:
            f.write("\n>>> Variant Calling: none")

        # Create main.sh
        subprocess.run(["mkdir", "-p", f"{work_dir}/scripts/"])
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

        cleanup(sample, toml_config, done, start)

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

    # set trimming, alignment, pseudoalignment, quantification and variants setting + other parameter
    toml_config["general"]["trimming"] = "True"

    if sequencing == "rna":
        toml_config["general"]["alignment"] = "star"
        toml_config["star"] = {}
        toml_config["star"]["outSAMtype1"] = "BAM"
        toml_config["star"]["outSAMtype2"] = "SortedByCoordinate"
        toml_config["star"]["twopassMode"] = "Basic"
        toml_config["star"]["outSJtype"] = "Standard"
        toml_config["star"]["quantMode"] = "GeneCounts"
        toml_config["fastp"] = {}
        toml_config["fastp"]["phred"] = 15
        toml_config["fastp"]["length"] = 50
        toml_config["general"]["quantification"] = "True"
        toml_config["featurecounts"] = {}
        toml_config["featurecounts"]["features"] = "gene"
        toml_config["featurecounts"]["attribute"] = "gene_id"
        toml_config["featurecounts"]["overlap"] = 1

    else:
        toml_config["general"]["alignment"] = "bwa"
        toml_config["general"]["quantification"] = "False"
        toml_config["fastp"] = {}
        toml_config["fastp"]["phred"] = 20
        toml_config["fastp"]["length"] = 36

    toml_config["general"]["pseudo"] = "False"
    toml_config["salmon"] = {}
    toml_config["salmon"]["minScoreFraction"] = 0.65

    toml_config["general"]["variants"] = "True"

    # FastQC
    toml_config["fastqc"] = {}
    toml_config["fastqc"]["kmers"] = 7

    # MarkDuplicates
    toml_config["markduplicates"] = {}
    toml_config["markduplicates"]["index"] = "true"
    toml_config["markduplicates"]["strategy"] = "SUM_OF_BASE_QUALITIES"
    toml_config["markduplicates"]["remove"] = "false"

    with open(work_dir + "/config_final.toml", "w") as f:
        toml.dump(toml_config, f)

    return toml_config


def fill_template(tool, toml_config, sample, cpu, mem, time_allocated, env, command):
    email = toml_config["general"]["email"]
    job = f"{work_dir}/scripts/{sample}_{tool}.slurm"
    infos = toml_config["general"]["output"] + "/" + sample + "/infos.txt"

    with open(
        f"{TOOL_PATH}main_pipelines/short-read/ShortReadSequencing/template.txt",
        "r",
    ) as f:
        slurm = f.read()
        slurm_filled = slurm.format(
            cpu, mem, time_allocated, sample, tool, email, env, command, infos
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

    cpu = "8"
    mem = "16"
    time_allocated = "00-01:00"
    env = "module load StdEnv/2023 python/3.11.5 fastqc/0.12.1\nsource /lustre09/project/6019267/shared/tools/main_pipelines/long-read/launch_pipeline_env/bin/activate"

    output = toml_config["general"]["output"] + "/" + sample + "/QC/fastQC"
    subprocess.run(["mkdir", "-p", output])

    temporary = toml_config["general"]["tmp"] + "/" + sample
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
            cpu,
            "--dir",
            temporary,
            "--kmers",
            str(toml_config["fastqc"]["kmers"]),
            Read1,
            Read2,
        ]

    else:
        Read = toml_config["general"]["fastq"] + "/" + sample + ".fastq.gz"

        command = [
            "fastqc",
            "-o",
            output,
            "--noextract",
            "--threads",
            cpu,
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


def fastp(sample, toml_config, done):
    tool = "fastp"

    cpu = "8"
    mem = "32"
    time_allocated = "00-01:00"
    env = "module load StdEnv/2023 fastp/1.0.1"

    output = toml_config["general"]["output"] + "/" + sample + "/Trimmed"
    subprocess.run(["mkdir", "-p", output])

    if toml_config["general"]["reads"] == "PE":
        I1 = toml_config["general"]["fastq"] + "/" + sample + "_R1.fastq.gz"
        I2 = toml_config["general"]["fastq"] + "/" + sample + "_R2.fastq.gz"
        O1 = output + "/" + sample + "_trimmed_R1.fastq.gz"
        O2 = output + "/" + sample + "_trimmed_R2.fastq.gz"

        command = [
            "fastp",
            "--in1",
            I1,
            "--in2",
            I2,
            "--out1",
            O1,
            "--out2",
            O2,
            "--thread",
            cpu,
            "--detect_adapter_for_pe",
            "--qualified_quality_phred",
            str(toml_config["fastp"]["phred"]),
            "--length_required",
            str(toml_config["fastp"]["length"]),
            "--html",
            output + "/" + sample + "_fastp_report.html",
            "--json",
            output + "/" + sample + "_fastp_report.json",
        ]
    else:
        I_SE = toml_config["general"]["fastq"] + "/" + sample + ".fastq.gz"
        O_SE = output + "/" + sample + "_trimmed.fastq.gz"
        command = [
            "fastp",
            "--in1",
            I_SE,
            "--out1",
            O_SE,
            "--thread",
            cpu,
            "--qualified_quality_phred",
            str(toml_config["fastp"]["phred"]),
            "--length_required",
            str(toml_config["fastp"]["length"]),
            "--html",
            output + "/" + sample + "_fastp_report.html",
            "--json",
            output + "/" + sample + "_fastp_report.json",
        ]

    command_str = " ".join(command)

    job = fill_template(
        tool, toml_config, sample, cpu, mem, time_allocated, env, command_str
    )

    if tool not in done:
        print(f"To-Do: {tool}")
        with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
            f.write(f"\n# Running {tool} for {sample}")
            f.write(f"\nfastp=$(sbatch --parsable {job})\n")
    else:
        print(f"Done: {tool}")


def star(sample, toml_config, done):
    tool = "STAR"

    cpu = "8"
    mem = "64"
    time_allocated = "00-11:00"
    env = "module load StdEnv/2023 star/2.7.11a"

    output = toml_config["general"]["output"] + "/" + sample + "/Aligned"
    subprocess.run(["mkdir", "-p", output])

    temporary = toml_config["general"]["tmp"] + "/" + sample + "/star_tmp"
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
            cpu,
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
            cpu,
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

    command_str1 = " ".join(command)
    command_str2 = (
        f"mv {output}/{sample}_Aligned.sortedByCoord.out.bam {output}/{sample}.bam"
    )
    command_str3 = f"mv {output}/{sample}_Log.final.out {output}/{sample}_summary_mapping_stats.out"
    command_str4 = (
        f"mv {output}/{sample}_Log.final.out {output}/{sample}_run_information.out"
    )
    command_str5 = f"rm {output}/{sample}_Log.progress.out"
    command_str6 = f"rm {output}/{sample}_Unmapped.out.mate1"
    command_str7 = f"rm {output}/{sample}_Unmapped.out.mate2"
    command_str8 = f"rm -r {output}/{sample}__STARpass1"
    command_str9 = f"mv {output}/{sample}__STARgenome/sjdbList.out.tab {output}/{sample}_sjdbList.out.tab"

    command_str = "\n".join(
        [
            command_str1,
            command_str2,
            command_str3,
            command_str4,
            command_str5,
            command_str6,
            command_str7,
            command_str8,
            command_str9,
        ]
    )

    job = fill_template(
        tool, toml_config, sample, cpu, mem, time_allocated, env, command_str
    )

    if tool not in done:
        print(f"To-Do: {tool}")
        if toml_config["general"]["trimming"] == "True":
            if "fastp" not in done:
                with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                    f.write(f"\n# Running {tool} for {sample}")
                    f.write(
                        f"\nstar=$(sbatch --parsable --dependency=afterok:$fastp {job})\n"
                    )
            else:
                with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                    f.write(f"\n# Running {tool} for {sample}")
                    f.write(f"\nstar=$(sbatch --parsable {job})\n")
        else:
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                f.write(f"\nstar=$(sbatch --parsable {job})\n")
    else:
        print(f"Done: {tool}")


def bwa(sample, toml_config, done):
    tool = "BWA-MEM2"

    cpu = "8"
    mem = "64"
    if toml_config["general"]["sequencing"] == "exome":
        time_allocated = "00-06:00"
    else:
        time_allocated = "00-23:00"
    env = "module load StdEnv/2023 bwa-mem2/2.2.1 samtools/1.22.1"

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
            "-t",
            cpu,
            "-v",
            "1",
            "-o",
            output + "/" + sample + ".sam",
            ref,
            I1_toAlign,
            I2_toAlign,
        ]
    else:
        command = [
            "bwa-mem2",
            "mem",
            "-t",
            cpu,
            "-v",
            "1",
            "-o",
            output + "/" + sample + ".sam",
            ref,
            I_toAlign,
        ]

    command_str1 = " ".join(command)
    command_str2 = f"samtools view -S --threads {cpu} -b {output}/{sample}.sam -o {output}/{sample}.bam"
    command_str3 = f"rm {output}/{sample}.sam"

    command_str = "\n".join([command_str1, command_str2, command_str3])

    job = fill_template(
        tool, toml_config, sample, cpu, mem, time_allocated, env, command_str
    )

    if tool not in done:
        print(f"To-Do: {tool}")
        if toml_config["general"]["trimming"] == "True":
            if "fastp" not in done:
                with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                    f.write(f"\n# Running {tool} for {sample}")
                    f.write(
                        f"\nbwa=$(sbatch --parsable --dependency=afterok:$fastp {job})\n"
                    )
            else:
                with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                    f.write(f"\n# Running {tool} for {sample}")
                    f.write(f"\nbwa=$(sbatch --parsable {job})\n")
        else:
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                f.write(f"\nbwa=$(sbatch --parsable {job})\n")
    else:
        print(f"Done: {tool}")


def salmon(sample, toml_config, done):
    tool = "Salmon"

    cpu = "8"
    mem = "32"
    if toml_config["general"]["sequencing"] == "exome":
        time_allocated = "00-06:00"
    else:
        time_allocated = "00-23:00"
    env = "module load StdEnv/2023 salmon/1.10.2"

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
            cpu,
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
            cpu,
            "--auxDir",
            "salmon_tmp",
            "--minScoreFraction",
            str(toml_config["salmon"]["minScoreFraction"]),
            "--unmatedReads",
            I_toAlign,
            "--output",
            output,
        ]

    command_str1 = " ".join(command)

    command_str2 = f"mv {output}/logs/salmon_quant.log {output}/{sample}_log.out"

    command_str3 = f"mv {output}/quant.sf {output}/{sample}_transcript_quant.sf"
    command_str4 = f"rm {output}/cmd_info.json"
    command_str5 = f"rm -r {output}/libParams"
    command_str6 = f"rm -r {output}/logs"
    command_str7 = f"rm -r {output}/salmon_tmp"

    command_str = "\n".join(
        [
            command_str1,
            command_str2,
            command_str3,
            command_str4,
            command_str5,
            command_str6,
            command_str7,
        ]
    )

    job = fill_template(
        tool, toml_config, sample, cpu, mem, time_allocated, env, command_str
    )

    if tool not in done:
        print(f"To-Do: {tool}")
        if toml_config["general"]["trimming"] == "True":
            if "fastp" not in done:
                with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                    f.write(f"\n# Running {tool} for {sample}")
                    f.write(
                        f"\nDEPS+=($(sbatch --parsable --dependency=afterok:$fastp {job}))\n"
                    )
            else:
                with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                    f.write(f"\n# Running {tool} for {sample}")
                    f.write(f"\nDEPS+=($(sbatch --parsable {job}))\n")
        else:
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                f.write(f"\nDEPS+=($(sbatch --parsable {job}))\n")
    else:
        print(f"Done: {tool}")


def samtools(sample, toml_config, done):
    tool = "Samtools"

    cpu = "8"
    mem = "32"
    if toml_config["general"]["sequencing"] == "rna":
        time_allocated = "00-06:00"
    elif toml_config["general"]["sequencing"] == "exome":
        time_allocated = "00-06:00"
    else:
        time_allocated = "00-11:00"
    env = "module load StdEnv/2023 samtools/1.22.1"

    in_out = toml_config["general"]["output"] + "/" + sample + "/Aligned"

    inBAM = in_out + "/" + sample + ".bam"
    bamCoord = in_out + "/" + sample + "_sortedCoordinate.bam"
    stats1 = in_out + "/" + sample + "_stats1.txt"
    stats2 = in_out + "/" + sample + "_stats2.txt"
    stats = in_out + "/" + sample + "_stats.txt"

    # Sort by coordinate
    command_str1 = f"samtools sort --threads {cpu} -m 4G {inBAM} -o {bamCoord}"

    # Index bam sorted by coordinates
    command_str2 = f"samtools index --threads {cpu} -b bamCoord -o {bamCoord}.bai"

    # alignment stats
    command_str3 = f"samtools stats {bamCoord} | grep ^SN | cut -f 2- > {stats}"
    # command_str4 = f"grep ^SN {stats1} > {stats2}"
    # command_str5 = f"cut -f 2- {stats2} > {stats}"

    # command_str6 = f"rm {stats1}"
    # command_str7 = f"rm {stats2}"
    command_str4 = f"rm {inBAM}"

    command_str = "\n".join(
        [
            command_str1,
            command_str2,
            command_str3,
            command_str4,
        ]
    )

    job = fill_template(
        tool, toml_config, sample, cpu, mem, time_allocated, env, command_str
    )

    if tool not in done:
        print(f"To-Do: {tool}")
        if (toml_config["general"]["sequencing"] == "RNA") & ("STAR" not in done):
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                f.write(
                    f"\nsamtools=$(sbatch --parsable --dependency=afterok:$star {job})\n"
                )
        elif (toml_config["general"]["sequencing"] in ["exome", "genome"]) & (
            "BWA-MEM2" not in done
        ):
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                f.write(
                    f"\nsamtools=$(sbatch --parsable --dependency=afterok:$bwa {job})\n"
                )
        else:
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                f.write(f"\nsamtools=$(sbatch --parsable {job})\n")
    else:
        print(f"Done: {tool}")


def bamqc(sample, toml_config, done):
    tool = "FastQC_bam"

    cpu = "8"
    mem = "16"
    time_allocated = "00-01:00"
    env = "module load StdEnv/2023 python/3.11.5 fastqc/0.12.1\nsource /lustre09/project/6019267/shared/tools/main_pipelines/long-read/launch_pipeline_env/bin/activate"

    output = toml_config["general"]["output"] + "/" + sample + "/QC/fastQC"
    subprocess.run(["mkdir", "-p", output])

    temporary = toml_config["general"]["tmp"] + "/" + sample
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
        cpu,
        "--dir",
        temporary,
        "--kmers",
        str(toml_config["fastqc"]["kmers"]),
        input,
    ]

    command_str = " ".join(command)

    job = fill_template(
        tool, toml_config, sample, cpu, mem, time_allocated, env, command_str
    )

    if tool not in done:
        print(f"To-Do: {tool}")
        if "Samtools" not in done:
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                f.write(
                    f"\nDEPS+=($(sbatch --parsable --dependency=afterok:$samtools {job}))\n"
                )
        else:
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                f.write(f"\nDEPS+=($(sbatch --parsable {job}))\n")
    else:
        print(f"Done: {tool}")


def markduplicates(sample, toml_config, done):
    tool = "MarkDuplicates"

    cpu = "4"
    mem = "16"
    time_allocated = "00-01:00"
    env = "module load StdEnv/2023 gatk/4.6.1.0"

    output = toml_config["general"]["output"] + "/" + sample + "/MarkDuplicates/"
    subprocess.run(["mkdir", "-p", output])

    temporary = toml_config["general"]["tmp"] + "/" + sample + "/md_tmp"
    subprocess.run(["mkdir", "-p", temporary])

    input = toml_config["general"]["output"] + "/" + sample + "/Aligned"
    bamCoord = input + "/" + sample + "_sortedCoordinate.bam"
    bam_RG = input + "/" + sample + "_sortedCoordinate_RG.bam"
    metrics = output + sample + "_duplicates_metrics.txt"
    records = output + sample + "_markDuplicates.bam"

    # Add RG tag to bam file
    command = [
        "samtools",
        "addreplacerg",
        "-r",
        f"@RG\\tID:{sample}\\tPL:Illumina\\tSM:{sample}\\tPU:{sample}",
        "-o",
        bam_RG,
        bamCoord,
    ]
    command_str1 = " ".join(command)

    command2 = [
        "gatk",
        "MarkDuplicates",
        "--VERBOSITY",
        "ERROR",
        "--INPUT",
        bam_RG,
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

    command_str2 = " ".join(command2)

    command_str3 = f"rm -r {temporary}"
    command_str4 = f"rm {bam_RG}"

    command_str = "\n".join(
        [
            command_str1,
            command_str2,
            command_str3,
            command_str4,
        ]
    )

    job = fill_template(
        tool, toml_config, sample, cpu, mem, time_allocated, env, command_str
    )

    if tool not in done:
        print(f"To-Do: {tool}")
        if "Samtools" not in done:
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                if toml_config["general"]["sequencing"] == "rna":
                    f.write(
                        f"\nmarkDuplicates=$(sbatch --parsable --dependency=afterok:$samtools {job})\n"
                    )
                else:
                    f.write(
                        f"\nDEPS+=($(sbatch --parsable --dependency=afterok:$samtools {job}))\n"
                    )
        else:
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                if toml_config["general"]["sequencing"] == "rna":
                    f.write(f"\nmarkDuplicates=$(sbatch --parsable {job})\n")
                else:
                    f.write(f"\nDEPS+=($(sbatch --parsable {job}))\n")
    else:
        print(f"Done: {tool}")


def featurecounts(sample, toml_config, done):
    tool = "FeatureCounts"

    cpu = "8"
    mem = "32"
    time_allocated = "00-01:00"
    env = "module load StdEnv/2023 subread/2.0.6"

    temporary = toml_config["general"]["tmp"] + "/" + sample
    output = toml_config["general"]["output"] + "/" + sample + "/FeatureCounts"
    subprocess.run(["mkdir", "-p", output])
    markduplicates_dir = (
        toml_config["general"]["output"] + "/" + sample + "/MarkDuplicates/"
    )

    input = markduplicates_dir + sample + "_markDuplicates.bam"

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
            cpu,
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
            cpu,
            "--tmpDir",
            temporary,
            "-a",
            gtf,
            "-o",
            output + "/" + sample + "_geneID.txt",
            input,
        ]

    command_str1 = " ".join(command)

    command_str2 = f"tail -n +2 {output}/{sample}_geneID.txt | cut -f1,7 | sed -i '1s|.*|gene_id\t{sample}|' > {output}/{sample}_counts.txt"

    command_str = "\n".join([command_str1, command_str2])

    job = fill_template(
        tool, toml_config, sample, cpu, mem, time_allocated, env, command_str
    )

    if tool not in done:
        print(f"To-Do: {tool}")
        if "MarkDuplicates" not in done:
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                f.write(
                    f"\nDEPS+=($(sbatch --parsable --dependency=afterok:$markDuplicates {job}))\n"
                )
        else:
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                f.write(f"\nDEPS+=($(sbatch --parsable {job}))\n")
    else:
        print(f"Done: {tool}")


def multiqc(sample, toml_config, done):
    tool = "MultiQC"

    cpu = "1"
    mem = "4"
    time_allocated = "00-01:00"
    env = "module load StdEnv/2023 python/3.11.5 apptainer/1.3.5"

    input = toml_config["general"]["output"] + "/" + sample + "/"
    output = toml_config["general"]["output"] + "/" + sample + "/QC/multiQC"
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

    command_str1 = f"apptainer run /lustre09/project/6019267/shared/tools/others/multiqc/multiqc.sif multiqc --file-list {output}/my_file_list.txt --force --filename {sample}_multiqc_report --outdir {output}/"
    command_str2 = f"rm {output}/my_file_list.txt"
    command_str3 = f"rm -r {output}/{sample}_multiqc_report_data"

    command_str = "\n".join([command_str1, command_str2, command_str3])

    job = fill_template(
        tool, toml_config, sample, cpu, mem, time_allocated, env, command_str
    )

    if tool not in done:
        print(f"To-Do: {tool}")
        with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
            f.write("\n# MultiQC")
            f.write('\nDEPENDENCY_LIST=$(IFS=:; echo "${DEPS[*]}")')
            f.write("\nif [ ${#DEPS[@]} -gt 0 ]; then")
            f.write(f"\n\tsbatch --dependency=afterok:$DEPENDENCY_LIST {job}")
            f.write("\nelse")
            f.write(f"\n\tsbatch {job}")
            f.write("\nfi\n")
    else:
        print(f"Done: {tool}")


def bcftools(sample, toml_config, done):
    tool = "BCFtools"

    cpu = "8"
    mem = "32"
    if toml_config["general"]["sequencing"] == "genome":
        time_allocated = "00-11:00"
    else:
        time_allocated = "00-03:00"
    env = "module load StdEnv/2023 bcftools/1.22"

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

    command_str = (
        f"bcftools mpileup --threads {cpu} -d 1000000 --max-idepth 1000000 -Ou -f {ref} {input} | "
        f"bcftools call --threads {cpu} -m -v -Ou | "
        f"bcftools reheader --samples {toml_config['general']['output']}/{sample}/sample.txt | "
        f"bcftools norm -f {ref} -m -any -Ou | "
        f"bcftools filter -i 'QUAL >= 10 && FORMAT/DP >= 5' -Oz -o {output}{sample}_bcftools.vcf.gz && "
        f"tabix -p vcf {output}{sample}_bcftools.vcf.gz"
    )

    job = fill_template(
        tool, toml_config, sample, cpu, mem, time_allocated, env, command_str
    )

    if tool not in done:
        print(f"To-Do: {tool}")
        if "Samtools" not in done:
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                f.write(
                    f"\nbcftools=$(sbatch --parsable --dependency=afterok:$samtools {job})\n"
                )
        else:
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                f.write(f"\nbcftools=$(sbatch --parsable {job})\n")
    else:
        print(f"Done: {tool}")


def freebayes(sample, toml_config, done):
    tool = "FreeBayes"

    cpu = "8"
    mem = "32"
    if toml_config["general"]["sequencing"] == "genome":
        time_allocated = "00-11:00"
    else:
        time_allocated = "00-03:00"
    env = "module load StdEnv/2023 freebayes/1.3.7 bcftools/1.22"

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

    command_str = (
        f"freebayes -f {ref} {input} | "
        f"bcftools reheader --samples {toml_config['general']['output']}/{sample}/sample.txt | "
        f"bcftools norm -f {ref} -m -any -Ou | "
        f"bcftools filter -i 'QUAL >= 10 && FORMAT/DP >= 5' -Oz -o {output}{sample}_freebayes.vcf.gz && "
        f"tabix -p vcf {output}{sample}_freebayes.vcf.gz"
    )

    job = fill_template(
        tool, toml_config, sample, cpu, mem, time_allocated, env, command_str
    )

    if tool not in done:
        print(f"To-Do: {tool}")
        if "Samtools" not in done:
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                f.write(
                    f"\nfreebayes=$(sbatch --parsable --dependency=afterok:$samtools {job})\n"
                )
        else:
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                f.write(f"\nfreebayes=$(sbatch --parsable {job})\n")
    else:
        print(f"Done: {tool}")


def bcftools_filter(sample, toml_config, done):
    tool = "BCFtools_filters"

    cpu = "4"
    mem = "16"
    if toml_config["general"]["sequencing"] == "genome":
        time_allocated = "00-11:00"
    else:
        time_allocated = "00-03:00"
    env = "module load StdEnv/2023 bcftools/1.22"

    output = toml_config["general"]["output"] + "/" + sample + "/Variants/"
    isec_dir = f"{output}isec_temp"

    subprocess.run(["mkdir", "-p", output])
    subprocess.run(["mkdir", "-p", isec_dir])

    command_str = (
        # 1. Split into temp directory (0000 = bcftools-only, 0001 = freebayes-only, 0002 = shared)
        f"bcftools isec -p {isec_dir} -O z "
        f"{output}{sample}_bcftools.vcf.gz {output}{sample}_freebayes.vcf.gz && "
        # 2. Concat all sites (-Ou) and pipe directly into final filter (-Oz)
        f"bcftools concat -a -Ou "
        f"{isec_dir}/0000.vcf.gz {isec_dir}/0001.vcf.gz {isec_dir}/0002.vcf.gz | "
        f"bcftools filter -i 'QUAL >= 10 && FORMAT/DP >= 5' -O v -o {output}{sample}_merged.vcf && "
        # 3. Cleanup temporary files
        f"rm -r {isec_dir} && "
        f"rm {toml_config['general']['output']}/{sample}/sample.txt"
    )

    job = fill_template(
        tool, toml_config, sample, cpu, mem, time_allocated, env, command_str
    )

    if tool not in done:
        print(f"To-Do: {tool}")
        active_deps = []
        if "bcftools" not in done:
            active_deps.append("$bcftools")
        if "freebayes" not in done:
            active_deps.append("$freebayes")

        if active_deps:
            dep_str = f"--dependency=afterok:{':'.join(active_deps)} "
        else:
            dep_str = ""
        with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
            f.write(f"\n# Running {tool} for {sample}")
            f.write(f"\nmerge=$(sbatch --parsable {dep_str}{job})\n")
    else:
        print(f"Done: {tool}")


def vep(sample, toml_config, done):
    tool = "VEP"

    cpu = "4"
    mem = "16"
    if toml_config["general"]["sequencing"] == "genome":
        time_allocated = "00-11:00"
    else:
        time_allocated = "00-03:00"
    env = "module load StdEnv/2023 apptainer/1.3.5"

    output = toml_config["general"]["output"] + "/" + sample + "/Variants/"
    vcf = f"{output}{sample}_merged.vcf"
    out_vcf = output + sample + "_vep.vcf"

    vep_cmd = [
        "apptainer",
        "run",
        "/lustre09/project/6019267/shared/tools/variants/annotation/VEP/vep.sif",
        "vep",
        "--offline",
        "--cache",
        "--dir_cache",
        "/lustre09/project/6019267/shared/tools/variants/annotation/VEP/vep_data/",
        "-i",
        vcf,
        "-o",
        out_vcf,
        "--vcf",
        "--force_overwrite",
        "--stats_file",
        output + sample + "_vep_summary.html",
        "--variant_class",
        "--domains",
        "--regulatory",
        "--canonical",
        "--individual_zyg",
        "all",
        "--hgvs",
        "--numbers",
        "--pubmed",
        "--uniprot",
        "--mane",
        "--protein",
        "--symbol",
        "--ccds",
        "--biotype",
        "--gene_phenotype",
        "--mirna",
        "--appris",
        "--tsl",
        "--pick",
        "--fork",
        cpu,
    ]

    command_str = " ".join(vep_cmd)

    job = fill_template(
        tool, toml_config, sample, cpu, mem, time_allocated, env, command_str
    )

    if tool not in done:
        print(f"To-Do: {tool}")
        if "BCFtools_filters" not in done:
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                f.write(
                    f"\nvep=$(sbatch --parsable --dependency=afterok:$merge {job})\n"
                )
        else:
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                f.write(f"\nvep=$(sbatch --parsable {job})\n")
    else:
        print(f"Done: {tool}")


def openCravat(sample, toml_config, done):
    tool = "openCravat"

    cpu = "8"
    mem = "128"
    if toml_config["general"]["sequencing"] == "genome":
        time_allocated = "01-11:00"
    else:
        time_allocated = "00-23:00"
    env = ""

    output = toml_config["general"]["output"] + "/" + sample + "/Variants/"
    vcf = output + sample + "_vep.vcf"

    genome = toml_config["general"]["reference"]

    if genome == "grch37":
        ref = "hg19"
    if genome == "grch38":
        ref = "hg38"

    oc = [
        "/lustre09/project/6019267/shared/tools/variants/annotation/openCravat_env/bin/oc",
        "run",
        vcf,
        "-a",
        "alfa",
        "alphamissense",
        "cadd",
        "clingen",
        "clinvar",
        "clinvar_acmg",
        "dbsnp",
        "denovo",
        "gerp",
        "geuvadis",
        "gnomad4",
        "gtex",
        "gwas_catalog",
        "hg19",
        "interpro",
        "mutation_assessor",
        "mutationtaster",
        "phastcons",
        "phylop",
        "polyphen2",
        "provean",
        "revel",
        "sift",
        "spliceai",
        "thousandgenomes",
        "aloft",
        "bayesdel",
        "clinpred",
        "fathmm",
        "genehancer",
        "hpo",
        "mitomap",
        "mutpred_indel",
        "omim",
        "pangolin",
        "repeat",
        "--mp",
        cpu,
        "--debug",
        "-l",
        ref,
        "-d",
        output,
        "-t",
        "text",
        "--cleanrun",
    ]

    command_str1 = " ".join(oc)

    command_str2 = f"python -u {TOOL_PATH}main_pipelines/short-read/ShortReadSequencing/oc_filtering.py --output {output} --sample {sample} --vcf {vcf}"

    command_str3 = (
        f"rm {output}{sample}_vep.vcf && "
        f"rm {output}{sample}_vep.vcf.sqlite && "
        f"rm {output}{sample}_vep.vcf.tsv && "
    )
    command_str = "\n".join([command_str1, command_str2, command_str3])

    job = fill_template(
        tool, toml_config, sample, cpu, mem, time_allocated, env, command_str
    )

    if tool not in done:
        print(f"To-Do: {tool}")
        if "VEP" not in done:
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                f.write(
                    f"\nDEPS+=($(sbatch --parsable --dependency=afterok:$vep {job}))\n"
                )
        else:
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                f.write(f"\nDEPS+=($(sbatch --parsable {job}))\n")
    else:
        print(f"Done: {tool}")


def snpeff(sample, toml_config, done):
    tool = "SnpEff"

    cpu = "4"
    mem = "32"
    time_allocated = "00-11:00"
    env = ""

    snpeff_path = (
        "/lustre09/project/6019267/shared/tools/main_pipelines/short-read/snpEff"
    )

    genome = toml_config["general"]["reference"]

    if genome == "grcz11":
        ref = "GRCz11.105"
    if genome == "grcm39":
        ref = "GRCm39.105"
    if genome == "wbcel235":
        ref = "WBcel235.105"

    path = toml_config["general"]["output"] + "/" + sample + "/Variants"

    command_str1 = (
        # 1. snpEff reads sample.vcf and streams output to stdout
        f"java -jar {snpeff_path}/snpEff.jar -noLog -c {snpeff_path}/snpEff.config "
        f"-stats {path}/{sample}_summary.html -csvStats {path}/{sample}_summary.csv "
        f"{ref} {path}/{sample}.vcf | "
        # 2. SnpSift varType reads from stdin (-) and streams output to stdout
        f"java -jar {snpeff_path}/SnpSift.jar varType -noLog - | "
        # 3. SnpSift extractFields reads from stdin (-) and writes TSV output to file
        f"java -jar {snpeff_path}/SnpSift.jar extractFields -noLog -s '|' -e '.' - "
        f"CHROM POS REF ALT QUAL HOM DP TYPE "
        f"'ANN[*].GENE' 'ANN[*].GENEID' 'ANN[*].FEATUREID' "
        f"'ANN[*].EFFECT' 'ANN[*].IMPACT' 'ANN[*].BIOTYPE' "
        f"'ANN[*].HGVS_C' 'ANN[*].HGVS_P' "
        f"> {path}/{sample}_annotated.txt"
    )

    command_str2 = f"python -u {TOOL_PATH}main_pipelines/short-read/ShortReadSequencing/snpeff_filtering.py --path {path} --sample {sample}"

    command_str = "\n".join([command_str1, command_str2])

    job = fill_template(
        tool, toml_config, sample, cpu, mem, time_allocated, env, command_str
    )

    if tool not in done:
        print(f"To-Do: {tool}")
        if "BCFtools_filters" not in done:
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                f.write(
                    f"\nDEPS+=($(sbatch --parsable --dependency=afterok:$merge {job}))\n"
                )
        else:
            with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
                f.write(f"\n# Running {tool} for {sample}")
                f.write(f"\nDEPS+=($(sbatch --parsable {job}))\n")
    else:
        print(f"Done: {tool}")


def cleanup(sample, toml_config, done, start):
    tool = "cleanup"

    cpu = "1"
    mem = "4"
    time_allocated = "00-03:00"
    env = ""

    command_str = f"python -u {TOOL_PATH}main_pipelines/short-read/ShortReadSequencing/cleanup.py --toml_config {toml_config} --start {start} --sample {sample}"

    job = fill_template(
        tool, toml_config, sample, cpu, mem, time_allocated, env, command_str
    )

    if tool not in done:
        print(f"To-Do: {tool}")
        with open(f"{work_dir}/scripts/{sample}.sh", "a") as f:
            f.write("\n# Cleanup")
            f.write('\nDEPENDENCY_LIST=$(IFS=:; echo "${DEPS[*]}")')
            f.write("\nif [ ${#DEPS[@]} -gt 0 ]; then")
            f.write(f"\n\tsbatch --dependency=afterok:$DEPENDENCY_LIST {job}")
            f.write("\nelse")
            f.write(f"\n\tsbatch {job}")
            f.write("\nfi\n")
    else:
        print(f"Done: {tool}")


if __name__ == "__main__":
    main()

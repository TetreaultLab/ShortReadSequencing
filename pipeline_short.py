import argparse
import re
import zipfile
import toml
from datetime import datetime
import subprocess
import sys
import pandas as pd
import numpy as np
import os
from pathlib import Path


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

    args = parser.parse_args()

    # Sample name
    sample = args.sample

    # Loading TOML config
    with open(args.config, "r") as f:
        toml_config = toml.load(f)

    # Start of pipeline
    start_str = ">>> {}-seq pipeline starting for {} at {}.".format(
        toml_config["general"]["sequencing"], sample, start
    )
    print(
        "=" * len(start_str) + "\n" + start_str + "\n" + "=" * len(start_str),
        file=sys.stdout,
    )

    # Creating output and tmp directories for sample
    output = toml_config["general"]["output"] + "/" + sample
    tmp = toml_config["general"]["temporary"] + "/" + sample
    subprocess.run(["mkdir", "-p", output])
    subprocess.run(["mkdir", "-p", tmp])
    print(f"\n>>> Output saved to {output}\n")

    # Open file for steps done
    steps = open(output + "/steps_done.txt", "a")
    steps.write("\nLoading ENV\n")
    steps.close()

    done = []
    with open(output + "/steps_done.txt", "r") as f:
        for line in f:
            done.append(line.strip())

    # Get tools and versions
    function_queue = []
    print(">>> Parameters:")

    genome = get_reference(toml_config["general"]["reference"], "")["fasta"]
    print(f"\t>>> Reference genome version: {genome}")

    # Quality control
    print("\t>>> Quality control: FastQC (v0.12.1)")
    if "FastQC" not in done:
        function_queue.append(fastqc)

    # Trimming
    if toml_config["general"]["trimming"] == "True":
        print("\t>>> Trimming: BBDuk (v39.06)")
        if "BBDuk" not in done:
            function_queue.append(bbduk)
    else:
        print("\t>>> Trimming: none")

    # Alignment
    if toml_config["general"]["alignment"] == "star":
        print("\t>>> Alignment: STAR (v2.7.11a)")
        if "STAR" not in done:
            function_queue.append(star)
    elif toml_config["general"]["alignment"] == "bwa":
        print("\t>>> Alignment: BWA-MEM2 (v2.2.1)")
        if "BWA-MEM2" not in done:
            function_queue.append(bwa)
    else:
        print("\t>>> Alignment: none")

    # Pseudo alignedment
    if toml_config["general"]["pseudo"] == "True":
        print("\t>>> Pseudoalignment: Salmon (v1.10.2)")
        if "Salmon" not in done:
            function_queue.append(salmon)
    else:
        print("\t>>> Pseudoalignment: none")

    # Sorting and indexing
    print("\t>>> Sorting/Indexing: Samtools (v1.18)")
    if "Samtools" not in done:
        function_queue.append(samtools)

    # Alignment QC
    print("\t>>> Quality control: FastQC (v0.12.1)")
    if "FastQC for bam" not in done:
        function_queue.append(bamqc)

    # MarkDuplicates
    # print("\t>>> MarkDuplicates: GATK (4.6.1.0) & Picard (v3.0.0)")
    # if "MarkDuplicates" not in done:
    #     function_queue.append(markduplicates)

    # Quantification
    if toml_config["general"]["quantification"] == "True":
        print("\t>>> Quantification: featureCounts (v2.0.6)")
        if "FeatureCounts" not in done:
            function_queue.append(featurecounts)
    else:
        print("\t>>> Quantification: none")

    # MultiQC
    print("\t>>> Quality control report: MultiQC (v1.31)")
    if "MultiQC" not in done:
        function_queue.append(multiqc)

    # Variant Calling
    if toml_config["general"]["variants"] == "True":
        print("\t>>> Variant Calling: BCFtools (v1.22) & FreeBayes (v1.37)")
        if "BCFtools" not in done:
            function_queue.append(bcftools)

        if "FreeBayes" not in done:
            function_queue.append(freebayes)

        # Variant filtering
        print("\t>>> Variant Filtering: BCFtools (v1.18)")
        if "BCFtools filters" not in done:
            function_queue.append(bcftools_filter)

        # Variant Annotation
        if (
            toml_config["general"]["reference"] == "grch37"
            or toml_config["general"]["reference"] == "grch38"
        ):
            # openCravat
            print("\t>>> Variant Annotation: openCravat (v2.17.0)")
            if "openCravat" not in done:
                function_queue.append(openCravat)
        else:
            # SnpEff
            print("\t>>> Variant Annotation: SnpEff + SnpSift (v5.2a)")
            if "SnpEff" not in done:
                function_queue.append(snpeff)
            if "formatting" not in done:
                function_queue.append(formatting)
    else:
        print("\t>>> Variant Calling: none")

    # Calling each steps
    for func in function_queue:
        try:
            func(sample, toml_config)
        except Exception as e:
            print(f"Error: {e}")
            exit(1)

    # Get log and extract all >>> lines from it and save to file "steps_summary.txt"
    # TO-DO

    # Move results to projects directory
    print("\n>>> Transferring from scratch to results")
    current_directory = os.getcwd()
    results = Path(str(current_directory).replace("/work", "/results/"))
    output = Path(toml_config["general"]["output"] + "/" + sample)
    results.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["rsync", "-avxH", "--no-g", "--no-p", "--partial", str(output), str(results)]
    subprocess.run(cmd, check=True)
    final = f"{str(results)}/{sample}"
    print(f"\n\n>>> Results found in {final}")

    end = get_time()
    total_time = end - start
    end_str = "\n>>> {}-seq pipeline for {} completed in {}.".format(
        toml_config["general"]["sequencing"], sample, total_time
    )
    print("=" * len(end_str) + "\n" + end_str + "\n" + "=" * len(end_str))


def get_time():
    now = datetime.now()
    return now


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


def check_fastqc_report(zip_path: Path):
    with zipfile.ZipFile(zip_path, "r") as z:
        file_list = z.namelist()
        summary_file = [f for f in file_list if f.endswith("summary.txt")][0]
        try:
            with z.open(summary_file) as f:
                content = f.read().decode("utf-8")
                print(content.rstrip())
        except KeyError:
            print(f"No summary.txt found in {zip_path}")
            return False

        # Search for WARN or FAIL
        if re.search(r"\b(FAIL)\b", content):
            print(
                ">>> FastQC failures detected.\n>>> Please check FastQC report (and adjust for trimming if necessary) before resubmitting (with --redo).\n"
            )
            # sys.exit(1)
        elif re.search(r"\b(WARN)\b", content):
            print(">>> FastQC warnings detected\n")
        else:
            print(">>> FastQC passed with no warnings or failures.\n")


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
            "6",
            "--dir",
            temporary,
            "--kmers",
            str(toml_config["fastqc"]["kmers"]),
            Read,
        ]

        command_str = " ".join(command)
        print(f">>> {command_str}\n")
        subprocess.run(command, check=True)

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("FastQC\n")

    # Check if any failures in FastQC report
    fastqc_reports = list(Path(output + "/").glob(sample + "*_fastqc.zip"))
    print()

    for report in fastqc_reports:
        check_fastqc_report(report)


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
    print(f">>> {command_str}\n")
    subprocess.run(command, check=True)

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("BBDuk\n")


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
    print(f">>> {command_str}\n")
    subprocess.run(command, check=True)

    # rename BAM file
    subprocess.run(
        [
            "mv",
            output + "/" + sample + "_Aligned.sortedByCoord.out.bam",
            output + "/" + sample + ".bam",
        ]
    )
    subprocess.run(
        [
            "mv",
            output + "/" + sample + "_Log.final.out",
            output + "/" + sample + "_summary_mapping_stats.out",
        ]
    )
    subprocess.run(
        [
            "mv",
            output + "/" + sample + "_Log.out",
            output + "/" + sample + "_run_information.out",
        ]
    )
    subprocess.run(["rm", output + "/" + sample + "_Log.progress.out"])
    subprocess.run(["rm", output + "/" + sample + "_Unmapped.out.mate1"])
    subprocess.run(["rm", output + "/" + sample + "_Unmapped.out.mate2"])
    subprocess.run(["rm", "-r", output + "/" + sample + "__STARpass1"])
    subprocess.run(
        [
            "mv",
            output + "/" + sample + "__STARgenome/sjdbList.out.tab",
            output + "/" + sample + "_sjdbList.out.tab",
        ]
    )

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("STAR\n")


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
    print(f">>> {command_str}\n")
    subprocess.run(command, check=True)

    # Change SAM to BAM format
    subprocess.run(
        [
            "samtools",
            "view",
            "-S",
            "-b",
            output + "/" + sample + ".sam",
            "-o",
            output + "/" + sample + ".bam",
        ],
        check=True,
    )
    # Remove the SAM file using subprocess
    sam_file = output + "/" + sample + ".sam"
    subprocess.run(["rm", sam_file], check=True)

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("BWA-MEM2\n")


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
    print(f">>> {command_str}\n")
    subprocess.run(command, check=True)

    subprocess.run(
        ["mv", output + "/logs/salmon_quant.log", output + "/" + sample + "_log.out"]
    )
    subprocess.run(
        ["mv", output + "/quant.sf", output + "/" + sample + "_transcript_quant.sf"]
    )
    subprocess.run(["rm", output + "/cmd_info.json"])
    subprocess.run(["rm", "-r", output + "/libParams"])
    subprocess.run(["rm", "-r", output + "/logs"])
    subprocess.run(["rm", "-r", output + "/salmon_tmp"])

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("Salmon\n")


def samtools(sample, toml_config):
    title("Samtools")
    in_out = toml_config["general"]["output"] + "/" + sample + "/Aligned"

    inBAM = in_out + "/" + sample + ".bam"
    bamCoord = in_out + "/" + sample + "_sortedCoordinate.bam"
    stats1 = in_out + "/" + sample + "_stats1.txt"
    stats2 = in_out + "/" + sample + "_stats2.txt"
    stats = in_out + "/" + sample + "_stats.txt"

    # Sort by coordinate
    subprocess.run(["samtools", "sort", inBAM, "-o", bamCoord])

    # Index bam sorted by coordinates
    subprocess.run(["samtools", "index", "-b", bamCoord, "-o", bamCoord + ".bai"])

    # alignment stats
    with open(stats1, "w") as outfile1:
        subprocess.run(["samtools", "stats", bamCoord], stdout=outfile1)

    with open(stats2, "w") as outfile2:
        subprocess.run(
            [
                "grep",
                "^SN",
                stats1,
            ],
            stdout=outfile2,
        )
    with open(stats, "w") as outfile:
        subprocess.run(["cut", "-f", "2-", stats2], stdout=outfile)

    subprocess.run(["rm", stats1])
    subprocess.run(["rm", stats2])
    subprocess.run(["rm", inBAM])

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("Samtools\n")


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
    print(f">>> {command_str}\n")
    subprocess.run(command, check=True)

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("FastQC for bam\n")

    # Check if any failures in FastQC report
    fastqc_reports = list(
        Path(output + "/").glob(sample + "*_sortedCoordinate_fastqc.zip")
    )
    print()

    for report in fastqc_reports:
        check_fastqc_report(report)


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

    command = [
        "freebayes",
        "-f",
        ref,
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
        output + sample + ".vcf",
        output + sample + "_merged.vcf.gz",
    ]

    command_filter_str = " ".join(command_filter)
    print(f">>> {command_filter_str}\n")
    subprocess.run(command_filter, check=True)

    # remove intermediate vcf files
    subprocess.run(["rm", output + sample + "_merged.vcf.gz"])
    subprocess.run(["rm", output + sample + "_freebayes.vcf"])
    subprocess.run(["rm", output + sample + "_bcftools.vcf"])

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("BCFtools filters\n")


def openCravat(sample, toml_config):
    title("Annotate with openCravat")

    threads = toml_config["general"]["threads"]
    output = toml_config["general"]["output"] + "/" + sample + "/Variants/"
    vcf = output + sample + ".vcf"

    genome = toml_config["general"]["reference"]

    if genome == "grch37":
        ref = "hg19"
    if genome == "grch38":
        ref = "hg38"

    oc = [
        "oc",
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
        "segway_blood",
        "segway_brain",
        "segway_muscle",
        "gencode",
        "calibrated_classification",
        "pathogenic",
        "rare_coding",
        "splicing",
        "--mp",
        str(threads),
        "--debug",
        "-l",
        ref,
        "-d",
        output,
        "-t",
        "text",
        "--cleanrun",
    ]

    command_str = " ".join(oc)
    print(f">>> {command_str}\n")

    current_directory = os.getcwd()
    with open(
        "/lustre09/project/6019267/shared/tools/main_pipelines/short-read/ShortReadSequencing/run_openCravat.bash",
        "r",
    ) as f:
        slurm = f.read()
        slurm_filled = slurm.format(command_str)

        with open(current_directory + "/run_openCravat_" + sample + ".bash", "w") as o:
            o.write(slurm_filled)
    subprocess.run(
        ["bash", current_directory + "/run_openCravat_" + sample + ".bash"],
        check=True,
    )

    def normalize_headers(h1, h2):
        # propagate h1 values forward
        last = ""
        h1_norm = []
        for x in h1:
            if x != "":
                last = x
            h1_norm.append(last)

        # combine
        cols = [f"{a}_{b}".replace(" ", "_") for a, b in zip(h1_norm, h2)]
        return cols

    def load_tsv_until_next_hash(path):
        with open(path, "r") as f:
            # Skip # metadata line(s)
            for line in f:
                if not line.startswith("#"):
                    header1 = line.rstrip("\n").split("\t")
                    break

            header2 = next(f).rstrip("\n").split("\t")

            # normalize headers so columns match correctly
            columns = normalize_headers(header1, header2)

            expected = len(columns)

            rows = []
            app = rows.append

            for line in f:
                if line.startswith("#"):
                    break
                fields = line.rstrip("\n").split("\t")

                if len(fields) != expected:
                    # Skip malformed row
                    print("malformed line: " + fields)

                app(fields)

        return pd.DataFrame(rows, columns=columns, dtype=str)

    df = load_tsv_until_next_hash(vcf + ".tsv")

    df.to_csv(output + sample + "_all.tsv", sep="\t", index=False, lineterminator="\n")

    cols = [
        "Variant_Annotation_Chrom",
        "Variant_Annotation_Position",
        "Variant_Annotation_End_Position",
        "Variant_Annotation_Ref_Base",
        "Variant_Annotation_Alt_Base",
        "VCF_Info_Phred",
        "VCF_Info_Zygosity",
        "VCF_Info_Alternate_reads",
        "VCF_Info_Total_reads",
        "Variant_Annotation_Gene",
        "Variant_Annotation_Transcript",
        "Variant_Annotation_Coding",
        "Variant_Annotation_Sequence_Ontology",
        "Variant_Annotation_Exon_Number",
        "Variant_Annotation_cDNA_change",
        "Variant_Annotation_Protein_Change",
        "Variant_Annotation_All_Mappings",
        "dbSNP_rsID",
        "1000_Genomes_AF",
        "ALFA:_Allele_Frequency_Aggregator_Global_AF",
        "gnomAD4_Global_AF",
        "gnomAD4_Global_AN",
        "gnomAD4_Global_AC",
        "gnomAD4_Non-Fin_Eur_AF",
        "gnomAD4_Non-Fin_Eur_AN",
        "gnomAD4_Non-Fin_Eur_AC",
        "CADD_Phred",
        "PolyPhen-2_HDIV_Rank_Score",
        "PolyPhen-2_HVAR_Rank_Score",
        "REVEL_Rank_Score",
        "SIFT_Rank_Score",
        "SIFT_Prediction",
        "SIFT_Confidence",
        "AlphaMissense_Protein",
        "AlphaMissense_Score",
        "AlphaMissense_Class",
        "GERP++_RS_Score",
        "Phast_Cons_Vert_Ranked_Score",
        "PhyloP_Vert_Ranked_Score",
        "BayesDel_Prediction_(AF)",
        "Mutation_Assessor_Functional_Impact",
        "MutationTaster_Prediction",
        "PROVEAN:_Protein_Variant_Effect_Analyzer_UniProt_Accession_Number",
        "PROVEAN:_Protein_Variant_Effect_Analyzer_Prediction",
        "ALoFT_Classification",
        "ALoFT_Confidence",
        "GeneHancer_GeneHancer_Type",
        "GeneHancer_Gene_Targets",
        "MutPred-Indel_Rank_score",
        "MutPred-Indel_Property",
        "ClinGen_Gene_Disease",
        "ClinVar_Clinical_Significance",
        "GTEx_eQTLs_Target_Gene",
        "GWAS_Catalog_Disease/Trait",
        "GWAS_Catalog_Odds_Ratio/Beta_Coeff",
        "GWAS_Catalog_P-value",
        "GWAS_Catalog_PMID",
        "OMIM_Entry_ID",
        "Pangolin_Splicing_Gain_Δscore",
        "Pangolin_Splicing_Loss_Δscore",
        "SpliceAI_Acceptor_Gain_Score",
        "SpliceAI_Acceptor_Loss_Score",
        "SpliceAI_Donor_Gain_Score",
        "SpliceAI_Donor_Loss_Score",
        "MITOMAP_Disease",
        "MITOMAP_MitoTip_Score",
        "MITOMAP_PubMed_ID",
        "Variant_Annotation_Samples",
    ]
    df_cols = df[cols]

    column_mapping = {
        "Variant_Annotation_Chrom": "chrom",
        "Variant_Annotation_Position": "start",
        "Variant_Annotation_End_Position": "end",
        "Variant_Annotation_Ref_Base": "ref",
        "Variant_Annotation_Alt_Base": "alt",
        "VCF_Info_Phred": "phred",
        "VCF_Info_Alternate_reads": "alt_reads",
        "VCF_Info_Total_reads": "total_reads",
        "VCF_Info_Zygosity": "zygosity",
        "Variant_Annotation_Gene": "gene",
        "Variant_Annotation_Transcript": "transcript",
        "Variant_Annotation_Coding": "coding",
        "Variant_Annotation_Sequence_Ontology": "sequence_ontology",
        "Variant_Annotation_Exon_Number": "exon_number",
        "Variant_Annotation_cDNA_change": "cDNA_change",
        "Variant_Annotation_Protein_Change": "protein_change",
        "Variant_Annotation_All_Mappings": "all_mappings",
        "dbSNP_rsID": "rsID",
        "1000_Genomes_AF": "1000_Genomes_AF",
        "ALFA:_Allele_Frequency_Aggregator_Global_AF": "ALFA_AF",
        "gnomAD4_Global_AF": "gnomAD4_global_AF",
        "gnomAD4_Global_AC": "gnomAD4_global_AC",
        "gnomAD4_Global_AN": "gnomAD4_global_AN",
        "gnomAD4_Non-Fin_Eur_AF": "gnomAD4_NFE_AF",
        "gnomAD4_Non-Fin_Eur_AC": "gnomAD4_NFE_AC",
        "gnomAD4_Non-Fin_Eur_AN": "gnomAD4_NFE_AN",
        "CADD_Phred": "CADD_score",
        "PolyPhen-2_HDIV_Rank_Score": "polyPhen2_HDIV_score",
        "PolyPhen-2_HVAR_Rank_Score": "polyPhen2_HVAR_score",
        "REVEL_Rank_Score": "REVEL_score",
        "SIFT_Rank_Score": "SIFT_score",
        "SIFT_Prediction": "SIFT_prediction",
        "SIFT_Confidence": "SIFT_confidence",
        "AlphaMissense_Protein": "AlphaMissense_protein",
        "AlphaMissense_Score": "AlphaMissense_score",
        "AlphaMissense_Class": "AlphaMissense_classification",
        "GERP++_RS_Score": "GERP++_score",
        "Phast_Cons_Vert_Ranked_Score": "PhastCons_score",
        "PhyloP_Vert_Ranked_Score": "PhyloP_score",
        "BayesDel_Prediction_(AF)": "BayesDel_prediction",
        "Mutation_Assessor_Functional_Impact": "MutationAssessor_functional_impact",
        "MutationTaster_Prediction": "MutationTaster_prediction",
        "PROVEAN:_Protein_Variant_Effect_Analyzer_UniProt_Accession_Number": "PROVEAN_uniprot",
        "PROVEAN:_Protein_Variant_Effect_Analyzer_Prediction": "PROVEAN_prediction",
        "ALoFT_Classification": "ALoFT_classification",
        "ALoFT_Confidence": "ALoFT_confidence",
        "GeneHancer_GeneHancer_Type": "GeneHancer_type",
        "GeneHancer_Gene_Targets": "GeneHancer_gene_targets",
        "MutPred-Indel_Rank_score": "MutPred-Indel_score",
        "MutPred-Indel_Property": "MutPred-Indel_proprety",
        "ClinGen_Gene_Disease": "ClinGen_gene_disease",
        "ClinVar_Clinical_Significance": "ClinVar_clinical_significance",
        "GTEx_eQTLs_Target_Gene": "GTEx_eQTLs_target_genes",
        "GWAS_Catalog_Disease/Trait": "GWAS_Catalog_trait",
        "GWAS_Catalog_Odds_Ratio/Beta_Coeff": "GWAS_Catalog_odds_ratio",
        "GWAS_Catalog_P-value": "GWAS_Catalog_pvalue",
        "GWAS_Catalog_PMID": "GWAS_Catalog_PMID",
        "OMIM_Entry_ID": "OMIM",
        "Pangolin_Splicing_Gain_Δscore": "Pangolin_splicing_gain_deltascore",
        "Pangolin_Splicing_Loss_Δscore": "Pangolin_splicing_loss_deltascore",
        "SpliceAI_Acceptor_Gain_Score": "SpliceAI_acceptor_gain_score",
        "SpliceAI_Acceptor_Loss_Score": "SpliceAI_acceptor_loss_score",
        "SpliceAI_Donor_Gain_Score": "SpliceAI_donor_gain_score",
        "SpliceAI_Donor_Loss_Score": "SpliceAI_donor_loss_score",
        "MITOMAP_Disease": "MITOMAP_disease",
        "MITOMAP_MitoTip_Score": "MITOMAP_MitoTip_score",
        "MITOMAP_PubMed_ID": "MITOMAP_PubMed_id",
        "Variant_Annotation_Samples": "sample",
    }

    df_small = df_cols.rename(columns=column_mapping)

    df_small["alt_reads"] = pd.to_numeric(df_small["alt_reads"], errors="coerce")
    df_small["total_reads"] = pd.to_numeric(df_small["total_reads"], errors="coerce")
    df_small["phred"] = pd.to_numeric(df_small["phred"], errors="coerce")

    df_filt = df_small[
        (df_small["alt_reads"] >= 2)
        & (df_small["total_reads"] >= 5)
        & (df_small["phred"] >= 20)
    ]

    print(">>> Filters:")
    print("\t>>> Quality Phred > 20")
    print("\t>>> Number of alt reads > 5")
    print("\t>>> Number of total reads > 5")

    all_var = len(df_cols.index)
    filtered_var = len(df_filt.index)
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

    df_filt.to_csv(
        output + sample + "_filtered.tsv", sep="\t", index=False, lineterminator="\n"
    )

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("openCravat\n")


def snpeff(sample, toml_config):
    title("SnpEff")

    snpeff = "/lustre09/project/6019267/shared/tools/main_pipelines/short-read/snpEff"

    genome = toml_config["general"]["reference"]

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
        path + "/" + sample + ".vcf",
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

    subprocess.run(["rm", path + "/" + sample + "_summary.csv"])
    subprocess.run(["rm", path + "/" + sample + "_snpeff.vcf"])

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("SnpEff\n")


def formatting(sample, toml_config):
    title("Format variants output")

    path = toml_config["general"]["output"] + "/" + sample + "/Variants"

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
                    "|".join([info[i] for info in infos]) for i in range(len(infos[0]))
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
    final.loc[(final["het_count"] > 1) & (final["Zygosity"] == "Het"), "Zygosity"] = (
        "Multiple-het"
    )

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


if __name__ == "__main__":
    main()

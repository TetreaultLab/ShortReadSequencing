import argparse
import toml
from datetime import datetime
import subprocess
import sys


def main():
    start = get_time()

    parser = argparse.ArgumentParser(
        prog="PipelineShort",
        description="ShortReadSequencing pipeline for RNA, Exome and Genome data.",
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
    print(f"\n>>> Output saved to '{output}'\n")

    # Get tools and versions
    dict_keys = toml_config.keys()

    function_queue = []
    print(">>> Parameters:")
    # Quality control
    vqc = subprocess.check_output(["cat", "../version_fastqc_N.txt"], text=True).strip()
    print(f"\t>>> Quality control: FastQC ({vqc})")
    function_queue.append(fastqc)

    # Trimming
    if "bbduk" in dict_keys:
        vt = subprocess.check_output(
            ["cat", "../version_bbduk_N.txt"], text=True
        ).strip()
        print(f"\t>>> Trimming: BBDuk (v{vt})")
        function_queue.append(bbduk)
    else:
        print("\t>>> Trimming: none")

    # Alignment
    if "star" in dict_keys:
        va = subprocess.check_output(
            ["cat", "../version_star_N.txt"], text=True
        ).strip()
        print(f"\t>>> Alignment: STAR (v{va})")
        function_queue.append(star)
    elif "bwa" in dict_keys:
        va = subprocess.check_output(
            ["cat", "../version_bwa-mem2_N.txt"], text=True
        ).strip()
        print(f"\t>>> Alignment: BWA-MEM2 (v{va})")
        function_queue.append(bwa)
    else:
        print("\t>>> Alignment: none")

    # Pseudoalignment
    if "salmon" in dict_keys:
        vp = subprocess.check_output(
            ["cat", "../version_salmon_N.txt"], text=True
        ).strip()
        print(f"\t>>> Pseudoalignment: Salmon (v{vp})")
        function_queue.append(salmon)
    else:
        print("\t>>> Pseudoalignment: none")

    # Sorting and indexing
    vsi = subprocess.check_output(
        ["cat", "../version_samtools_N.txt"], text=True
    ).strip()
    print(f"\t>>> Sorting/Indexing: Samtools (v{vsi})")
    function_queue.append(samtools)

    # MarkDuplicates
    vgatk = subprocess.check_output(["cat", "../version_gatk_N.txt"], text=True).strip()

    vpic = subprocess.check_output(
        ["cat", "../version_picard_N.txt"], text=True
    ).strip()
    print(f"\t>>> MarkDuplicates: GATK ({vgatk}) & Picard (v{vpic})")
    function_queue.append(markduplicates)

    # Quantification
    if "featurecounts" in dict_keys:
        vq = subprocess.check_output(
            ["cat", "../version_featurecounts_N.txt"], text=True
        ).strip()
        print(f"\t>>> Quantification: featureCounts ({vq})")
        function_queue.append(featurecounts)
    else:
        print("\t>>> Quantification: none")

    # MultiQC
    vmqc = subprocess.check_output(
        ["cat", "../version_multiqc_N.txt"], text=True
    ).strip()
    print(f"\t>>> Quality control report: multiQC (v{vmqc})")
    function_queue.append(multiqc)

    # Calling each steps
    for func in function_queue:
        func(sample, toml_config)

    end = get_time()
    total_time = end - start
    end_str = ">>> {}-seq pipeline for {} completed in {}.".format(
        toml_config["general"]["sequencing"], sample, total_time
    )
    print("=" * len(end_str) + "\n" + end_str + "\n" + "=" * len(end_str))


def get_time():
    now = datetime.now()
    return now


def title(message):
    print(f"\n>>> Running {message} [{get_time()}]")


def get_file_align(toml_config, output, sample):
    files = {}
    if "bbduk" in toml_config.keys():
        files = {
            "I1align": toml_config["general"]["output"]
            + "/"
            + sample
            + "/BBDuk/"
            + sample
            + "_trimmed_R1.fastq.gz",
            "I2align": toml_config["general"]["output"]
            + "/"
            + sample
            + "/BBDuk/"
            + sample
            + "_trimmed_R2.fastq.gz",
            "Ialign": toml_config["general"]["output"]
            + "/"
            + sample
            + "/BBDuk/"
            + sample
            + "_trimmed.fastq.gz",
            "Oalign": output + "/" + sample + "_trimmed",
        }

    else:
        files = {
            "I1align": toml_config["general"]["fastq"] + "/" + sample + "_R1.fastq.gz",
            "I2align": toml_config["general"]["fastq"] + "/" + sample + "_R2.fastq.gz",
            "Ialign": toml_config["general"]["fastq"] + "/" + sample + ".fastq.gz",
            "Oalign": output + "/" + sample,
        }

    return files


def get_reference(ref, tool):
    index_path = "/lustre03/project/6019267/shared/tools/REFERENCE_GENOMES/"
    gtf_path = "/lustre03/project/6019267/shared/tools/ANNOTATED_GTF/"
    reference: {}
    match ref:
        case "hg19":
            reference = {
                "index": index_path + "human/hg19_GRCh37/index_" + tool,
                "gtf": gtf_path + "human/hg19_GRCh37/Homo_sapiens.GRCh37.87.gtf.gz",
            }
        case "hg38":
            reference = {
                "index": index_path + "human/hg38_GRCh38/index_" + tool,
                "gtf": gtf_path + "human/hg38_GRCh38/Homo_sapiens.GRCh38.110.gtf.gz",
            }

        case "mm39":
            reference = {
                "index": index_path + "mouse/mm39_GRCm39/index_" + tool,
                "gtf": gtf_path + "mouse/mm39_GRCm39/Mus_musculus.GRCm39.110.gtf.gz",
            }

        case "ce11":
            reference = {
                "index": index_path + "worm/ce11_WBcel235/index_" + tool,
                "gtf": gtf_path
                + "worm/ce11_WBcel235/Caenorhabditis_elegans.WBcel235.110.gtf.gz",
            }

        case "danRer11":
            reference = {
                "index": index_path + "zebrafish/danRer11_GRCz11/index_" + tool,
                "gtf": gtf_path
                + "zebrafish/danRer11_GRCz11/Danio_rerio.GRCz11.110.gtf.gz",
            }
    return reference


def fastqc(sample, toml_config):
    title("FastQC")

    output = toml_config["general"]["output"] + "/" + sample + "/FastQC"
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
    subprocess.run(command)


def bbduk(sample, toml_config):
    title("BBDuk")

    output = toml_config["general"]["output"] + "/" + sample + "/BBDuk"
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
            "k=" + str(toml_config["bbduk"]["kmers"]),
            "qtrim=" + toml_config["bbduk"]["qtrim"],
            "trimq=" + str(toml_config["bbduk"]["trimq"]),
            "minlength=" + str(toml_config["bbduk"]["minlength"]),
            "mlf=" + str(toml_config["bbduk"]["mlf"]),
            "minavgquality=" + str(toml_config["bbduk"]["minavgquality"]),
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
            "k=" + str(toml_config["bbduk"]["kmers"]),
            "qtrim=" + toml_config["bbduk"]["qtrim"],
            "trimq=" + str(toml_config["bbduk"]["trimq"]),
            "minlength=" + str(toml_config["bbduk"]["minlength"]),
            "mlf=" + str(toml_config["bbduk"]["mlf"]),
            "minavgquality=" + str(toml_config["bbduk"]["minavgquality"]),
        ]
    command_str = " ".join(command)
    print(f">>> {command_str}\n")
    subprocess.run(command)


def star(sample, toml_config):
    title("STAR")
    output = toml_config["general"]["output"] + "/" + sample + "/STAR"
    subprocess.run(["mkdir", "-p", output])

    temporary = toml_config["general"]["temporary"] + "/" + sample + "/star_tmp"
    subprocess.run(["rm", "-r", temporary])

    ref = get_reference(toml_config["general"]["reference"], "star")["index"]

    files = get_file_align(toml_config, output, sample)
    I1align = files["I1align"]
    I2align = files["I2align"]
    Ialign = files["Ialign"]
    Oalign = files["Oalign"]

    if toml_config["general"]["reads"] == "PE":
        command = [
            "STAR",
            "--runMode",
            "alignReads",
            "--runThreadN",
            str(toml_config["general"]["threads"]),
            "--genomeDir",
            ref,
            "--outFileNamePrefix",
            Oalign,
            "--outTmpDir",
            temporary,
            "--readFilesCommand",
            "zcat",
            "--readFilesIn",
            I1align,
            I2align,
            "--outSAMtype",
            toml_config["star"]["outSAMtype1"],
            toml_config["star"]["outSAMtype2"],
            "--twopassMode",
            toml_config["star"]["twopassMode"],
            "--outWigType",
            toml_config["star"]["outWigType"],
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
            "--genomeDir",
            ref,
            "--outFileNamePrefix",
            Oalign,
            "--outTmpDir",
            temporary,
            "--readFilesCommand",
            "zcat",
            "--readFilesIn",
            Ialign,
            "--outSAMtype",
            toml_config["star"]["outSAMtype1"],
            toml_config["star"]["outSAMtype2"],
            "--twopassMode",
            toml_config["star"]["twopassMode"],
            "--outWigType",
            toml_config["star"]["outWigType"],
            "--outSJtype",
            toml_config["star"]["outSJtype"],
            "--quantMode",
            toml_config["star"]["quantMode"],
        ]
    command_str = " ".join(command)
    print(f">>> {command_str}\n")
    subprocess.run(command)


def bwa(sample, toml_config):
    title("BWA-MEM2")

    output = toml_config["general"]["output"] + "/" + sample + "/BWA-MEM2"
    subprocess.run(["mkdir", "-p", output])

    ref = get_reference(toml_config["general"]["reference"], "star")["index"]

    files = get_file_align(toml_config, output, sample)
    I1align = files["I1align"]
    I2align = files["I2align"]
    Ialign = files["Ialign"]
    Oalign = files["Oalign"]

    if toml_config["general"]["reads"] == "PE":
        command = []
    else:
        command = []
    command_str = " ".join(command)
    print(f">>> {command_str}\n")
    subprocess.run(command)


def salmon(sample, toml_config):
    title("Salmon")


def samtools(sample, toml_config):
    title("Samtools")


def markduplicates(sample, toml_config):
    title("MarkDuplicates")


def featurecounts(sample, toml_config):
    title("FeatureCounts")


def multiqc(sample, toml_config):
    title("MutliQC")


if __name__ == "__main__":
    main()

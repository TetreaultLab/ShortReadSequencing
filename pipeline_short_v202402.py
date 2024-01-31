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
    function_queue = []
    print(">>> Parameters:")

    # # Quality control
    # print("\t>>> Quality control: FastQC (v0.12.1)")
    # function_queue.append(fastqc)

    # # Trimming
    # if toml_config["general"]["trimming"] == "bbduk":
    #     print("\t>>> Trimming: BBDuk (v39.06)")
    #     function_queue.append(bbduk)
    # else:
    #     print("\t>>> Trimming: none")

    # # Alignment
    # if toml_config["general"]["alignment"] == "star":
    #     print("\t>>> Alignment: STAR (v2.7.11a)")
    #     function_queue.append(star)
    # elif toml_config["general"]["alignment"] == "bwa":
    #     print("\t>>> Alignment: BWA-MEM2 (v2.2.1)")
    #     function_queue.append(bwa)
    # else:
    #     print("\t>>> Alignment: none")

    # # Pseudo alignedment
    # if toml_config["general"]["pseudo"] == "salmon":
    #     print("\t>>> Pseudo_alignedment: Salmon (v1.10.2)")
    #     function_queue.append(salmon)
    # else:
    #     print("\t>>> Pseudo_alignedment: none")

    # # Sorting and indexing
    # print("\t>>> Sorting/Indexing: Samtools (v1.18)")
    # function_queue.append(samtools)

    # # MarkDuplicates
    # print("\t>>> MarkDuplicates: GATK (4.4.0.0) & Picard (v3.0.0)")
    # function_queue.append(markduplicates)

    # # Quantification
    # if toml_config["general"]["quantification"] == "featurecounts":
    #     print("\t>>> Quantification: featureCounts (v2.0.6)")
    #     function_queue.append(featurecounts)
    # else:
    #     print("\t>>> Quantification: none")

    # # MultiQC
    # print("\t>>> Quality control report: MultiQC (v1.18)")
    # function_queue.append(multiqc)

    # # Variant Calling
    # # BAM to VCF
    print("\t>>> Variant Calling: BCFtools (v1.18)")
    function_queue.append(bcftools)

    # # Variant Effect Predictor
    # if toml_config["general"]["variant"] == "vep":
    #     print("\t>>> Variant Calling: Variant Effect Predictor (VEP) (v2.0.6)")
    #     function_queue.append(vep)
    # else:
    #     print("\t>>> Variant Calling: none")

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
    index_path = "/lustre03/project/6019267/shared/tools/REFERENCE_GENOMES/"
    gtf_path = "/lustre03/project/6019267/shared/tools/ANNOTATED_GTF/"
    reference: {}
    match ref:
        case "hg19":
            reference = {
                "fasta": index_path
                + "human/hg19_GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa",
                "index": index_path + "human/hg19_GRCh37/index_" + tool,
                "gtf": gtf_path + "human/hg19_GRCh37/Homo_sapiens.GRCh37.87.gtf.gz",
                "gff3": gtf_path + "human/hg19_GRCh37/Homo_sapiens.GRCh37.87.gff3.gz",
            }
        case "hg38":
            reference = {
                "fasta": index_path
                + "human/hg38_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
                "index": index_path + "human/hg38_GRCh38/index_" + tool,
                "gtf": gtf_path + "human/hg38_GRCh38/Homo_sapiens.GRCh38.110.gtf.gz",
                "gff3": gtf_path + "human/hg38_GRCh38/Homo_sapiens.GRCh38.110.gff3.gz",
            }

        case "mm39":
            reference = {
                "fasta": index_path
                + "mouse/mm39_GRCm39/Mus_musculus.GRCm39.dna.primary_assembly.fa",
                "index": index_path + "mouse/mm39_GRCm39/index_" + tool,
                "gtf": gtf_path + "mouse/mm39_GRCm39/Mus_musculus.GRCm39.110.gtf.gz",
                "gff3": gtf_path + "mouse/mm39_GRCm39/Mus_musculus.GRCm39.110.gff3.gz",
            }

        case "ce11":
            reference = {
                "fasta": index_path
                + "worm/ce11_WBcel235/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa",
                "index": index_path + "worm/ce11_WBcel235/index_" + tool,
                "gtf": gtf_path
                + "worm/ce11_WBcel235/Caenorhabditis_elegans.WBcel235.110.gtf.gz",
                "gff3": gtf_path
                + "worm/ce11_WBcel235/Caenorhabditis_elegans.WBcel235.110.gff3.gz",
            }

        case "danRer11":
            reference = {
                "fasta": index_path
                + "zebrafish/danRer11_GRCz11/Danio_rerio.GRCz11.dna.primary_assembly.fa",
                "index": index_path + "zebrafish/danRer11_GRCz11/index_" + tool,
                "gtf": gtf_path
                + "zebrafish/danRer11_GRCz11/Danio_rerio.GRCz11.110.gtf.gz",
                "gff3": gtf_path
                + "zebrafish/danRer11_GRCz11/Danio_rerio.GRCz11.110.gff3.gz",
            }
    return reference


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
            "2",
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
            "--genomeDir",
            ref,
            "--outFileNamePrefix",
            O_aligned,
            "--outTmpDir",
            temporary,
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
            "--genomeDir",
            ref,
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
            "--outWigType",
            toml_config["star"]["outWigType"],
            "--outWigStrand",
            toml_config["star"]["outWigStrand"],
            "--outSJtype",
            toml_config["star"]["outSJtype"],
            "--quantMode",
            toml_config["star"]["quantMode"],
        ]
    command_str = " ".join(command)
    print(f">>> {command_str}\n")
    subprocess.run(command)

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


def bwa(sample, toml_config):
    title("BWA-MEM2")

    output = toml_config["general"]["output"] + "/" + sample + "/Aligned"
    subprocess.run(["mkdir", "-p", output])

    ref = (
        get_reference(toml_config["general"]["reference"], "bwa-mem2")["index"]
        + "/bwa-mem2"
    )

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
            ref,
            I_toAlign,
        ]

    command_str = " ".join(command)
    print(f">>> {command_str}\n")
    subprocess.run(command)

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
        ]
    )


def salmon(sample, toml_config):
    title("Salmon")

    output = toml_config["general"]["output"] + "/" + sample + "/Salmon"
    subprocess.run(["mkdir", "-p", output])

    ref = get_reference(toml_config["general"]["reference"], "salmon")["index"]
    gtf = get_reference(toml_config["general"]["reference"], "salmon")["gtf"]

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
            "--index",
            ref,
            "--threads",
            str(toml_config["general"]["threads"]),
            "--auxDir",
            "salmon_tmp",
            "--geneMap",
            gtf,
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
            "--geneMap",
            gtf,
            "--minScoreFraction",
            str(toml_config["salmon"]["minScoreFraction"]),
            "--unmatedReads",
            I_toAlign,
            "--output",
            output,
        ]

    command_str = " ".join(command)
    print(f">>> {command_str}\n")
    subprocess.run(command)


def samtools(sample, toml_config):
    title("Samtools")
    input = toml_config["general"]["output"] + "/" + sample + "/Aligned"
    output = toml_config["general"]["output"] + "/" + sample + "/Samtools"
    subprocess.run(["mkdir", "-p", output])

    inBAM = input + "/" + sample + ".bam"
    bamName = output + "/" + sample + "_sortedName.bam"
    bamCoord = output + "/" + sample + "_sortedCoordinate.bam"
    idxstats = output + "/" + sample + "_sortedCoordinate_idxstats.txt"

    # Sort by name
    subprocess.run(["samtools", "sort", "-n", inBAM, "-o", bamName])

    # Sort by coordinate
    subprocess.run(["samtools", "sort", inBAM, "-o", bamCoord])

    # Index bam sorted by coordinates
    subprocess.run(["samtools", "index", "-b", bamCoord, "-o", bamCoord + ".bai"])

    # idxstats
    subprocess.run(["samtools", "idxstats", bamCoord, ">", idxstats])


def markduplicates(sample, toml_config):
    title("MarkDuplicates")

    output = toml_config["general"]["output"] + "/" + sample + "/MarkDuplicates/"
    subprocess.run(["mkdir", "-p", output])

    temporary = toml_config["general"]["temporary"] + "/" + sample + "/md_tmp"
    subprocess.run(["mkdir", "-p", temporary])

    input = toml_config["general"]["output"] + "/" + sample + "/Samtools"
    bamCoord = input + "/" + sample + "_sortedCoordinate.bam"
    metrics = output + sample + "_duplicates_metrics.txt"
    records = output + sample + "_duplicates.bam"

    command = [
        "gatk",
        "MarkDuplicates",
        "--INPUT",
        bamCoord,
        "--METRICS_FILE",
        metrics,
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
    subprocess.run(command)


def featurecounts(sample, toml_config):
    title("FeatureCounts")

    temporary = toml_config["general"]["temporary"] + "/" + sample
    output = toml_config["general"]["output"] + "/" + sample + "/FeatureCounts"
    subprocess.run(["mkdir", "-p", output])

    input = (
        toml_config["general"]["output"]
        + "/"
        + sample
        + "/Samtools/"
        + sample
        + "_sortedCoordinate.bam"
    )

    gtf = get_reference(toml_config["general"]["reference"], "")["gtf"]
    print(gtf)

    if toml_config["general"]["reads"] == "PE":
        command = [
            "featureCounts",
            "-t",
            toml_config["featurecounts"]["features"],
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
            output + "/" + sample + "_counts.txt",
            input,
        ]
    else:
        command = [
            "featureCounts",
            "-t",
            toml_config["featurecounts"]["features"],
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
            output + "/" + sample,
            input,
        ]

    command_str = " ".join(command)
    print(f">>> {command_str}\n")
    subprocess.run(command)


def multiqc(sample, toml_config):
    title("MutliQC")
    input = toml_config["general"]["output"] + "/" + sample + "/"
    output = toml_config["general"]["output"] + "/" + sample + "/QC/multiQC/"
    subprocess.run(["mkdir", "-p", output])

    with open(output + "/my_file_list.txt", "w") as f:
        f.write(input + "QC/fastQC/\n")
        if toml_config["general"]["trimming"] != "none":
            f.write(input + "Trimmed/\n")
        if toml_config["general"]["alignment"] != "none":
            f.write(input + "Aligned/\n")
        if toml_config["general"]["pseudo"] != "none":
            f.write(input + "Salmon/\n")
        if toml_config["general"]["quantification"] != "none":
            f.write(input + "FeatureCounts/\n")
        f.write(input + "Samtools/\n")
        f.write(input + "MarkDuplicates/\n")
    f.close()

    subprocess.run(
        [
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


def bcftools(sample, toml_config):
    title("BCFtools")
    input = (
        toml_config["general"]["output"]
        + "/"
        + sample
        + "/Samtools/"
        + sample
        + "_sortedCoordinate.bam"
    )
    output = toml_config["general"]["output"] + "/" + sample + "/BCFtools/"
    subprocess.run(["mkdir", "-p", output])

    ref = get_reference(toml_config["general"]["reference"], "")["fasta"]
    gff3 = get_reference(toml_config["general"]["reference"], "")["gff3"]

    mpileup = [
        "bcftools",
        "mpileup",
        "--threads",
        str(toml_config["general"]["threads"]),
        "-d",
        "10000",
        "-Ob",
        "-o",
        output + sample + ".bcf.gz",
        "-f",
        ref,
        input,
    ]

    call = [
        "bcftools",
        "call",
        "-m",
        "-Ob",
        "-o",
        output + sample + "_calls.bcf.gz",
        output + sample + ".bcf.gz",
    ]

    consequence = [
        "bcftools",
        "csq",
        "-f",
        ref,
        "-g",
        gff3,
        output + sample + "_calls.bcf.gz",
        "-Oz",
        "-o",
        output + sample + "_var.vcf",
    ]

    command_1 = " ".join(mpileup)
    print(f">>> {command_1}\n")
    subprocess.run(mpileup)

    command_2 = " ".join(call)
    print(f">>> {command_2}\n")
    subprocess.run(call)

    command_3 = " ".join(consequence)
    print(f">>> {command_3}\n")
    subprocess.run(consequence)


def vep(sample, toml_config):
    title("VEP")

    vep_cache = "/lustre03/project/6019267/shared/tools/VARIATION_VEP"
    ref = get_reference(toml_config["general"]["reference"], "")["fasta"]
    gff3 = get_reference(toml_config["general"]["reference"], "")["gff3"]

    input = (
        toml_config["general"]["output"]
        + "/"
        + sample
        + "/BCFtools/"
        + sample
        + "_var.vcf"
    )
    output = (
        toml_config["general"]["output"]
        + "/"
        + sample
        + "/Variants/"
        + sample
        + "_all_variants.txt"
    )

    command = [
        "vep",
        "--cache",
        "--dir_cache",
        vep_cache,
        "--species",
        toml_config["vep"]["species"],
        "--gff",
        gff3,
        "--fasta",
        ref,
        "--cache_version",
        "110",
        "--everything",
        "--tab",
        "-i",
        input,
        "-o",
        output,
    ]

    command_str = " ".join(command)
    print(f">>> {command_str}\n")
    subprocess.run(command)


if __name__ == "__main__":
    main()

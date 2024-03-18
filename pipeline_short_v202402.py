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

    # Open file for steps done
    steps = open(output + "/steps_done.txt", "a")
    steps.write("\nLoading ENV\n")
    steps.close()

    done = []
    with open(output + "/steps_done.txt", "r") as f:
        for line in f:
            done.append(line.strip())

    print(done)

    # Get tools and versions
    function_queue = []
    print(">>> Parameters:")

    # Quality control
    print("\t>>> Quality control: FastQC (v0.12.1)")
    if "FastQC" not in done:
        function_queue.append(fastqc)

    # Trimming
    if toml_config["general"]["trimming"] == "bbduk":
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
    if toml_config["general"]["pseudo"] == "salmon":
        print("\t>>> Pseudo_alignedment: Salmon (v1.10.2)")
        if "Salmon" not in done:
            function_queue.append(salmon)
    else:
        print("\t>>> Pseudo_alignedment: none")

    # Sorting and indexing
    print("\t>>> Sorting/Indexing: Samtools (v1.18)")
    if "Samtools" not in done:
        function_queue.append(samtools)

    # Alignment QC
    if "FastQC for bam" not in done:
        function_queue.append(bamqc)

    # MarkDuplicates
    print("\t>>> MarkDuplicates: GATK (4.4.0.0) & Picard (v3.0.0)")
    if "MarkDuplicates" not in done:
        function_queue.append(markduplicates)

    # Quantification
    if toml_config["general"]["quantification"] == "featurecounts":
        print("\t>>> Quantification: featureCounts (v2.0.6)")
        if "FeatureCounts" not in done:
            function_queue.append(featurecounts)
    else:
        print("\t>>> Quantification: none")

    # MultiQC
    print("\t>>> Quality control report: MultiQC (v1.18)")
    if "MutliQC" not in done:
        function_queue.append(multiqc)

    # Variant Calling : BAM to VCF
    print("\t>>> Variant Calling: BCFtools (v1.18)")
    if "BCFtools" not in done:
        function_queue.append(bcftools)

    # Variant Calling : SnpEff (annotation)
    if toml_config["general"]["variant"] == "snpeff":
        print("\t>>> Variant Calling: SnpEff + SnpSift (v5.2a)")
        if "SnpEff" not in done:
            function_queue.append(snpeff)
    else:
        print("\t>>> Variant Calling: none")

    # Calling each steps
    for func in function_queue:
        func(sample, toml_config)

    steps.close()

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
    path = "/lustre03/project/6019267/shared/tools/PIPELINES/References/"
    reference: {}  # type: ignore
    match ref:
        case "grch37":
            reference = {
                "fasta": path + "Homo_sapiens.GRCh37.87.dna.primary_assembly.fa",
                "index": path + "index_" + tool + "/" + ref,
                "gtf": path + "Homo_sapiens.GRCh37.87.gtf",
                "gff3": path + "Homo_sapiens.GRCh37.87.gff3.gz",
            }
        case "grch38":
            reference = {
                "fasta": path + "Homo_sapiens.GRCh38.105.dna.primary_assembly.fa",
                "index": path + "index_" + tool + "/" + ref,
                "gtf": path + "Homo_sapiens.GRCh38.105.gtf",
                "gff3": path + "Homo_sapiens.GRCh38.105.gff3.gz",
            }

        case "grcm39":
            reference = {
                "fasta": path + "Mus_musculus.GRCm39.105.dna.primary_assembly.fa",
                "index": path + "index_" + tool + "/" + ref,
                "gtf": path + "Mus_musculus.GRCm39.105.gtf",
                "gff3": path + "Mus_musculus.GRCm39.105.gff3.gz",
            }

        case "wbcel235":
            reference = {
                "fasta": path + "Caenorhabditis_elegans.WBcel235.105.dna.toplevel.fa",
                "index": path + "index_" + tool + "/" + ref,
                "gtf": path + "Caenorhabditis_elegans.WBcel235.105.gtf",
                "gff3": path + "Caenorhabditis_elegans.WBcel235.105.gff3.gz",
            }

        case "grcz11":
            reference = {
                "fasta": path + "Danio_rerio.GRCz11.105.dna.primary_assembly.fa",
                "index": path + "index_" + tool + "/" + ref,
                "gtf": path + "Danio_rerio.GRCz11.105.gtf",
                "gff3": path + "Danio_rerio.GRCz11.105.gff3.gz",
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

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("FastQC\n")


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
            "34105090798",
            "--genomeDir",
            ref,
            "--outFileNamePrefix",
            O_aligned,
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

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("STAR\n")


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
    subprocess.run(command)

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
    subprocess.run(command)

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("FastQC for bam\n")


def markduplicates(sample, toml_config):
    title("MarkDuplicates")

    output = toml_config["general"]["output"] + "/" + sample + "/MarkDuplicates/"
    subprocess.run(["mkdir", "-p", output])

    temporary = toml_config["general"]["temporary"] + "/" + sample + "/md_tmp"
    subprocess.run(["mkdir", "-p", temporary])

    input = toml_config["general"]["output"] + "/" + sample + "/Aligned"
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
    print(gtf)

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
    subprocess.run(command)

    # Work on counts output
    with open(output + "/" + sample + "_0.counts", "w") as outfile:
        subprocess.run(
            ["tail", "-n", "+2", output + "/" + sample + "_geneID.txt"], stdout=outfile
        )

    with open(output + "/" + sample + "_1.counts", "w") as outfile:
        subprocess.run(
            ["cut", "-f1,7", output + "/" + sample + "_0.counts"], stdout=outfile
        )

    subprocess.run(
        [
            "sed",
            "-i",
            "1s|.*|gene_id\t" + sample + "|",
            output + "/" + sample + "_1.counts",
        ]
    )

    with open(output + "/" + sample + "_2.counts", "w") as outfile:
        subprocess.run(["sort", output + "/" + sample + "_1.counts"], stdout=outfile)

    with open(output + "/" + sample + "_3.counts", "w") as outfile:
        subprocess.run(
            [
                "join",
                output + "/" + sample + "_2.counts",
                "/lustre04/scratch/mlab/pipeline2024/ShortReadSequencing/"
                + toml_config["general"]["reference"]
                + "_gene_association.txt",
            ],
            stdout=outfile,
        )
    # bring last line to the top
    with open(output + "/" + sample + "_4.counts", "w") as outfile:
        subprocess.run(
            ["sed", "$!H;1h;$!d;G", output + "/" + sample + "_3.counts"], stdout=outfile
        )

    with open(output + "/" + sample + "_biotype_counts.txt", "w") as outfile:
        subprocess.run(
            [
                "awk",
                'BEGIN {FS=" "; OFS="\t"} {print $3, $2, $1, $4}',
                output + "/" + sample + "_4.counts",
            ],
            stdout=outfile,
        )

    with open(output + "/" + sample + ".counts", "w") as outfile:
        subprocess.run(
            [
                "awk",
                'BEGIN {FS=" "; OFS="\t"} {print $3, $2}',
                output + "/" + sample + "_4.counts",
            ],
            stdout=outfile,
        )

    subprocess.run(["rm", output + "/" + sample + "_0.counts"])
    subprocess.run(["rm", output + "/" + sample + "_1.counts"])
    subprocess.run(["rm", output + "/" + sample + "_2.counts"])
    subprocess.run(["rm", output + "/" + sample + "_3.counts"])
    subprocess.run(["rm", output + "/" + sample + "_4.counts"])

    # Steps done
    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("FeatureCounts\n")


def multiqc(sample, toml_config):
    title("MutliQC")
    input = toml_config["general"]["output"] + "/" + sample + "/"
    output = toml_config["general"]["output"] + "/" + sample + "/QC/multiQC/"
    subprocess.run(["mkdir", "-p", output])

    with open(output + "/my_file_list.txt", "w") as f:
        f.write(input + "QC/fastQC/\n")
        if toml_config["general"]["trimming"] != "none":
            f.write(input + "Trimmed/\n")
        if toml_config["general"]["pseudo"] != "none":
            f.write(input + "Salmon/\n")
        if toml_config["general"]["quantification"] != "none":
            f.write(input + "FeatureCounts/\n")
        f.write(input + "Aligned/\n")
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

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("MutliQC\n")


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

    mpileup = [
        "bcftools",
        "mpileup",
        "--threads",
        str(toml_config["general"]["threads"]),
        "-d",
        "10000",
        "-o",
        output + sample + ".vcf",
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
        output + sample + "_all.vcf",
        output + sample + ".vcf",
    ]

    command_1 = " ".join(mpileup)
    print(f">>> {command_1}\n")
    subprocess.run(mpileup)

    command_2 = " ".join(call)
    print(f">>> {command_2}\n")
    subprocess.run(call)

    with open(
        toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    ) as steps:
        steps.write("BCFtools\n")


def snpeff(sample, toml_config):
    title("SnpEff")

    snpeff = (
        "/lustre03/project/6019267/shared/tools/PIPELINES/ShortReadSequencing/snpEff"
    )

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
        ref == "WBcel235.105"

    path = toml_config["general"]["output"] + "/" + sample + "/Variants"

    cmd_snpeff = [
        "java",
        "-Xmx8g",
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
        path + "/" + sample + "_all.vcf",
    ]

    command_str1 = " ".join(cmd_snpeff)
    print(f">>> {command_str1}\n")

    cmd_vartype = [
        "java",
        "-Xmx8g",
        "-jar",
        snpeff + "/SnpSift.jar",
        "varType",
        "-noLog",
        path + "/" + sample + "_snpeff.vcf",
    ]

    command_str2 = " ".join(cmd_vartype)
    print(f">>> {command_str2}\n")

    cmd_extract = [
        "java",
        "-Xmx8g",
        "-jar",
        snpeff + "/SnpSift.jar",
        "extractFields",
        "-noLog",
        "-s",
        "|",
        "-e",
        ".",
        path + "/" + sample + "_vartype.vcf",
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "QUAL",
        "MQBZ",
        "RPBZ",
        "SCBZ",
        "DP4",
        "HOM",
        "AN",
        "AC",
        "VARTYPE",
        "ANN[*].GENE",
        "ANN[*].FEATUREID",
        "ANN[*].EFFECT",
        "ANN[*].IMPACT",
        "ANN[*].BIOTYPE",
        "ANN[*].HGVS_C",
        "ANN[*].HGVS_P",
    ]

    command_str3 = " ".join(cmd_extract)
    print(f">>> {command_str3}\n")

    with open(path + "/" + sample + "_snpeff.vcf", "w") as outfile:
        subprocess.run(
            cmd_snpeff,
            stdout=outfile,
        )

    with open(path + "/" + sample + "_vartype.vcf", "w") as outfile:
        subprocess.run(
            cmd_vartype,
            stdout=outfile,
        )

    with open(path + "/" + sample + "_all.txt", "w") as outfile:
        subprocess.run(
            cmd_extract,
            stdout=outfile,
        )

    # with open(
    #     toml_config["general"]["output"] + "/" + sample + "/steps_done.txt", "a"
    # ) as steps:
    #     steps.write("SnpEff\n")


if __name__ == "__main__":
    main()

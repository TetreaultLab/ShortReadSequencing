import argparse
import toml
from datetime import datetime
import subprocess


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
    start_str = "### {}-seq pipeline starting for {} at {} ###".format(
        toml_config["general"]["sequencing"], sample, start
    )
    print("#" * len(start_str) + "\n" + start_str + "\n" + "#" * len(start_str))

    # Get tools and versions
    dict_keys = toml_config.keys()

    function_queue = []
    # Quality control
    vqc = subprocess.check_output(
        ["cat ../version_fastqc_N.txt"], text=True, shell=True
    ).strip()
    print(f"Quality control: FastQC ({vqc})")

    # Trimming
    if "bbduk" in dict_keys:
        vt = subprocess.check_output(
            ["cat ../version_bbduk_N.txt"], text=True, shell=True
        ).strip()
        print(f"Trimming: BBDuk (v{vt})")
        function_queue.append(bbduk)
    else:
        print("Trimming: none")

    # Alignment
    if "star" in dict_keys:
        va = subprocess.check_output(["STAR --version"], text=True, shell=True).strip()
        print(f"Alignment: STAR (v{va})")
    elif "bwa" in dict_keys:
        va = subprocess.check_output(
            ["bwa-mem2 version"], text=True, shell=True
        ).strip()
        print(f"Alignment: BWA-MEM2 (v{va})")
    else:
        print("Alignment: none")

    # Pseudoalignment
    if "salmon" in dict_keys:
        vp = (
            subprocess.check_output(["salmon --version"], text=True, shell=True)
            .strip()
            .split(" ")[1]
        )
        print(f"Pseudoalignment: Salmon (v{vp})")
    else:
        print("Pseudoalignment: none")

    # Quantification
    if "featurecounts" in dict_keys:
        vq = subprocess.check_output(
            ["cat ../version_featurecounts_N.txt"], text=True, shell=True
        ).strip()
        print(f"Quantification: featureCounts ({vq})")
    else:
        print("Quantification: none")

    # Sorting and indexing
    vsi = subprocess.check_output(
        ["samtools version | head -n 1 | awk '{print $2}'"], text=True, shell=True
    ).strip()
    print(f"Sorting/Indexing: Samtools (v{vsi})")

    # Calling initial quality control
    fastqc()

    # Calling each steps
    for func in function_queue:
        func(toml_config)
    # Steps
    # 01. Quality control
    # 02. Trimming
    # 03. Quality control
    # 04. Alignment
    # 05. Quality control
    # 06. Pseudoalignment
    # 07. Quality control
    # 08. Sorting/indexing
    # 09. If PE: markduplicates
    # 10. Quantification

    # reads = toml_config['general']['reads']
    # fastq = toml_config['general']['fastq']

    # if reads == "PE":
    #     R1 = fastq + "/" + "" + sample + "_R1.fastq.gz"
    #     R2 = fastq + "/" + "" + sample + "_R2.fastq.gz"
    # elif reads == "SE":
    #     R = fastq + "/" + "" + sample + ".fastq.gz"

    end = get_time()
    total_time = end - start
    end_str = "### {}-seq pipeline for {} completed in {} ###".format(
        toml_config["general"]["sequencing"], sample, total_time
    )
    print("#" * len(end_str) + "\n" + end_str + "\n" + "#" * len(end_str))


def get_time():
    now = datetime.now()
    return now


def title(message):
    print(f"\n>>> Running {message} [{get_time()}]\n")


def get_file_names():
    print()


def fastqc():
    title("FastQC")


def bbduk(toml_config):
    title("BBDuk")
    print(toml_config["bbduk"])
    # subprocess.call(["bbduk.sh", "--version"], text=True)


if __name__ == "__main__":
    main()

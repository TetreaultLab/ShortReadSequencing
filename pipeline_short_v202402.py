import argparse
import toml
from datetime import datetime

def main():
    start = get_time()

    parser = argparse.ArgumentParser(
                    prog='PipelineShort',
                    description='ShortReadSequencing pipeline for RNA, Exome and Genome data.')

    parser.add_argument('--sample', type=str, required=True, help='Sample name.')
    parser.add_argument('--config', type=str, required=True, help='Project config file, including path.')

    args = parser.parse_args()
    
    # Sample name
    sample = args.sample

    # Loading TOML config
    with open(args.config, 'r') as f:
        toml_config = toml.load(f)

    # Start of pipeline
    start_str = "### {}-seq pipeline starting for {} at {} ###".format(toml_config['general']['sequencing'], sample, start.strftime("%d/%m/%Y %H:%M:%S"))
    print("#"*len(start_str) + "\n" + start_str + "\n" + "#"*len(start_str))
    
    # Get tools
    dict_keys = toml_config.keys()
    print(dict_keys)

    print()

    # reads = toml_config['general']['reads']
    # fastq = toml_config['general']['fastq']

    # if reads == "PE":
    #     R1 = fastq + "/" + "" + sample + "_R1.fastq.gz"
    #     R2 = fastq + "/" + "" + sample + "_R2.fastq.gz"
    # elif reads == "SE":
    #     R = fastq + "/" + "" + sample + ".fastq.gz"
    
    end = get_time()
    total_time = end - start
    end_str = "### {}-seq pipeline for {} completed in {} ###".format(toml_config['general']['sequencing'], sample, total_time)
    print("#"*len(end_str) + "\n" + end_str + "\n" + "#"*len(end_str))


def get_time():
    now = datetime.now()
    return now

def get_file_names():
    print()


if __name__ == "__main__":
    main()
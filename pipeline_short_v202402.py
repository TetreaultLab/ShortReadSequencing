import argparse
import toml
import os

# Loading TOML config
parser = argparse.ArgumentParser(
                    prog='PipelineShort',
                    description='ShortReadSequencing pipeline for RNA, Exome and Genome data.')

parser.add_argument('--config', type=str, required=True, help='Project config file, including path.')

args = parser.parse_args()

with open(args.config, 'r') as f:
    toml_config = toml.load(f)

print(toml_config)
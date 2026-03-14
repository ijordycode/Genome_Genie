import os
import argparse
import pandas as pd
import yaml
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument("--ext", type=str, required=True)
parser.add_argument("--read1", type=str, required=True)
parser.add_argument("--read2", type=str, required=True)

## args for which pipeline to run
parser.add_argument("--vc", default=False, action='store_true', required=False)
parser.add_argument("--gq", default=False, action='store_true', required=False)
parser.add_argument("--deseq", default=False, action='store_true', required=False)
args = parser.parse_args()

cwd_directory = os.getcwd()

## Read config yaml file
yaml_file = "config.yaml"
with open(yaml_file, "r") as f:
    config = yaml.safe_load(f)
    config['fastq_ext'] = args.ext
    config['read1'] = args.read1
    config['read2'] = args.read2
    config['run_variant_calling'] = args.vc
    config['run_gene_qunatification'] = args.gq
    config['run_deseq2'] = args.deseq

# 3. Write the updated data back to the file
with open(yaml_file, 'w') as file:
    # Use safe_dump for safer output and specify sort_keys=False to preserve order
    yaml.safe_dump(config, file)



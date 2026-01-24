import os
import argparse
import pandas as pd
import yaml
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument("--ext", type=str, required=False)
parser.add_argument("--read1", type=str, required=False)
parser.add_argument("--read2", type=str, required=False)
args = parser.parse_args()

cwd_directory = os.getcwd()

## Read config yaml file
yaml_file = "config.yaml"
with open(yaml_file, "r") as f:
    config = yaml.safe_load(f)

## Add timestamp to yaml file
config['timestamp'] = datetime.now().strftime('output_%Y_%m_%d_%H_%M_%S')
with open('config.yaml', 'w') as f:
    yaml.dump(config, f)

## Make timestamp directory
full_output_path = os.path.join(cwd_directory, "outputs", config['timestamp'])
os.mkdir(full_output_path)

def generate_samplesheet():
    samples = []
    for file in os.listdir(os.path.join(cwd_directory, "fastq_files")):
        samplename = file.replace(config['fastq_ext'], "")
        for read in [config['read1'], config['read2']]:
            samplename = samplename.replace(read,"")

        samples.append(samplename) if samplename not in samples else None

    df = pd.DataFrame({"samples": samples, "condition":None})
    df.to_csv(os.path.join(full_output_path, "samplesheet.csv"), index=False)

generate_samplesheet()

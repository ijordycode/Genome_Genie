import os
import argparse
import pandas as pd
import yaml

parser = argparse.ArgumentParser()
parser.add_argument("--ext", type=str, required=False)
parser.add_argument("--read1", type=str, required=False)
parser.add_argument("--read2", type=str, required=False)
args = parser.parse_args()

cwd_directory = os.getcwd()

yaml_file = "config.yaml"
with open(yaml_file, "r") as f:
    config = yaml.safe_load(f)

def generate_samplesheet():
    samples = []
    for file in os.listdir(os.path.join(cwd_directory, "Fastq_files")):
        samplename = file.replace(config['fastq_ext'], "")
        for read in [config['read1'], config['read2']]:
            samplename = samplename.replace(read,"")

        samples.append(samplename) if samplename not in samples else None

    df = pd.DataFrame({"samples": samples, "condition":None})
    output_dir = os.path.join(cwd_directory, "Outputs/samplesheet.csv")
    df.to_csv(os.path.join(cwd_directory, "Outputs", "samplesheet.csv"), index=False)

generate_samplesheet()




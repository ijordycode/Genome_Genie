import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--ext", type=str, required=True)
parser.add_argument("--read1", type=str, required=True)
parser.add_argument("--read2", type=str, required=True)
args = parser.parse_args()

cwd_directory = os.getcwd()
extenstion = args.ext

read_ext = [args.read1, args.read2]

samples = []
for file in os.listdir(os.path.join(cwd_directory, "Fastq_files")):
    samplename = file.replace(args.ext, "")
    for read in read_ext:
        samplename = samplename.replace(read,"")

    samples.append(samplename) if samplename not in samples else None

df = pd.DataFrame({"samples": samples})
output_dir = os.path.join(cwd_directory, "Outputs/samplesheet.csv")
df.to_csv(os.path.join(cwd_directory, "Outputs", "samplesheet.csv"), index=False)

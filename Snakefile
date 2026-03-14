import os
import pandas as pd
from datetime import datetime


configfile: '/BIP/src/config.yaml'

configfile: 'settings.yaml'

## Directories 
reference_dirs = config['reference_dirs']
output_dir = os.path.join(config['output_dir'], config['timestamp'])
trimmed_fastq_files = os.path.join(output_dir, config['fastq_dir'])
sam_files = os.path.join(output_dir, config['sam_files_output'])
bam_files = os.path.join(output_dir, config['bam_files_output'])
vcf_files = os.path.join(output_dir, config['vcf_files_output'])
count_files = os.path.join(output_dir, config['count_files_output'])

samplesheet_path = os.path.join(output_dir, "samplesheet.csv")
samples_df = pd.read_csv(samplesheet_path)
sample_list = samples_df['samples'].tolist()

if config.get("run_variant_calling", True):
    include: "subworkflow/variant_calling.smk"

if config.get("run_gene_quantification", True):
    include: "subworkflow/gene_counts.smk"

rule all:
    input:
        expand(os.path.join(sam_files, "{sample}.sam"), sample=sample_list),
        expand( os.path.join(bam_files, "{sample}.bam.bai"), sample=sample_list),
        expand(os.path.join(vcf_files, "{sample}.vcf"), sample=sample_list) if config.get("run_variant_calling") else [],
        expand(os.path.join(count_files, "{sample}_counts.txt"), sample=sample_list) if config.get("run_gene_quantification") else [],
        expand(os.path.join(count_files, "deseq2_results.csv")) if config.get("run_deseq2") else []

rule trim_fastq_files:
    input:
        fastq1 = config["fastq_dir"] + "/{sample}" + config['read1'] + config["fastq_ext"],
        fastq2 = config["fastq_dir"] + "/{sample}" + config['read2'] + config["fastq_ext"]
    output:
        filtered_qc_fastq1 = os.path.join(trimmed_fastq_files, "{sample}" + config['read1'] + config["fastq_ext"]),
        filtered_qc_fastq2 = os.path.join(trimmed_fastq_files, "{sample}" + config['read2'] + config["fastq_ext"])
    params:
        json_qc_file = os.path.join(trimmed_fastq_files, "fastp.json"),
        html_qc_file = os.path.join(trimmed_fastq_files, "fastp.html")
    message: "Running fastp filtering on {wildcards.sample}"
    #log: directories["filterqc_summary_directory"] + "/logs/{sample}.log"
    shell:
        "fastp -i {input.fastq1} -I {input.fastq2} -o {output.filtered_qc_fastq1} -O {output.filtered_qc_fastq2} -h {params.html_qc_file} -j {params.json_qc_file}" # > {log} 2>&1


rule mapping:
    input:
        fastq1 = os.path.join(trimmed_fastq_files, "{sample}" + config['read1'] + config["fastq_ext"]),
        fastq2 = os.path.join(trimmed_fastq_files, "{sample}" + config['read2'] + config["fastq_ext"])
    output:
        sam = os.path.join(sam_files, "{sample}.sam"),
        summary = os.path.join(sam_files, "{sample}_summary.txt"),
        metrics = os.path.join(sam_files, "{sample}_metrics.tsv")
    params:
        genome_reference_index = reference_dirs["genome_reference_index"]
    message: "Running Hisat2 alignment on {wildcards.sample}"
    threads: config["cores"]
    #log: directories["alignment_summary_directory"] + "/logs/{sample}.log"
    shell:
        "hisat2 -x {params.genome_reference_index} -p {threads} --new-summary --summary-file {output.summary} --met-file {output.metrics} -1 {input.fastq1} -2 {input.fastq2} -S {output.sam}" #" > {log} 2>&1"

rule convert_sam_to_bam:
    input:
        os.path.join(sam_files, "{sample}.sam")
    output:
        os.path.join(bam_files, "{sample}.bam")
    message: "Converting {wildcards.sample}.sam to {wildcards.sample}.bam"
    threads: config["cores"]
    shell:
        "samtools view -b -S -@ {threads} {input} | samtools sort -o {output}"

rule index_bam:
    input:
        os.path.join(bam_files, "{sample}.bam")
    output:
        os.path.join(bam_files, "{sample}.bam.bai")
    message: "Indexing {wildcards.sample}.bam"
    shell:
        "samtools index {input}"
        
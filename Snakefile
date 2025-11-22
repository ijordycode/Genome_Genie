import os
import pandas as pd


configfile: 'config.yaml'

directories = config['directories']
samplesheet_path = directories['samplesheet']
samples_df = pd.read_csv(samplesheet_path)
sample_list = samples_df['samples'].tolist()

if config.get("run_variant_calling", True):
    include: "Subworkflow/variant_calling.smk"

rule all:
    input:
        expand(os.path.join(directories["fastq_files_output"], "{sample}" + config['read1'] + config["fastq_ext"]), sample=sample_list),
        expand(os.path.join(directories["fastq_files_output"], "{sample}" + config['read2'] + config["fastq_ext"]), sample=sample_list),
        expand(os.path.join(directories['sam_files_output'], "{sample}.sam"), sample=sample_list),
        expand( os.path.join(directories["bam_files_output"], "{sample}.bam.bai"), sample=sample_list),
        expand(os.path.join(directories["vcf_files_output"], "{sample}.vcf"), sample=sample_list) if config.get("run_variant_calling") else [],
        expand(os.path.join(directories["count_files_output"], "{sample}_counts.txt"), sample=sample_list) if config.get("run_gene_quantification") else []

rule trim_fastq_files:
    input:
        fastq1 = config["fastq_dir"] + "/{sample}" + config['read1'] + config["fastq_ext"],
        fastq2 = config["fastq_dir"] + "/{sample}" + config['read2'] + config["fastq_ext"]
    output:
        filtered_qc_fastq1 = os.path.join(directories["fastq_files_output"], "{sample}" + config['read1'] + config["fastq_ext"]),
        filtered_qc_fastq2 = os.path.join(directories["fastq_files_output"], "{sample}" + config['read2'] + config["fastq_ext"]),
        #json_qc_file = os.path.join(directories["fastq_files_output"], "fastp.json"),
        #html_qc_file = os.path.join(directories["fastq_files_output"], "fastp.html")
    params:
        None
    message: "Running fastp filtering on {wildcards.sample}"
    #log: directories["filterqc_summary_directory"] + "/logs/{sample}.log"
    shell:
        "fastp -i {input.fastq1} -I {input.fastq2} -o {output.filtered_qc_fastq1} -O {output.filtered_qc_fastq2} -h {output.html} -j {output.json}" # > {log} 2>&1


rule mapping:
    input:
        fastq1 = os.path.join(directories["fastq_files_output"], "{sample}" + config['read1'] + config["fastq_ext"]),
        fastq2 = os.path.join(directories["fastq_files_output"], "{sample}" + config['read2'] + config["fastq_ext"])
    output:
        sam = os.path.join(directories["sam_files_output"], "{sample}.sam"),
        summary = os.path.join(directories["sam_files_output"], "{sample}_summary.txt"),
        metrics = os.path.join(directories["sam_files_output"], "{sample}_metrics.tsv")
    params:
        genome_reference_index = directories["genome_reference_index"]
    message: "Running Hisat2 alignment on {wildcards.sample}"
    threads: config["cores"]
    #log: directories["alignment_summary_directory"] + "/logs/{sample}.log"
    shell:
        "hisat2 -x {params.genome_reference_index} -p {threads} --new-summary --summary-file {output.summary} --met-file {output.metrics} -1 {input.fastq1} -2 {input.fastq2} -S {output.sam}" #" > {log} 2>&1"

rule convert_sam_to_bam:
    input:
        os.path.join(directories["sam_files_output"], "{sample}.sam")
    output:
        os.path.join(directories["bam_files_output"], "{sample}.bam")
    message: "Converting {wildcards.sample}.sam to {wildcards.sample}.bam"
    threads: config["cores"]
    shell:
        "samtools view -b -S -@ {threads} {input} | samtools sort -o {output}"

rule index_bam:
    input:
        os.path.join(directories["bam_files_output"], "{sample}.bam")
    output:
        os.path.join(directories["bam_files_output"], "{sample}.bam.bai")
    message: "Indexing {wildcards.sample}.bam"
    shell:
        "samtools index {input}"


"""

rule vcf_filtering:
    input:
        directories["vcf_directory"] + "/{sample}.vcf"
    output:
        directories["filtered_vcf_directory"] + "/{sample}.vcf"
    message: "Filtering {wildcards.sample}.vcf"
    run:
        parser = vcf_parser(str(input))
        options = config["vcf_filtering_options"]
        if "filter_chrom" in options and \
            (filter_chrom := options["filter_chrom"]):
            parser.slice_on_chrom(*filter_chrom)
        if "filter_info" in options and \
            (filter_info := options["filter_info"]):
            for key, ops in filter_info.items():
                for op in ("min", "max"):
                    if ops[op] != None:
                        parser.slice_on_info(key, op, ops[op])
        parser.dump(str(output))

rule effect_prediction:
    input:
        directories["filtered_vcf_directory"] + "/{sample}.vcf"
    output:
        data = directories["variant_effect_directory"] + "/{sample}.txt",
        summary = directories["variant_effect_summary_directory"] + "/{sample}.html"
    message: "Predicting variant effects for {wildcards.sample}"
    threads: 1
    log: directories["variant_effect_summary_directory"] + "/logs/{sample}.log"
    shell:
        "vep -i {input} -o {output.data} --format vcf --sf {output.summary} --database --symbol --tab --check_existing --fork {threads} > {log} 2>&1"

rule generate_pdf_summary:
    input:
        directories["variant_effect_directory"] + "/{sample}.txt"
    output:
        directories["pdf_summary_directory"] + "/{sample}.pdf"
    params:
        summaries_database = directories["summaries_database"],
        extra_options = get_options(config["vep_filtering_options"])
    message: "Generating PDF summary for {wildcards.sample}."
    log: directories["pdf_summary_directory"] + "/logs/{sample}.log"
    shell:
        "python3 helploid/generate_pdf.py -i {input} -o {output} -s {wildcards.sample} -db {params.summaries_database} {params.extra_options} > {log} 2>&1" """
        
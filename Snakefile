import os
import pandas as pd


configfile: 'config.yaml'

directories = config['directories']
samplesheet_path = directories['samplesheet']
samples_df = pd.read_csv(samplesheet_path)
sample_list = samples_df['samples'].tolist()

rule all:
    input:
        expand(os.path.join(directories["fastq_files_output"], "{sample}" + config['read1'] + config["fastq_ext"]), sample=sample_list),
        expand(os.path.join(directories["fastq_files_output"], "{sample}" + config['read2'] + config["fastq_ext"]), sample=sample_list),
        expand(os.path.join(directories['sam_files_output'], "{sample}.sam"), sample=sample_list)

rule trim_fastq_files:
    input:
        fastq1 = config["fastq_dir"] + "/{sample}" + config['read1'] + config["fastq_ext"],
        fastq2 = config["fastq_dir"] + "/{sample}" + config['read2'] + config["fastq_ext"]
    output:
        filtered_qc_fastq1 = os.path.join(directories["fastq_files_output"], "{sample}" + config['read1'] + config["fastq_ext"]),
        filtered_qc_fastq2 = os.path.join(directories["fastq_files_output"], "{sample}" + config['read2'] + config["fastq_ext"])
    params:
        None
    message: "Running fastp filtering on {wildcards.sample}"
    #log: directories["filterqc_summary_directory"] + "/logs/{sample}.log"
    shell:
        "fastp -i {input.fastq1} -I {input.fastq2} -o {output.filtered_qc_fastq1} -O {output.filtered_qc_fastq2}" # -h {output.html} -j {output.json} > {log} 2>&1


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

"""
def get_sample_suffix(read: int):
    if read == 1:
        return config["reads"][0]+'.fastq.gz'
    elif read == 2:
        return config["reads"][1]+'.fastq.gz'

rule all:
    input:
        expand("{dir}/{sample}{read}_fastqc.html", sample=config["samples"], read=config["reads"], dir=directories["pre_filterqc_directory"]),
        expand("{dir}/{sample}{read}_fastqc.html", sample=config["samples"], read=config["reads"], dir=directories["post_filterqc_directory"]),
        expand("{dir}/{sample}.pdf", sample=config["samples"], dir=directories["pdf_summary_directory"])


rule filterqc_fastq_files:
    input:
        fastq1 = directories["raw_fastq_files"] + "/{sample}" + get_sample_suffix(1),
        fastq2 = directories["raw_fastq_files"] + "/{sample}" + get_sample_suffix(2)
    output:
        filtered_qc_fastq1 = directories["filterqc_fastq_directory"] + "/{sample}" + get_sample_suffix(1),
        filtered_qc_fastq2 = directories["filterqc_fastq_directory"] + "/{sample}" + get_sample_suffix(2),
        html = directories["filterqc_summary_directory"] + "/{sample}.html",
        json = directories["filterqc_summary_directory"] + "/{sample}.json"
    params:
        extra_options = get_options(config["fastp"])
    message: "Running fastp filtering on {wildcards.sample}"
    log: directories["filterqc_summary_directory"] + "/logs/{sample}.log"
    shell:
        "fastp -i {input.fastq1} -I {input.fastq2} -o {output.filtered_qc_fastq1} -O {output.filtered_qc_fastq2} -h {output.html} -j {output.json} {params.extra_options} > {log} 2>&1"

 rule align_sample:
    input:
        fastq1 = directories["filterqc_fastq_directory"] + "/{sample}" + get_sample_suffix(1),
        fastq2 = directories["filterqc_fastq_directory"] + "/{sample}" + get_sample_suffix(2)
    output:
        sam = directories["aligned_sam_directory"] + "/{sample}.sam",
        summary = directories["alignment_summary_directory"] + "/{sample}_summary.txt",
        metrics = directories["alignment_summary_directory"] + "/{sample}_metrics.tsv"
    params:
        genome_reference_index = directories["genome_reference_index"]
    message: "Running Hisat2 alignment on {wildcards.sample}"
    threads: config["cores"]
    log: directories["alignment_summary_directory"] + "/logs/{sample}.log"
    shell:
        "hisat2 -x {params.genome_reference_index} -p {threads} --new-summary --summary-file {output.summary} --met-file {output.metrics} -1 {input.fastq1} -2 {input.fastq2} -S {output.sam} > {log} 2>&1"

rule convert_sam_to_bam:
    input:
        directories["aligned_sam_directory"] + "/{sample}.sam"
    output:
        directories["aligned_bam_directory"] + "/{sample}.bam"
    message: "Converting {wildcards.sample}.sam to {wildcards.sample}.bam"
    threads: config["cores"]
    shell:
        "samtools view -b -S -@ {threads} {input} | samtools sort -o {output}"

rule index_bam:
    input:
        directories["aligned_bam_directory"] + "/{sample}.bam"
    output:
        directories["aligned_bam_directory"] + "/{sample}.bam.bai"
    message: "Indexing {wildcards.sample}.bam"
    shell:
        "samtools index {input}"

rule variant_calling:
    input:
        directories["aligned_bam_directory"] + "/{sample}.bam.bai",
        bam = directories["aligned_bam_directory"] + "/{sample}.bam"
    output:
        directories["vcf_directory"] + "/{sample}.vcf"
    params:
        reference = directories["reference_genome_file"],
        extra_options = get_options(config["freebayes"])
    message: "Performing variant calling for {wildcards.sample}"
    log: directories["variant_calling_summary_directory"] + "/logs/{sample}.log"
    shell:
        "freebayes -f {params.reference} -v {output} {params.extra_options} {input.bam} > {log} 2>&1"

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
        
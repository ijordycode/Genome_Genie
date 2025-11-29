import os 


rule gene_counts:
    input:
        os.path.join(directories["bam_files_output"], "{sample}.bam.bai"),
        bam = os.path.join(directories["bam_files_output"], "{sample}.bam")
    output:
        counts = os.path.join(directories["count_files_output"], "{sample}.txt")
    params:
        gtf_file = directories["gtf_file"],
    message: "Performing gene countuing for {wildcards.sample}"
    #log: directories["variant_calling_summary_directory"] + "/logs/{sample}.log"
    shell:
        "featureCounts -p -t exon -g gene_id -a {params.gtf_file} -o {output.counts} {input.bam}"

rule run_deseq2:
    input:
        expand(os.path.join(directories["count_files_output"], "{sample}.txt"), sample=sample_list)
    output:
        os.path.join(directories["count_files_output"], "deseq2_results.csv")
    conda: "rbio"
    message: "Performing Deseq"
    #log: directories["variant_calling_summary_directory"] + "/logs/{sample}.log"
    shell:
        "Rscript Scripts/Deseq.R"


import os 


rule gene_counts:
    input:
        os.path.join(bam_files, "{sample}.bam.bai"),
        bam = os.path.join(bam_files, "{sample}.bam")
    output:
        counts = os.path.join(count_files, "{sample}.txt")
    params:
        gtf_file = reference_dirs["gtf_file"],
    message: "Performing gene countuing for {wildcards.sample}"
    #log: reference_dirs["variant_calling_summary_directory"] + "/logs/{sample}.log"
    shell:
        "featureCounts -p -t exon -g gene_id -a {params.gtf_file} -o {output.counts} {input.bam}"

rule run_deseq2:
    input:
        expand(os.path.join(count_files, "{sample}.txt"), sample=sample_list)
    output:
        os.path.join(count_files, "deseq2_results.csv")
    conda: "rbio"
    message: "Performing Deseq"
    #log: reference_dirs["variant_calling_summary_directory"] + "/logs/{sample}.log"
    shell:
        "Rscript Scripts/Deseq.R"


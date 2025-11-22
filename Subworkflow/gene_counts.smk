import os 


rule gene_counts:
    input:
        os.path.join(directories["bam_files_output"], "{sample}.bam.bai"),
        bam = os.path.join(directories["bam_files_output"], "{sample}.bam")
    output:
        os.path.join(directories["vcf_files_output"], "{sample}.vcf")
    params:
        gtf_file = directories["gtf_file"],
    message: "Performing gene countuing for {wildcards.sample}"
    #log: directories["variant_calling_summary_directory"] + "/logs/{sample}.log"
    shell:
        "featureCounts -t exon -g gene_id -a {params.gtf_file} -o {sample}_counts.txt {input.bam}"
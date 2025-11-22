import os 


rule variant_calling:
    input:
        os.path.join(directories["bam_files_output"], "{sample}.bam.bai"),
        bam = os.path.join(directories["bam_files_output"], "{sample}.bam")
    output:
        os.path.join(directories["vcf_files_output"], "{sample}.vcf")
    params:
        fasta_reference = directories["fasta_reference"],
    message: "Performing variant calling for {wildcards.sample}"
    #log: directories["variant_calling_summary_directory"] + "/logs/{sample}.log"
    shell:
        "freebayes -f {params.fasta_reference} -v {output} {input.bam}" #" > {log} 2>&1"
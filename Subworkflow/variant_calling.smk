import os 


rule variant_calling:
    input:
        os.path.join(bam_files, "{sample}.bam.bai"),
        bam = os.path.join(bam_files, "{sample}.bam")
    output:
        os.path.join(vcf_files, "{sample}.vcf")
    params:
        fasta_reference = reference_dirs["fasta_reference"],
    message: "Performing variant calling for {wildcards.sample}"
    #log: reference_dirs["variant_calling_summary_directory"] + "/logs/{sample}.log"
    shell:
        "freebayes -f {params.fasta_reference} -v {output} {input.bam}" #" > {log} 2>&1"
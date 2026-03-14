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

rule snpeff_annotation:
    input:
        os.path.join(vcf_files, "{sample}.vcf")
    output:
        os.path.join(vcf_files, "{sample}_annotated.vcf")
    params:
        snpeff_database = reference_dirs["snpeff_dir"],
    message: "Performing variant calling annotation for {wildcards.sample}"
    #log: reference_dirs["variant_calling_summary_directory"] + "/logs/{sample}.log"
    shell:
        "java -jar $(echo $CONDA_PREFIX/share/snpeff-*/snpEff.jar) -Xmx8g -dataDir {params.snpeff_database} {input} > {output}"



rule target:
    input:
        expand("data/snp_matrix/{lg}.012", lg=["LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10"])


# convert the vcf file to a matrix of 0, 1, 2 (homozygous ref, het, homozygous alternate)

rule snp_mat:
    input:
        "data/aurantiacus.vcf.gz"
    output:
        "data/snp_matrix/{lg}.012"
    log:
        "logs/snp_mat/{lg}.log"
    params:
        chrom="{lg}",
        out="data/snp_matrix/{lg}"
    resources:
        mem_mb_per_cpu=6000,
        partition="Standard"
    shell:
        """
        (vcftools --gzvcf {input} --012 --chr {params.chrom} \
        --out {params.out}) 2> {log}
        """

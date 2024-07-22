


rule target:
    input:
        expand("data/freqs/{popn}.frq", popn=["red", "yellow", "skt"])

# extract the red and yellow samples into separate vcfs to get allele frequencies

rule red_yellow:
    input:
        "data/aurantiacus.vcf.gz"
    output:
        "data/{popn}.vcf.gz"
    log:
        "logs/red_yellow/{popn}.log"
    resources:
        mem_mb_per_cpu=6000,
        partition="Standard"
    shell:
        """
        (bcftools view -Oz -S data/samps/{wildcards.popn}.samps {input} \
        -o {output}) 2> {log}
        """

# get frequencies of each site in the red and yellow groups

rule freqs:
    input:
        "data/{popn}.vcf.gz"
    output:
        "data/freqs/{popn}.frq"
    log:
        "logs/freqs/{popn}.log"
    params:
        "data/freqs/{popn}"
    resources:
        mem_mb_per_cpu=6000,
        partition="Standard"
    shell:
        """
        (vcftools --gzvcf {input} --freq --out {params}) 2> {log}
        """

# this runs a script to filter down to just divergent sites, the difference decided by "min" and "max"
# this returns a list with two columns, chromosome and position that will be able to be used by 
# bcftools to subset the final vcf

rule get_diffs:
    input:
        red="data/freqs/red.frq",
        yellow="data/freqs/yellow.frq",
        script="scripts/filter_freqs.sh"
    output:
        "data/freqs/diffs.list"
    log:
        "logs/get_diffs/log"
    params:
        min=".2",
        max=".8"
    resources:
        mem_mb_per_cpu=6000,
        partition="Standard"
    shell:
        """
        ({input.script} {params.min} {params.max} {input.red} {input.yellow} {output}) 2> {log}
        """

## this rule subsets the vcfs (red, yellow, and skt) so that only the diverged variants remain

rule vcf_subset:
    input:
        vcf_in="data/{popn}.vcf.gz",
        sites="data/freqs/diffs.list"
    output:
        "data/{popn}_diff.vcf.gz"
    log:
        "logs/vcf_subset/{popn}.log"
    resources:
        mem_mb_per_cpu=6000,
        partition="Standard"
    shell:
        """
        (bcftools view -Oz -R {input.sites} {input.vcf_in} \
        -o {output}) 2> {log}
        """


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

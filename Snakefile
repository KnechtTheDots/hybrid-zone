
# target changes depending on what the desired output is

rule target:
    input:
        expand("data/{popn}.vcf.gz", popn=["red","yellow","skt"])


# create a set of lists for reach population

rule get_samps:
    input:
        "data/aurantiacus.vcf.gz"
    output:
        expand("data/samps/{popn}.samps",popn=["red","yellow","skt"])
    log:
        "logs/get_samps/get.log"
    resources:
        mem_mb_per_cpu=6000,
        partition="Standard"
    shell:
        """
        (bcftools query -l {input} | awk '/PUN-R|CRS|UCSD_22/' > data/samps/red.samps 
        bcftools query -l {input} | awk '/INJ_22|PCT_2016|PUN-Y/' > data/samps/yellow.samps 
        bcftools query -l {input} | awk '/SKT/' > data/samps/skt.samps) 2> {log}
        """


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
        -o {output}
        bcftools index -t {output}) 2> {log}
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
        "div_sites/{popn}_diff.vcf.gz"
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


# convert the vcf file to a matrix of 0, 1, 2 (homozygous reference, het, homozygous alternate)

rule snp_mat:
    input:
        "div_sites/{popn}_diff.vcf.gz"
    output:
        "data/snp_matrix/{popn}.012"
    log:
        "logs/snp_mat/{popn}.log"
    params:
        out="data/snp_matrix/{popn}"
    resources:
        mem_mb_per_cpu=6000,
        partition="Standard"
    shell:
        """
        (vcftools --gzvcf {input} --012 --out {params.out}) 2> {log}
        """

# create a matrix that is scored not in comparison to the reference, but in comparison to 
# the ancestry of the genotype (red = 0, het = .5, yellow = 1)

rule ancestry_matrix:
    input:
        "data/snp_matrix/red.012",
        "data/snp_matrix/skt.012",
        script="scripts/make_ancestry_matrix.R"
    output:
        "data/ancestry_matrix.csv"
    log:
        "logs/ancestry_matrix/log"
    resources:
        mem_mb_per_cpu=6000,
        partition="Standard"
    shell:
        """
        (Rscript {input.script}) 2> {log}
        """

# correlation test of all pairwise sites between each set of chromosomes, this is
# stupid inefficient right now since each pair is being compared twice, but it runs pretty fast, I can 
# figure out a better way later

rule cor_test:
    input:
        data="data/ancestry_matrix.csv",
        diffs="data/freqs/diffs.list",
        script="scripts/cor_test.R"
    output:
        cors="cor_stats/correlations_{lg1}_{lg2}.csv",
        pvals="cor_stats/p_vals_{lg1}_{lg2}.csv"
    log:
        "logs/cor_test/{lg1}_{lg2}.log"
    params:
        lg1="{lg1}",
        lg2="{lg2}"
    resources:
        mem_mb_per_cpu=6000,
        partition="Standard"
    script:
        "scripts/cor_test.R"

# takes in each correlation matrix (and more importantly the associated p-values (-log10(pval)))
# and plots them

rule plot_cors:
    input:
        cors="cor_stats/correlations_{lg1}_{lg2}.csv",
        pvals="cor_stats/p_vals_{lg1}_{lg2}.csv",
        diffs="data/freqs/diffs.list"
    output:
        "figures/{lg1}_{lg2}_pval.jpeg"
    params:
        lg1="{lg1}",
        lg2="{lg2}"
    resources:
        mem_mb_per_cpu=6000,
        partition="Standard"
    script:
        "scripts/cor_plots.R"

# filters sites by correlation p-value and makes a plot with the entire genome on each axis and points
# at the sites with large -log10(p_value)

rule plot_all_sites:
    input:
        expand("cor_stats/p_vals_{lg1}_{lg2}.csv", lg1=["LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10"],
                lg2=["LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10"])
    output:
        "figures/all_sig.jpeg"
    log:
        "logs/all_sites/fig.log"
    resources:
        mem_mb_per_cpu=6000,
        partition="Standard"
    shell:
        """
        (Rscript scripts/all_sig_sites_plot.R) 2> {log}
        """
        

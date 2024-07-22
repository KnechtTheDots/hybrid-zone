## read the vcf into the environment

vcf <- VariantAnnotation::readVcf("data/aurantiacus.vcf.gz",
                                  "LG1")

## convert genotypes to matrix of 0,1,2 (homozygous ref, het, homozygous alternate)

snp_mat <- VariantAnnotation::genotypeToSnpMatrix(vcf)

## snp mat values are in raw, convert them to numeric

d <- snp_mat$genotypes@.Data

dd <- matrix(as.numeric(d), ncol = ncol(d), nrow = nrow(d))

colnames(dd) <- colnames(d)
rownames(dd) <- rownames(d)

# create csv with output

write.csv(dd, "data/snp_mat/LG1.csv")
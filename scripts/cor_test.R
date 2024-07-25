
# start with an if else statement, if the two chromosomes are the same, just write a nonsense
# data frame to not waste time, if they aren't, load in the data and run the analysis
# or if the chromosomes have already been evaluated (if the second lg is smaller than the first)


if(snakemake@params[['lg1']] == snakemake@params[['lg2']]){
  d <- data.frame("Comparison with self", "!")
  write.csv(d, snakemake@output[['cors']])
  write.csv(d, snakemake@output[['pvals']])
}else{

  
  # load in the ancestry matrix with red alleles coded as 0, hets as .5 and
  # yellow alleles coded as 1
  
  ancestry_matrix <- read.csv(snakemake@input[['data']])
  
  # load in the diffs data frame, the first column is the name of the linkage group
  diffs <- read.table(snakemake@input[['diffs']])
  
  # set of id's to index with later
  ids <- 1:ncol(ancestry_matrix)
  
  # pull out the ids of the first and second linkage groups
  lg1 <- ids[diffs[,1]==snakemake@params[['lg1']]]
  lg2 <- ids[diffs[,1]==snakemake@params[['lg2']]]
  
  
  # define two matrices to hold the correlations and pvalues
  cor_mat <- matrix(nrow = length(lg1), ncol = length(lg2))
  p_mat <- matrix(nrow = length(lg1), ncol = length(lg2))
  
  # calculate the correlations and p-values
  for(i in 1:length(lg1)){
    for(j in 1:length(lg2)){
      m <- cor.test(ancestry_matrix[,lg1[i]], ancestry_matrix[,lg2[j]])
      cor_mat[i,j] <- as.numeric(m$estimate)
      p_mat[i,j] <- -log10(m$p.value)
    }
  }
  
  
  # write the matrices to csv's
  
  write.csv(cor_mat, snakemake@output[['cors']], row.names = F, col.names = F)
  
  write.csv(p_mat, snakemake@output[['pvals']], row.names = F, col.names = F)

}
# load in the ancestry matrix with red alleles coded as 0, hets as .5 and
# yellow alleles coded as 1

ancestry_matrix <- read.csv("data/ancestry_matrix.csv")


# define the matrix to put correlations in, define the matrix (p_mat) to put
# -log10(p value) in

cor_mat <- matrix(nrow = ncol(ancestry_matrix), ncol = ncol(ancestry_matrix))
p_mat <- matrix(nrow = ncol(ancestry_matrix), ncol = ncol(ancestry_matrix))

# iterate over each column and run a cor.test with every other column I checked
# the ancestry matrix and there are no -1's, so I don't have to worry about that

for(i in 1:ncol(ancestry_matrix)){
  for(j in i:ncol(ancestry_matrix)){
    if(i==j){
      cor_mat[i,j] <- 1
      p_mat[i,j] <- -log10(1e-1000) # this is just to get infinite
    }else{
      m <- cor.test(ancestry_matrix[,i], ancestry_matrix[,j])
      cor_mat[i,j] <- cor_mat[j,i] <- as.numeric(m$estimate)
      p_mat[i,j] <- p_mat[j,i] <- -log10(m$p.value)
    }
  }
}

# write the matrices to csv's

write.csv(cor_mat, "cor_stats/correlations.csv", row.names = F, col.names = F)

write.csv(p_mat, "cor_stats/p_vals.csv", row.names = F, col.names = F)
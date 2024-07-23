# read in the red population snp matrix with 0,1,2 as reference homozygote, 
# het, and alternate homozygote. -1 indicates an NA

red_mat <- read.table("data/snp_matrix/red.012")

# the first column is rownames, so remove that

red_mat <- red_mat[,2:ncol(red_mat)]

# red mat has individuals are rows and sites as columns, so iterate over the columns
# to calculate the allele frequency in the reds, removing any that are -1 because
# this is an NA

freqs <- apply(red_mat, 2, function(x) sum(x[x!=-1])/(2*length(x[x!=-1])))

# if the frequency is less than or equal to .2 set it to 0, else set it to 2

ancestry <- as.numeric(ifelse(freqs <= .2, 0, 2))

# read in the skt matrix table

skt <- read.table("data/snp_matrix/skt.012")

# remove first column

skt <- skt[,2:ncol(skt)]

# recode skt so that 0 is the red ancestry, .5 is a het, and 1 is yellow ancestry

ancestry_matrix <- matrix(ncol = ncol(skt), nrow = nrow(skt))

for(i in 1:nrow(skt)){
  for(j in 1:ncol(skt)){
    ancestry_matrix[i,j] <- ifelse(skt[i,j]==ancestry[j], 0, 
                                   ifelse(skt[i,j]==1, .5,
                                          ifelse(skt[i,j]==-1, -1, 1)))
  }
}

# write the output to a csv

write.csv(ancestry_matrix, "data/ancestry_matrix.csv", row.names = F)
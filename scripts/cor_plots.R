# load in packages to make the plots
library(tidyverse)
library(cowplot)

# I don't compare chromosomes with themselves, so in order to deal with it and not 
# make snakemake crash, if lg1 is the same as lg2 I just make a histogram quick, else
# go on to making the plot

if(snakemake@params[['lg1']]==snakemake@params[['lg2']]){
  a <- data.frame(a = rnorm(100)) %>% 
    ggplot(aes(x = a)) +
    geom_histogram() +
    geom_vline(xintercept = 50)
  
  ggsave(snakemake@output[[1]], 
         plot = a,
         device = "jpeg")
  
}else{
  
  # cors is a matrix of correlations between all pairwise sites of lg1 and lg2, pvals
  # is the matrix of the -log10(pval) of those correlations
  cors <- read.csv(snakemake@input[['cors']])
  pvals <- read.csv(snakemake@input[['pvals']])
  
  # diffs is all sites that are diverged between red and yellow, I will use it to 
  # pull out the physical location of the site on each chromosome.
  diffs <- read.table(snakemake@input[['diffs']])
  
  # the first column is the chromosome "LG_" the second column is the physical location
  # on the chromosome
  colnames(diffs) <- c("chr", "site")
  
  # pull out the locations of all the diverged sites on lg1
  lg1 <- diffs$site[diffs$chr==snakemake@params[['lg1']]]
  # same for lg2
  lg2 <- diffs$site[diffs$chr==snakemake@params[['lg2']]]
  
  # make the columnames the location of the site
  colnames(cors) <- lg2
  colnames(pvals) <- lg2
  
  # plot the -log10(p-values) on lg1
  p_plot1 <- pvals %>% 
    # add a column that has the locations of lg1
    mutate(site_1 = lg1) %>% 
    # pivot so that each row is a site and site_1 is the location on lg1, site_2 is the location on lg2
    pivot_longer(1:length(lg2), names_to = "site_2", 
                 values_to = "pvals") %>% 
    # make sure everything is numeric, it was reading them wierd
    mutate(site_1 = as.numeric(site_1),
           site_2 = as.numeric(site_2),
           pvals = as.numeric(pvals)) %>% 
    # plot
    ggplot(aes(x = site_1, y = pvals)) +
    geom_point()
  
  # same as above but plotting the other chromosome
  p_plot2 <- pvals %>% 
    mutate(site_1 = lg1) %>% 
    pivot_longer(1:length(lg2), names_to = "site_2",
                 values_to = "pvals") %>% 
    mutate(site_1 = as.numeric(site_1),
           site_2 = as.numeric(site_2),
           pvals = as.numeric(pvals)) %>% 
    ggplot(aes(y = site_2, x = pvals)) +
    geom_point()
  
  # use cowplot to put them into a panel together, probably not the best way but ok for now
  p <- plot_grid(p_plot2, p_plot1, labels = c(snakemake@params[['lg2']],
                                              snakemake@params[['lg1']]))
  
  # save the plot
  ggsave(snakemake@output[[1]],
         plot = p,
         device = "jpeg", width = 12, height = 8)

}
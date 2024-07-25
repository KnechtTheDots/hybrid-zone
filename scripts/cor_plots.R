library(tidyverse)
library(cowplot)

if(snakemake@params[['lg1']]==snakemake@params[['lg2']]){
  a <- data.frame(a = rnorm(100)) %>% 
    ggplot(aes(x = a)) +
    geom_histogram() +
    geom_vline(xintercept = 50)
  
  ggsave(snakemake@output[[1]], 
         plot = a,
         device = "jpeg")
  
}else{

  cors <- read.csv(snakemake@input[['cors']])
  pvals <- read.csv(snakemake@input[['pvals']])
  diffs <- read.table(snakemake@input[['diffs']])
  
  colnames(diffs) <- c("chr", "site")
  
  lg1 <- diffs$site[diffs$chr==snakemake@params[['lg1']]]
  lg2 <- diffs$site[diffs$chr==snakemake@params[['lg2']]]
  
  colnames(cors) <- lg2
  colnames(pvals) <- lg2
  
  # cor_plot <- cors %>% 
  #   mutate(site_1 = lg1) %>% 
  #   pivot_longer(1:length(lg2), names_to = "site_2", values_to = "cors") %>% 
  #   mutate(site_1 = as.numeric(site_1),
  #          site_2 = as.numeric(site_2),
  #          cors = as.numeric(cors)) %>% 
  #   ggplot(aes(x = site_1, y = site_2)) +
  #   geom_point(aes(color = cors))
  # 
  # ggsave(snakemake@output[[1]], 
  #        plot = cor_plot,
  #        device = "jpeg")
  
  p_plot1 <- pvals %>% 
    mutate(site_1 = lg1) %>% 
    pivot_longer(1:length(lg2), names_to = "site_2",
                 values_to = "pvals") %>% 
    mutate(site_1 = as.numeric(site_1),
           site_2 = as.numeric(site_2),
           pvals = as.numeric(pvals)) %>% 
    ggplot(aes(x = site_1, y = pvals)) +
    geom_point()
  
  p_plot2 <- pvals %>% 
    mutate(site_1 = lg1) %>% 
    pivot_longer(1:length(lg2), names_to = "site_2",
                 values_to = "pvals") %>% 
    mutate(site_1 = as.numeric(site_1),
           site_2 = as.numeric(site_2),
           pvals = as.numeric(pvals)) %>% 
    ggplot(aes(y = site_2, x = pvals)) +
    geom_point()
  
  p <- plot_grid(p_plot2, p_plot1, labels = c(snakemake@params[['lg2']],
                                              snakemake@params[['lg1']]))
  
  ggsave(snakemake@output[[1]],
         plot = p,
         device = "jpeg", width = 12, height = 8)

}
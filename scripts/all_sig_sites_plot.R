library(tidyverse)


diffs <- read.table("data/freqs/diffs.list")
colnames(diffs) <- c("chr", "site")

lgs <- paste0("LG", 1:10)
sites_genome <- c(0, 20551961, 24346716, 20798979, 18557064,
                  25557270, 15765608, 17530972, 20267432,
                  14510045)

sites_tot <- cumsum(sites_genome)

d <- data.frame()
num <- 1:10
lg1 <- lg2 <- c()
for(i in lgs){
  lg1 <- as.numeric(diffs$site[diffs$chr==i]) + sites_tot[num[lgs==i]]
  pvals <- data.frame()
  for(j in lgs){
    if(i==j) next
    lg2 <- as.numeric(diffs$site[diffs$chr==j]) + sites_tot[num[lgs==j]]
    pvals <- read.csv(paste0("cor_stats/p_vals_",i,"_",j,".csv"))
    colnames(pvals) <- lg2
    pvals <- pvals %>% 
      mutate(site_1 = lg1) %>% 
      pivot_longer(1:length(lg2),
                   names_to = "site_2",
                   values_to = "pval") %>% 
      mutate(site_1 = as.numeric(site_1),
             site_2 = as.numeric(site_2),
             pval = as.numeric(pval)) %>% 
      filter(pval > 20) 
    d <- rbind(d, pvals)
  }
}

colnames(d) <- c("site_1", "site_2", "pval")

p <- d %>% 
  ggplot(aes(x = site_1, y = site_2)) +
  geom_point() +
  geom_hline(yintercept = sites_tot, color = "grey") +
  geom_vline(xintercept = sites_tot, color = "grey") +
  geom_abline(intercept = 0, slope = 1, color = "blue") +
  theme_minimal()

ggsave("figures/all_sig.jpeg", plot = p, device = "jpeg",
       width = 12, height = 8)
#import libraries
library(tidyverse)

fp <- "GSM2230760_marker_genes.csv"
#import the data
full_res <- read_csv(fp)
#select the top 10 adj. P-value genes from each cluster.
top10clusters <- full_res %>% group_by(cluster) %>% slice_min(p_val_adj, n= 10) %>% pull(gene)
#write values out to CSV file for enrichment in DAVID
top10 <- write.csv(top10clusters, "top10_gene_list.csv", row.names = FALSE, quote = FALSE)
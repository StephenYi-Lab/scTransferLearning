library(data.table)
library(tidyverse)
library(biomaRt)
setwd("/stor/work/SongYi/maya_y/scTransferLearning/scTransferLearning/Code/")

# https://ftp.ensembl.org/pub/current_mysql/ensembl_stable_ids_109/
gene_data <- fread("../Data/stable_id_lookup.txt")

is_human <- str_detect(gene_data$V1,"^ENSG[:digit:]")
human_genes <- gene_data[is_human,]

human_genes <- human_genes$V1

ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
# attributes <- listAttributes(ensembl)
# listFilters(ensembl)
final <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
               filters=c("ensembl_gene_id"),
               values = human_genes,
               mart = ensembl)

colnames(final) <- c("Gene.stable.ID","Gene.name")

write_csv(final,"../Data/hsaENSEMBL-GENES.txt")
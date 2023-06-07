############################
### hsaENSEMBL-GENES.txt ###
############################
# wget https://ftp.ensembl.org/pub/current_mysql/ensembl_stable_ids_109/
# converts ensembl_stable_ids to hsaENSEMBL-GENES.txt
source("make_ids.R")


########################################
### RAW.UMI.counts.BC.cell.lines.rds ###
########################################
# data information: https://figshare.com/articles/dataset/Single_Cell_Breast_Cancer_cell-line_Atlas/15022698?file=28893384
# wget https://figshare.com/ndownloader/files/28893384



################
### MCF7.csv ###
################
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114459/suppl/GSE114459_Polyak.csv.gz
# gunzip GSE114459_Polyak.csv.gz
# mv ./GSE114459_Polyak.csv ./MCF7.csv


############
### T47D ###
############
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4285nnn/GSM4285803/suppl/GSM4285803_scRNA_RawCounts.csv.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4285nnn/GSM4285803/suppl/GSM4285803_scRNA_metaInfo.csv.gz
# gunzip GSM4285803_scRNA_RawCounts.csv.gz
# gunzip GSM4285803_scRNA_metaInfo.csv.gz


### BT474 ###
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150949/suppl/GSE150949_pooled_watermelon.count.matrix.csv.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150949/suppl/GSE150949_pooled_watermelon.metadata.matrix.csv.gz
# gunzip GSE150949_pooled_watermelon.count.matrix.csv.gz
# gunzip GSE150949_pooled_watermelon.metadata.matrix.csv.gz


### count_matrix_sparse.mtx
### count_matrix_genes.tsv
### count_matrix_barcodes.tsv
### metadata.csv
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz
# tar -xf GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz



###################################################
### Original screen_All tissues_fitted data.csv ###
###################################################
# https://figshare.com/articles/dataset/Original_screen_drug_combination_data/16843597
# wget https://figshare.com/ndownloader/files/34006655
# mv ./34006655 "Original screen_All tissues_fitted data.csv"


###############
### SCP 542 ###
###############
# from:
# https://singlecell.broadinstitute.org/single_cell/study/SCP542/pan-cancer-cell-line-heterogeneity#study-download
# create a login and use "bulk download" 


###############################################
### CCLE_NP24.2009_Drug_data_2012.02.20.csv ###
###############################################
# https://www.maayanlab.net/LINCS/DCB/
# download source data
# dirs: drug senstivity data and cell-line mutation data/CCLE/CCLE_NP24.2009_Drug_data_2012.02.20.csv
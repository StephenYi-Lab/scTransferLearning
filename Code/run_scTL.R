library(tidyverse)
library(Seurat)
# wd is assumed to be in the code folder of this repo so data can be accessed properly 

################################################################################


# helper function for use_Seurat_in_scTL()
# maps cells to cell lines and returns cell_type_pred_knn and cell_type_pred_knn_prob
# takes:
#     yourData- counts
#     yourMetaData- metadata
#     orig_seurat- the original seurat object
use_scTL <- function(yourData,yourMetaData,orig_seurat){
  
  # following documentation from scTransferLearning repo
  load('../Results/refBRCA.RData')
  Q <- mapQuery(yourData, metadata_query = yourMetaData, ref_obj = clBRCA)
  Q <- knnPredict(Q, clBRCA, clBRCA$meta_data$cellLine, k = 5)
  
  # add metadata to seurat object
  orig_seurat <- AddMetaData(orig_seurat,col.name = "cell_type_pred_knn",Q$meta_data$cell_type_pred_knn) 
  orig_seurat <- AddMetaData(orig_seurat,col.name = "cell_type_pred_knn_prob",Q$meta_data$cell_type_pred_knn_prob) 
  
  # return original seurat object with attached meta data of predicted cell line 
  return(orig_seurat)
}


################################################################################


# helper function for use_Seurat_in_scTL()
clean_prop_data <- function(seurat_post_scTL){
  
  # take metadata from seurat_post_scTL object
  # organizes for downstream use
  seurat_post_scTL@meta.data %>% 
    group_by(cell_type_pred_knn) %>% 
    summarise(count = n()) %>% 
    as.data.frame() %>% 
    column_to_rownames("cell_type_pred_knn") %>% 
    mutate(prop = count/sum(count)) -> cell_props
  
  
  # reads framework for clean data
  # can be moved outside the function if needed for memory, but this file is <1kB
  input_tumorComposition <- read_rds("../Data/input_tumorComposition.rds")
  # makes copy
  input_tumorComposition_x <- input_tumorComposition
  # assigns input values into framework
  input_tumorComposition_x[rownames(cell_props)] <- cell_props$prop
  # return filled framework
  return(input_tumorComposition_x)
}


################################################################################


# takes in a seurat object and maps cells to cell lines
# with mapped information 
use_Seurat_in_scTL <- function(seurat_obj){
  
  # checking assumptions
  if(class(seurat_obj)!="Seurat"){
    stop("Your input class was not a Seurat object\nYour object class: ",class(seurat_obj))
  }
  
  
  # assumes active assay
  yourData_seurat_obj <- GetAssayData(seurat_obj)
  yourMetaData_seurat_obj <- seurat_obj[[]]
  
  # see helper function
  seurat_obj <- use_scTL(yourData_seurat_obj,
                         yourMetaData_seurat_obj,
                         seurat_obj)
  

  
  # see helper function doc above
  out_data <- clean_prop_data(seurat_obj)
  
  # returns cleaned data
  return(out_data)
}



# end of script
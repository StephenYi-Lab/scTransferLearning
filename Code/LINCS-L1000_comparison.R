BiocManager::install("ccdata")

library(ccdata)

# Moderated unbiased effect sizes values for 
# all 1309 drugs in the Connectivity Map build 02.
# A matrix where columns (1309) correspond to drugs and 
#                rows (13832) to gene symbols
data(cmap_es)


# Variances of unbiased effect sizes values for 
# all 1309 drugs in the Connectivity Map build 02.
# A matrix where columns correspond to drugs and rows to gene symbols
data(cmap_var)


#############################################
#############################################
#############################################



CCLE_drug_data <- read_csv("../Data/CCLE_NP24.2009_Drug_data_2012.02.20.csv")
data <- data.frame()
for( row in 1:nrow(CCLE_drug_data) ){  # row <- 1
  Cell_line <- CCLE_drug_data$`Primary Cell Line Name`[row]
  drug <- CCLE_drug_data$Compound[row]

  doses <- str_split(CCLE_drug_data$`Doses (uM)`[row],",")[[1]]
  activities <- str_split(CCLE_drug_data$`Activity Data (median)`[row],",")[[1]]
  activities_sd <- str_split(CCLE_drug_data$`Activity SD`[row],",")[[1]]
  
  if(length(doses)==length(activities)){
    out <- paste(Cell_line,drug,doses,activities,activities_sd,sep = ",")
  }else{
    stop("lengths don't match")
  }
  
  out <- separate(as.data.frame(out),
                  col="out",
                  into = c("Cell_line","drug","doses","activities","activities_sd"),
                  sep = ",")
  data <- rbind(data,out)
  
  
}

data$doses <- as.numeric(data$doses)
data$activities <- as.numeric(data$activities)
data$activities_sd <- as.numeric(data$activities_sd)


data %>% str()




data[tumorComposition %in% data$Cell_line,]


data$activities %>% mean()

# CV with LINCS data
# show low accuracy

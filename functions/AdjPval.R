#!/usr/bin/env Rscript

AdjustedPval<-function(filePath){
  # function calculates Adjusted P values. It calculates bonferroni, BH,
  # holm, hochberge, hommel and BY, assuming P-value is in column "P.value".
  options(scipen = 999)
    
  data<-read.table(filePath, header = TRUE) # read file
  sortedData =data[order(data$p_value),] # order by P-Value
  sortedData$Bonferroni=p.adjust(sortedData$p_value,method = "bonferroni")
  sortedData$Holm =p.adjust(sortedData$p_value,method = "holm")
  sortedData$BH =p.adjust(sortedData$p_value,method = "BH")
  sortedData$Hochberg =p.adjust(sortedData$p_value,method = "hochberg")
  sortedData$Hommel =p.adjust(sortedData$p_value,method = "hommel")
  sortedData$BY =p.adjust(sortedData$p_value,method = "BY")
  sortedData$fdr =p.adjust(sortedData$p_value,method = "fdr")

  write.table(sortedData,saveFile,row.names = F,sep="\t",quote = F) # Save to new file, by adding new columns of Adj. P-values
}

###############################################################################

# Pass input and outpute filenames as arguments
args <- commandArgs(trailingOnly = TRUE)
filePath <- args[1]
saveFile <- args[2]

dataWithAdjPval<-AdjustedPval(filePath)

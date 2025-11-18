#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2025
#
#

library(openxlsx)
# read in iBAQ, get columnMean write to rds
myiBAQ <- read.xlsx("../resources/IBAQ_SA6850_myContrst.xlsx", sheet = "Sheet1")
str(myiBAQ)

qCols <- grep(pattern = "SA6850_", x = colnames(myiBAQ))
iBDF <- data.frame(protein_Id = myiBAQ[,1])
iBDF$mean_iBAQ <- rowMeans(myiBAQ[,qCols], na.rm = TRUE)

# read in results from DEA write to rds
myDEAres <- read.xlsx("../resources/DE_WUSA6850_myContrst.xlsx", sheet = "diff_exp_analysis_wide")

# join DEA n iBAQ
myDEAres <- left_join(x = myDEAres, y = iBDF)

# write to rds
write_rds(x = myDEAres, file = "../resources/SA6850_prXallWide.rds")

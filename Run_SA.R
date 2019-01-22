##### Run sensitivity analyses
rm(list=ls(all=TRUE))

# NAME THE RESULTS FILE
results_file <- paste0("SensAnalysis_", gsub(" ", "_", paste0(substr(date(), 5, 10), "_", substr(date(), 21, 24))))

# loading packages
library(xlsx) # working with excel
library(lmerTest) # mixed models
library(boot) # bootstrapping
library(quantreg) # quantile regression

##### set working directory
setwd("//ads.bris.ac.uk/filestore/BRMS/Research/Metabolic profiles in UPBEAT/RScripts/SensitivityAnalyses/")

##### source functions for each of the various sensitivity analyses
source('Functions_forSA.R')

##### read in main result and write to excel file
MainResult <- read.csv("MainResult.csv")

write.xlsx(MainResult, file=paste0(results_file, ".xlsx"), sheetName="Main Analyses")

##### MLM variations
SA_MLM_variations("SD", RemoveOutliers=FALSE, results_file)
SA_MLM_variations("SD", RemoveOutliers=TRUE, results_file)
SA_MLM_variations("IQR", RemoveOutliers=FALSE, results_file)
SA_MLM_variations("IQR", RemoveOutliers=TRUE, results_file)


##### ttest bootstrapped variations
ttest_bootstrap_variations("SD", RemoveOutliers=FALSE, results_file)
ttest_bootstrap_variations("SD", RemoveOutliers=TRUE, results_file)
ttest_bootstrap_variations("IQR", RemoveOutliers=FALSE, results_file)
ttest_bootstrap_variations("IQR", RemoveOutliers=TRUE, results_file)
MADS(RemoveOutliers=FALSE, results_file)
MADS(RemoveOutliers=TRUE, results_file)


##### Quantile regression
qreg_diffcentile(0.5, RemoveOutliers=FALSE, results_file=results_file)
qreg_diffcentile(0.5, RemoveOutliers=TRUE, MADprop=3.5, results_file=results_file)
qreg_diffcentile(0.75, RemoveOutliers=FALSE, results_file=results_file)
qreg_diffcentile(0.75, RemoveOutliers=TRUE, MADprop=3.5, results_file=results_file)


##### Proportion of outliers
Proportion_outliers(2.24, results_file=results_file)
Proportion_outliers(3.5, results_file=results_file)

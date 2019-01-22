### functions for SA 

## Remove outliers - method 1, removing top and bottom 1%
remove_outliers_1pc <- function(data_touse, variable, timepoint, treatment){
  if (treatment=NA){
    quant <- quantile(data_touse[which(data_touse$time_ind==timepoint), variable], na.rm=TRUE, probs=c(0.01, 0.99), type=2) # type determined by http://data.princeton.edu/stata/markdown/quantiles.htm
    data_touse[which(data_touse[, variable]<quant["1%"] & data_touse$time_ind==timepoint), variable] <- NA
    data_touse[which(data_touse[, variable]>quant["99%"] & data_touse$time_ind==timepoint), variable] <- NA
  } else {
    quant <- quantile(data_touse[which(data_touse$time_ind==timepoint & data_touse$Treatment==treatment), variable], na.rm=TRUE, probs=c(0.01, 0.99), type=2) # type determined by http://data.princeton.edu/stata/markdown/quantiles.htm
    data_touse[which(data_touse[, variable]<quant["1%"] & data_touse$time_ind==timepoint & data_touse$Treatment==treatment), variable] <- NA
    data_touse[which(data_touse[, variable]>quant["99%"] & data_touse$time_ind==timepoint & data_touse$Treatment==treatment), variable] <- NA
  }
  return(data_touse)
}

## function to run MLM models
SA_MLM_variations <- function(SD_IQR, RemoveOutliers, results_file){
  
  ### load and prepare data
  data_touse <- read.csv("//ads.bris.ac.uk/filestore/BRMS/Research/Metabolic profiles in UPBEAT/Data/upbeat_withmet_data_longSTATA_v2.csv", na.strings='.')
  
  ## force columns into factors
  make_factor <- c("BMI", "Ethnicity", "Parity", "Centre")
  for (m in make_factor){
    data_touse[, m] <- as.factor(data_touse[, m])
  }
  
  # create new time variable
  data_touse$time=data_touse$time_point/4
  
  # create new interaction variable of treatment*time
  data_touse$inter=data_touse$Treatment*data_touse$time
  
  # name the metabolites
  varlist <- colnames(data_touse[11:168])
  
  # remove outliers if required
  if (RemoveOutliers==TRUE){
    for (i in 1:length(varlist)){
      data_touse <- remove_outliers_1pc(data_touse, varlist[i], timepoint=15, treatment=NA)
      data_touse <- remove_outliers_1pc(data_touse, varlist[i], timepoint=28, treatment=0)
      data_touse <- remove_outliers_1pc(data_touse, varlist[i], timepoint=28, treatment=1)
      data_touse <- remove_outliers_1pc(data_touse, varlist[i], timepoint=36, treatment=0)
      data_touse <- remove_outliers_1pc(data_touse, varlist[i], timepoint=36, treatment=1)
    }
  }
  
  # scale data using SD or IQR accordingly
  w <- which(data_touse$time_ind==15)
  if (SD_IQR=="SD"){
    ### convert data using SD  
    for (i in 1:length(varlist)){
      sd_var <- sd(data_touse[w, varlist[i]], na.rm=TRUE)
      data_touse[, varlist[i]] <- data_touse[, varlist[i]] / sd_var
    }
  } else if (SD_IQR=="IQR"){
    ### convert data using IQR
    for (i in 1:length(varlist)){
      quant <- quantile(data_touse[w, varlist[i]], na.rm=TRUE, probs=c(0.25, 0.75), type=2) # type determined by http://data.princeton.edu/stata/markdown/quantiles.htm
      data_touse[, varlist[i]] <- data_touse[, varlist[i]] /((quant["75%"]-quant["25%"])/1.34)
    }
  } 
  
  # create results table
  results <- data.frame(matrix(NA, nrow=length(varlist), ncol=4), row.names = varlist)
  colnames(results) <- c("zestmlm", "coeffmlm", "secoeffmlm", "p")
  
  # run the mixed models
  # note different metabolites need different parameters (see kate's do file)
  # can we use actual model from paper instead of this version - should be the same
  for (i in 1:length(varlist)){
    
    #### run the model
    mod <- lmer(as.formula(paste0(varlist[i], " ~ BMI + Ethnicity + Parity + Age + Centre + time_point + Treatment:time_point + (time_point|participant_id)")), data=data_touse)

    #### save results
    su <- summary(mod)
    results[varlist[i], "zestmlm"] <- su$coefficients["time_point:Treatment", "Estimate"] / su$coefficients["time_point:Treatment", "Std. Error"]
    results[varlist[i], "coeffmlm"] <- su$coefficients["time_point:Treatment", "Estimate"] 
    results[varlist[i], "secoeffmlm"] <- su$coefficients["time_point:Treatment", "Std. Error"]
    
  }
  
  ######## add p-value
  results$p[results$zestmlm>0] <- 2*(1-pnorm(results$zestmlm[results$zestmlm>0]))
  results$p[results$zestmlm<=0] <- 2*pnorm(results$zestmlm[results$zestmlm<=0])
  
  Results_tosave <- results[, c("coeffmlm", "secoeffmlm", "p")]
  names(Results_tosave) <- c("Coeff", "SE", "P-value")
  
  ## write to file
  if (RemoveOutliers==TRUE){Out = ", No Outliers"} else {Out = ""}
  write.xlsx(Results_tosave, file=paste0(results_file, ".xlsx"), sheetName=paste0("MLM, ", SD_IQR, Out), append=TRUE)
  
}

# lm function used in the bootstrap
lm_diff_treatment <- function(formula, data, indices){
  d <- data[indices,] # allows boot to select sample 
  fit <- lm(formula, data=d)
  su_fit <- summary(fit)
  ret <- c("coeff"=su_fit$coefficients["Treatment", "Estimate"], 
           "se"=su_fit$coefficients["Treatment", "Std. Error"], 
           "z"=su_fit$coefficients["Treatment", "Estimate"]/su_fit$coefficients["Treatment", "Std. Error"])
  return(ret)
}

# function for ttest bootstrap models
ttest_bootstrap_variations <- function(SD_IQR, RemoveOutliers, results_file){
  
  ### load and prepare data
  data_touse <- read.csv("//ads.bris.ac.uk/filestore/BRMS/Research/Metabolic profiles in UPBEAT/Data/upbeat_withmet_data_longSTATA_v2.csv", na.strings='.')
  
  # delete time_ind=28
  data_touse <- data_touse[-which(data_touse$time_ind==28), ]
  
  # create new timediff variable
  timediff=(36-15)/4
  
  # name the metabolites
  varlist <- colnames(data_touse[11:168])
  
  # remove outliers if required
  if (RemoveOutliers==TRUE){
    for (i in 1:length(varlist)){
      data_touse <- remove_outliers_1pc(data_touse, varlist[i], timepoint=15, treatment=NA)
      data_touse <- remove_outliers_1pc(data_touse, varlist[i], timepoint=36, treatment=0)
      data_touse <- remove_outliers_1pc(data_touse, varlist[i], timepoint=36, treatment=1)
    }
  }
  
  # scale data using SD, IQR
  w <- which(data_touse$time_ind==15)
  if (SD_IQR=="SD"){
    ### convert data using sd
    for (i in 1:length(varlist)){
      
      sd_var <- sd(data_touse[w, varlist[i]], na.rm=TRUE)
      data_touse[, paste0(varlist[i], "sd")] <- data_touse[, varlist[i]] / sd_var
      
      for (p_id in unique(data_touse$participant_id)){
        sd_at15 <- data_touse[which(data_touse$time_ind==15 & data_touse$participant_id==p_id), paste0(varlist[i], "sd")]
        if (is.na(sd_at15)==FALSE){
          data_touse[which(data_touse$time_ind==36 & data_touse$participant_id==p_id), paste0(varlist[i], "diffsd")] <- 
            (data_touse[which(data_touse$time_ind==36 & data_touse$participant_id==p_id), paste0(varlist[i], "sd")] - sd_at15) / timediff
        }
      }
    } 
    diff_using="diffsd"
  } else if (SD_IQR=="IQR"){
    ### convert data using IQR
    for (i in 1:length(varlist)){
      
      quant <- quantile(data_touse[w, varlist[i]], na.rm=TRUE, probs=c(0.25, 0.75), type=2) # type determined by http://data.princeton.edu/stata/markdown/quantiles.htm
      data_touse[, paste0(varlist[i], "iqr")] <- data_touse[, varlist[i]] /((quant["75%"]-quant["25%"])/1.34)
      
      for (p_id in unique(data_touse$participant_id)){
        iqr_at15 <- data_touse[which(data_touse$time_ind==15 & data_touse$participant_id==p_id), paste0(varlist[i], "iqr")]
        if (is.na(iqr_at15)==FALSE){
          data_touse[which(data_touse$time_ind==36 & data_touse$participant_id==p_id), paste0(varlist[i], "diffiqr")] <- 
            (data_touse[which(data_touse$time_ind==36 & data_touse$participant_id==p_id), paste0(varlist[i], "iqr")] - iqr_at15) / timediff
        }
      }
    }
    diff_using="diffiqr"
  }
  
  # reduce dataset to just time_ind 36
  data_touse <- data_touse[data_touse$time_ind==36, ]
  
  # create results table
  results <- data.frame(matrix(NA, nrow=length(varlist), ncol=4), row.names = varlist)
  colnames(results) <- c("z", "coeff", "se", "p")
  
  
  for (i in 1:length(varlist)){
    
    res <- boot(data=data_touse, statistic=lm_diff_treatment, R=1000, formula = as.formula(paste0(varlist[i], diff_using, " ~ Treatment")))
    # function lm_diff_treatment defined in Functions_forSA.R file
    # returns t1=coeff, t2=sd, t3=z
    
    #### save results
    results[varlist[i], "z"] <- res$t0["z"]
    results[varlist[i], "coeff"] <- res$t0["coeff"]
    results[varlist[i], "se"] <- res$t0["se"]
  }
  
  # calculate p value
  results$p[results$z>0] <- 2*(1-pnorm(results$z[results$z>0]))
  results$p[results$z<=0] <- 2*pnorm(results$z[results$z<=0])
  
  Results_tosave <- results[, c("coeff", "se", "p")]
  names(Results_tosave) <- c("Coeff", "SE", "P-value")
  
  ## write to file
  if (RemoveOutliers==TRUE){Out = ", NoOutlier"} else {Out = ""}
  write.xlsx(Results_tosave, file=paste0(results_file, ".xlsx"), sheetName=paste0("P-t-test bstrap, ", SD_IQR, Out), append=TRUE)
  
}

# MADS method
MADS <- function(RemoveOutliers, results_file){
  
  ### load and prepare data
  data_touse <- read.csv("//ads.bris.ac.uk/filestore/BRMS/Research/Metabolic profiles in UPBEAT/Data/upbeat_withmet_data_longSTATA_v2.csv", na.strings='.')
  
  # create new timediff variable
  timediff=(36-15)/4
  
  # name the metabolites
  varlist <- colnames(data_touse[11:168])
  
  # remove outliers if required
  if (RemoveOutliers==TRUE){
    MADprop=3.5
    for (i in 1:length(varlist)){
      ### remove outliers
      med_var <- vector(len=nrow(data_touse))
      
      # time 15
      med_var[which(data_touse$time_ind==15)] <- median(data_touse[which(data_touse$time_ind==15), varlist[i]], na.rm=TRUE) 
      devvar <- abs(data_touse[, varlist[i]] - med_var)[which(data_touse$time_ind==15)]
      data_touse[which(data_touse$time_ind==15), paste0(varlist[i], "mad")] <- median(devvar, na.rm=TRUE)
      
      # time 28 and treatment 0
      med_var[which(data_touse$time_ind==28 & data_touse$Treatment==0)] <- median(data_touse[which(data_touse$time_ind==28 & data_touse$Treatment==0), varlist[i]], na.rm=TRUE) 
      devvar <- abs(data_touse[, varlist[i]] - med_var)[which(data_touse$time_ind==28 & data_touse$Treatment==0)]
      data_touse[which(data_touse$time_ind==28 & data_touse$Treatment==0), paste0(varlist[i], "mad")] <- median(devvar, na.rm=TRUE)
      
      # time 28 and treatment 1
      med_var[which(data_touse$time_ind==28 & data_touse$Treatment==1)] <- median(data_touse[which(data_touse$time_ind==28 & data_touse$Treatment==1), varlist[i]], na.rm=TRUE) 
      devvar <- abs(data_touse[, varlist[i]] - med_var)[which(data_touse$time_ind==28 & data_touse$Treatment==1)]
      data_touse[which(data_touse$time_ind==28 & data_touse$Treatment==1), paste0(varlist[i], "mad")] <- median(devvar, na.rm=TRUE)
      
      # time 36 and treatment 0
      med_var[which(data_touse$time_ind==36 & data_touse$Treatment==0)] <- median(data_touse[which(data_touse$time_ind==36 & data_touse$Treatment==0), varlist[i]], na.rm=TRUE) 
      devvar <- abs(data_touse[, varlist[i]] - med_var)[which(data_touse$time_ind==36 & data_touse$Treatment==0)]
      data_touse[which(data_touse$time_ind==36 & data_touse$Treatment==0), paste0(varlist[i], "mad")] <- median(devvar, na.rm=TRUE)
      
      # time 36 and treatment 1
      med_var[which(data_touse$time_ind==36 & data_touse$Treatment==1)] <- median(data_touse[which(data_touse$time_ind==36 & data_touse$Treatment==1), varlist[i]], na.rm=TRUE) 
      devvar <- abs(data_touse[, varlist[i]] - med_var)[which(data_touse$time_ind==36 & data_touse$Treatment==1)]
      data_touse[which(data_touse$time_ind==36 & data_touse$Treatment==1), paste0(varlist[i], "mad")] <- median(devvar, na.rm=TRUE)
      
      # mark outliers
      data_touse[ , paste0(varlist[i], "outlier")] <- 0
      data_touse[which(data_touse[ , varlist[i]] >= (med_var + MADprop*data_touse[, paste0(varlist[i], "mad")])), paste0(varlist[i], "outlier")] <- 1
      data_touse[which(data_touse[ , varlist[i]] <= (med_var - MADprop*data_touse[, paste0(varlist[i], "mad")])), paste0(varlist[i], "outlier")] <- 1
      data_touse[which(is.na(data_touse[ , varlist[i]])) , paste0(varlist[i], "outlier")] <- NA
      
      # remove outliers
      data_touse[which(data_touse[ , paste0(varlist[i], "outlier")]==1), varlist[i]] <- NA
    } 
  }
  
  # scale data 
  w <- which(data_touse$time_ind==15)
  ### convert data using mean observed deviation
  for (i in 1:length(varlist)){
    
    med_var_o <- median(data_touse[which(data_touse$time_ind==15), varlist[i]], na.rm=TRUE) 
    med_dev_var_o <- median(abs(data_touse[which(data_touse$time_ind==15), varlist[i]] - med_var_o), na.rm=TRUE)
    
    data_touse[, paste0(varlist[i], "madsc")] <- data_touse[, varlist[i]] / (med_dev_var_o*1.4826)
    
    for (p_id in unique(data_touse$participant_id)){
      mad_at15 <- data_touse[which(data_touse$time_ind==15 & data_touse$participant_id==p_id), paste0(varlist[i], "madsc")]
      if (is.na(mad_at15)==FALSE){
        data_touse[which(data_touse$time_ind==36 & data_touse$participant_id==p_id), paste0(varlist[i], "diffmad")] <-
          (data_touse[which(data_touse$time_ind==36 & data_touse$participant_id==p_id), paste0(varlist[i], "madsc")] - mad_at15) / timediff
      }
    }
    
  }
  diff_using="diffmad"
  
  # reduce dataset to just time_ind 36
  data_touse <- data_touse[data_touse$time_ind==36, ]
  
  # create results table
  results <- data.frame(matrix(NA, nrow=length(varlist), ncol=4), row.names = varlist)
  colnames(results) <- c("z", "coeff", "se", "p")
  
  
  for (i in 1:length(varlist)){
    
    res <- boot(data=data_touse, statistic=lm_diff_treatment, R=100, formula = as.formula(paste0(varlist[i], diff_using, " ~ Treatment")))
    # function lm_diff_treatment defined in Functions_forSA.R file
    # returns t1=coeff, t2=sd, t3=z
    
    #### save results
    results[varlist[i], "z"] <- res$t0["z"]
    results[varlist[i], "coeff"] <- res$t0["coeff"]
    results[varlist[i], "se"] <- res$t0["se"]
  }
  
  # calculate p value
  results$p[results$z>0] <- 2*(1-pnorm(results$z[results$z>0]))
  results$p[results$z<=0] <- 2*pnorm(results$z[results$z<=0])
  
  Results_tosave <- results[, c("coeff", "se", "p")]
  names(Results_tosave) <- c("Coeff", "SE", "P-value")
  
  ## write to file
  if (RemoveOutliers==TRUE){Out = ", NoOutlier"} else {Out = ""}
  write.xlsx(Results_tosave, file=paste0(results_file, ".xlsx"), sheetName=paste0("P-t-test bstrap, MAD", Out), append=TRUE)
  
}

# Qreg method
qreg_diffcentile <- function(centile, RemoveOutliers, MADprop, results_file){
  
  ### load and prepare data
  data_touse <- read.csv("//ads.bris.ac.uk/filestore/BRMS/Research/Metabolic profiles in UPBEAT/Data/upbeat_withmet_data_longSTATA_v2.csv", na.strings='.')
  
  # create new timediff variable
  timediff=(36-15)/4
  
  # name the metabolites
  varlist <- colnames(data_touse[11:168])
  
  # remove outliers if required using MADs
  if (RemoveOutliers==TRUE){
    for (i in 1:length(varlist)){
      ### remove outliers
      med_var <- vector(len=nrow(data_touse))
      
      # time 15
      med_var[which(data_touse$time_ind==15)] <- median(data_touse[which(data_touse$time_ind==15), varlist[i]], na.rm=TRUE) 
      devvar <- abs(data_touse[, varlist[i]] - med_var)[which(data_touse$time_ind==15)]
      data_touse[which(data_touse$time_ind==15), paste0(varlist[i], "mad")] <- median(devvar, na.rm=TRUE)
      
      # time 28 and treatment 0
      med_var[which(data_touse$time_ind==28 & data_touse$Treatment==0)] <- median(data_touse[which(data_touse$time_ind==28 & data_touse$Treatment==0), varlist[i]], na.rm=TRUE) 
      devvar <- abs(data_touse[, varlist[i]] - med_var)[which(data_touse$time_ind==28 & data_touse$Treatment==0)]
      data_touse[which(data_touse$time_ind==28 & data_touse$Treatment==0), paste0(varlist[i], "mad")] <- median(devvar, na.rm=TRUE)
      
      # time 28 and treatment 1
      med_var[which(data_touse$time_ind==28 & data_touse$Treatment==1)] <- median(data_touse[which(data_touse$time_ind==28 & data_touse$Treatment==1), varlist[i]], na.rm=TRUE) 
      devvar <- abs(data_touse[, varlist[i]] - med_var)[which(data_touse$time_ind==28 & data_touse$Treatment==1)]
      data_touse[which(data_touse$time_ind==28 & data_touse$Treatment==1), paste0(varlist[i], "mad")] <- median(devvar, na.rm=TRUE)
      
      # time 36 and treatment 0
      med_var[which(data_touse$time_ind==36 & data_touse$Treatment==0)] <- median(data_touse[which(data_touse$time_ind==36 & data_touse$Treatment==0), varlist[i]], na.rm=TRUE) 
      devvar <- abs(data_touse[, varlist[i]] - med_var)[which(data_touse$time_ind==36 & data_touse$Treatment==0)]
      data_touse[which(data_touse$time_ind==36 & data_touse$Treatment==0), paste0(varlist[i], "mad")] <- median(devvar, na.rm=TRUE)
      
      # time 36 and treatment 1
      med_var[which(data_touse$time_ind==36 & data_touse$Treatment==1)] <- median(data_touse[which(data_touse$time_ind==36 & data_touse$Treatment==1), varlist[i]], na.rm=TRUE) 
      devvar <- abs(data_touse[, varlist[i]] - med_var)[which(data_touse$time_ind==36 & data_touse$Treatment==1)]
      data_touse[which(data_touse$time_ind==36 & data_touse$Treatment==1), paste0(varlist[i], "mad")] <- median(devvar, na.rm=TRUE)
      
      data_touse[ , paste0(varlist[i], "outlier")] <- 0
      data_touse[which(data_touse[ , varlist[i]] >= (med_var + MADprop*data_touse[, paste0(varlist[i], "mad")])), paste0(varlist[i], "outlier")] <- 1
      data_touse[which(data_touse[ , varlist[i]] <= (med_var - MADprop*data_touse[, paste0(varlist[i], "mad")])), paste0(varlist[i], "outlier")] <- 1
      
      # remove outliers
      data_touse[data_touse[ , paste0(varlist[i], "outlier")]==1, varlist[i]] <- NA
      
    } 
  }
  
  # scale data
  w <- which(data_touse$time_ind==15)
  mad_vars <- vector(len=length(varlist))
  for (i in 1:length(varlist)){
    
    med_var_o <- median(data_touse[which(data_touse$time_ind==15), varlist[i]], na.rm=TRUE) 
    med_dev_var_o <- median(abs(data_touse[which(data_touse$time_ind==15), varlist[i]] - med_var_o), na.rm=TRUE)
    
    data_touse[, paste0(varlist[i], "madsc")] <- data_touse[, varlist[i]] / (med_dev_var_o*1.4826)
    
    for (p_id in unique(data_touse$participant_id)){
      mad_at15 <- data_touse[which(data_touse$time_ind==15 & data_touse$participant_id==p_id), paste0(varlist[i], "madsc")]
      if (is.na(mad_at15)==FALSE){
        data_touse[which(data_touse$time_ind==36 & data_touse$participant_id==p_id), paste0(varlist[i], "diffmad")] <- 
          (data_touse[which(data_touse$time_ind==36 & data_touse$participant_id==p_id), paste0(varlist[i], "madsc")] - mad_at15) / timediff
      }
    }
  }
  
  
  # reduce dataset to just time_ind 36
  data_touse <- data_touse[data_touse$time_ind==36, ]
  
  # create results table
  results <- data.frame(matrix(NA, nrow=length(varlist), ncol=4), row.names = varlist)
  colnames(results) <- c("z", "coeff", "se", "p")
  
  for (i in 1:length(varlist)){
    
    # quantile regression
    res <- rq(as.formula(paste0(varlist[i], "diffmad" , "~ Treatment")), data=data_touse, tau = centile)
    su=summary(res, se="iid")
    
    #### save results
    results[varlist[i], "z"] <- su$coefficients["Treatment", "t value"] # su$coefficients["Treatment", "Value"]/su$coefficients["Treatment", "Std. Error"]
    results[varlist[i], "coeff"] <- su$coefficients["Treatment", "Value"]
    results[varlist[i], "se"] <- su$coefficients["Treatment", "Std. Error"]
  }
  
  # calculate p value
  results$p[results$z>0] <- 2*(1-pnorm(results$z[results$z>0]))
  results$p[results$z<=0] <- 2*pnorm(results$z[results$z<=0])
  
  Results_tosave <- results[, c("coeff", "se", "p")]
  names(Results_tosave) <- c("Coeff", "SE", "P-value")
  
  ## write to file
  if (RemoveOutliers==TRUE){Out = ", No Outliers"} else {Out = ""}
  if (centile==0.50){Cent = "Median"} else {Cent = paste0(centile*100, "th Centile")}
  write.xlsx(Results_tosave, file=paste0(results_file, ".xlsx"), sheetName=paste0("Qreg ", Cent, Out), append=TRUE)
  
}

# report the proportion of outliers
Proportion_outliers <- function(MADprop, results_file){
  
  ### load and prepare data
  data_touse <- read.csv("//ads.bris.ac.uk/filestore/BRMS/Research/Metabolic profiles in UPBEAT/Data/upbeat_withmet_data_longSTATA_v2.csv", na.strings='.')
  
  # name the metabolites
  varlist <- colnames(data_touse[11:168])
  
  prop_outliers <- data.frame(matrix(NA, nrow=length(varlist), ncol=6))
  colnames(prop_outliers) <- c("Proportion of values >3.5*MAD", 
                               "Proportion outliers at first timepoint", 
                               "Proportion of outliers in untreated at second timepoint", "Proportion of outliers in treated at second timepoint", 
                               "Proportion of outliers in untreated at third timepoint", "Proportion of outliers in treated at third timepoint")
  rownames(prop_outliers) <- varlist
  for (i in 1:length(varlist)){
    ### remove outliers
    med_var <- vector(len=nrow(data_touse))
    
    # time 15
    med_var[which(data_touse$time_ind==15)] <- median(data_touse[which(data_touse$time_ind==15), varlist[i]], na.rm=TRUE) 
    devvar <- abs(data_touse[, varlist[i]] - med_var)[which(data_touse$time_ind==15)]
    data_touse[which(data_touse$time_ind==15), paste0(varlist[i], "mad")] <- median(devvar, na.rm=TRUE)
    
    # time 28 and treatment 0
    med_var[which(data_touse$time_ind==28 & data_touse$Treatment==0)] <- median(data_touse[which(data_touse$time_ind==28 & data_touse$Treatment==0), varlist[i]], na.rm=TRUE) 
    devvar <- abs(data_touse[, varlist[i]] - med_var)[which(data_touse$time_ind==28 & data_touse$Treatment==0)]
    data_touse[which(data_touse$time_ind==28 & data_touse$Treatment==0), paste0(varlist[i], "mad")] <- median(devvar, na.rm=TRUE)
    
    # time 28 and treatment 1
    med_var[which(data_touse$time_ind==28 & data_touse$Treatment==1)] <- median(data_touse[which(data_touse$time_ind==28 & data_touse$Treatment==1), varlist[i]], na.rm=TRUE) 
    devvar <- abs(data_touse[, varlist[i]] - med_var)[which(data_touse$time_ind==28 & data_touse$Treatment==1)]
    data_touse[which(data_touse$time_ind==28 & data_touse$Treatment==1), paste0(varlist[i], "mad")] <- median(devvar, na.rm=TRUE)
    
    # time 36 and treatment 0
    med_var[which(data_touse$time_ind==36 & data_touse$Treatment==0)] <- median(data_touse[which(data_touse$time_ind==36 & data_touse$Treatment==0), varlist[i]], na.rm=TRUE) 
    devvar <- abs(data_touse[, varlist[i]] - med_var)[which(data_touse$time_ind==36 & data_touse$Treatment==0)]
    data_touse[which(data_touse$time_ind==36 & data_touse$Treatment==0), paste0(varlist[i], "mad")] <- median(devvar, na.rm=TRUE)
    
    # time 36 and treatment 1
    med_var[which(data_touse$time_ind==36 & data_touse$Treatment==1)] <- median(data_touse[which(data_touse$time_ind==36 & data_touse$Treatment==1), varlist[i]], na.rm=TRUE) 
    devvar <- abs(data_touse[, varlist[i]] - med_var)[which(data_touse$time_ind==36 & data_touse$Treatment==1)]
    data_touse[which(data_touse$time_ind==36 & data_touse$Treatment==1), paste0(varlist[i], "mad")] <- median(devvar, na.rm=TRUE)
    
    # indicate outliers
    data_touse[ , paste0(varlist[i], "outlier")] <- 0
    data_touse[which(data_touse[ , varlist[i]] >= (med_var + MADprop*data_touse[, paste0(varlist[i], "mad")])), paste0(varlist[i], "outlier")] <- 1
    data_touse[which(data_touse[ , varlist[i]] <= (med_var - MADprop*data_touse[, paste0(varlist[i], "mad")])), paste0(varlist[i], "outlier")] <- 1
    data_touse[which(is.na(data_touse[ , varlist[i]])) , paste0(varlist[i], "outlier")] <- NA
    
    ### save proportion outliers
    tab_tmp <- table(data_touse[ , paste0(varlist[i], "outlier")], data_touse$time_ind, useNA = "ifany")
    
    # without NA included in totals
    tab_tmp_total <- tab_tmp["0", ] + tab_tmp["1", ] # WITHOUT NA
    timepoint_treatment <- table(data_touse$time_ind[which(!is.na(data_touse[ , varlist[i]]))], data_touse$Treatment[which(!is.na(data_touse[ , varlist[i]]))], useNA="ifany") # WITHOUT NA
    
    # with NA included in totals
    #tab_tmp_total <- colSums(tab_tmp) # WITH NA
    #timepoint_treatment <- table(data_touse$time_ind, data_touse$Treatment) # WITH NA
    
    prop_outliers[varlist[i], paste0("Proportion of values >", MADprop, "*MAD")] <- sum(tab_tmp["1", ]) / sum(tab_tmp_total) 
    prop_outliers[varlist[i], "Proportion outliers at first timepoint"] <- sum(tab_tmp["1", "15"]) / tab_tmp_total["15"]
    prop_outliers[varlist[i], "Proportion of outliers in untreated at second timepoint"] <-  
      table(data_touse[which(data_touse$time_ind==28 & data_touse$Treatment==0), paste0(varlist[i], "outlier")])["1"] / timepoint_treatment["28", "0"]
    prop_outliers[varlist[i], "Proportion of outliers in treated at second timepoint"] <- 
      table(data_touse[which(data_touse$time_ind==28 & data_touse$Treatment==1), paste0(varlist[i], "outlier")])["1"] / timepoint_treatment["28", "1"]
    prop_outliers[varlist[i], "Proportion of outliers in untreated at third timepoint"] <- 
      table(data_touse[which(data_touse$time_ind==36 & data_touse$Treatment==0), paste0(varlist[i], "outlier")])["1"] / timepoint_treatment["36", "0"]
    prop_outliers[varlist[i], "Proportion of outliers in treated at third timepoint"] <- 
      table(data_touse[which(data_touse$time_ind==36 & data_touse$Treatment==1), paste0(varlist[i], "outlier")])["1"] / timepoint_treatment["36", "1"]
    # remove NA
    prop_outliers[varlist[i], which(is.na(prop_outliers[varlist[i], ]))]<-0
  } 
  
  
  ## write to file
  write.xlsx(prop_outliers, file=paste0(results_file, ".xlsx"), sheetName=paste0("Proportion values > ", MADprop, " MAD "), append=TRUE)
}
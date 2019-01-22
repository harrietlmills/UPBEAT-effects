##### code to run a lme model (lmer) and output results for every metabolite

### clear workspace
rm(list = ls())

### load packages
library(lmerTest)
library(arm)

### define path to main folder (contains subfolders called "RScripts" and "Data")
mainpath <- "MainFolder/"

### define path to Results folder (contains subfolder called "Raw")
resultspath <- paste0(mainpath, "/Results/Raw")

### load cleaned data (long rather than wide version)
data_long <- readRDS(file=paste0(mainpath, "Data/data_long.rds"))

### ENSURE CORRECT UNITS OF TIME - HERE WE USE 4 WEEK UNITS
data_long[, c("time_point")] = data_long[, c("time_point")] / 4 #timeunit_4weeks = TRUE

### name the metabolites by extracting the relevant column names
METABOLITES = colnames(data_long)[which(colnames(data_long)=="xxlvldlp"):dim(data_long)[2]]

### rescale those which are very small to be around 0
data_long_orig = data_long
limit = 0.5
nonNA = colSums(!is.na(data_long[, METABOLITES]))
tmp = sapply(METABOLITES, function(e) sum(data_long[, e]<limit, na.rm=TRUE))
Scaled_Met = METABOLITES[which(tmp>(nonNA[METABOLITES]/4))] # if more than 25% are small
Scaled_by = rep(1, length(Scaled_Met)); names(Scaled_by) = Scaled_Met
for (s in 1:length(Scaled_Met)){
  new_col = data_long[, Scaled_Met[s]]
  while (sum(new_col<limit, na.rm=TRUE)>(nonNA[Scaled_Met[s]]/4)){
    new_col = 10 * new_col
    Scaled_by[s] = 10 * Scaled_by[s]
  }
  data_long[, Scaled_Met[s]] = new_col
}
Scaled_by_all = rep(1, length(METABOLITES))
names(Scaled_by_all) = METABOLITES
Scaled_by_all[Scaled_Met] = Scaled_by

## for the one metabolite that didn't converge when scaled, we unscale it
data_long[, "idlce"] = data_long[, "idlce"] / Scaled_by_all["idlce"]
Scaled_by_all["idlce"] = 1

############ fit models to all the data
### create data.frame to save the fit results in and define the columns we want
Fit_Results_lme4 = data.frame(matrix(NA, nrow = length(METABOLITES), ncol = 12))
rownames(Fit_Results_lme4) = METABOLITES
colnames(Fit_Results_lme4) = c("Intercept", "Int 2.5CI", "Int 97.5CI", "Slope", "Slope 2.5CI", "Slope 97.5CI",
                               "TI*time", "TI*time 2.5CI", "TI*time 97.5CI", 
                               "pvalue", "Var_residuals", "Var_intercept")
ERRORS_lme4 = CI_error = vector(len=length(METABOLITES))
parm_CI = c("(Intercept)", "time_point", "time_point:TreatmentIntervention") 

### create data.frame to save the predictions in, with confidence intervals and standardised versions
summary_data_predicted = data.frame(matrix(NA, nrow=length(METABOLITES), ncol=18))
colnames(summary_data_predicted) = c("Mean_16wk", "Median_16wk", "SD_16wk", "Mean_36wk", "Median_36wk", "SD_36wk", 
                                     "Mean_diff", "SE_Mean_diff", "Lower_95CI_diff", "Upper_95CI_diff", 
                                     "Mean_diff_16weekSD_units", "SE_Mean_diff_16weekSD_units", "Lower_95CI_diff_16weekSD_units", "Upper_95CI_diff_16weekSD_units",
                                     "Mean_diff_SD_units", "SE_Mean_diff_SD_units", "Lower_95CI_diff_SD_units", "Upper_95CI_diff_SD_units")
rownames(summary_data_predicted) = METABOLITES

for (m in 1:length(METABOLITES)){
  metabolite = METABOLITES[m]
  print(metabolite)
  # in a tryCatch to catch any warning / error messages from models
  lmer_mod <- tryCatch(lmer(as.formula(paste0(metabolite, " ~ BMI + Ethnicity + Parity + Age + Centre + time_point + Treatment:time_point + (time_point|participant_id)")), data=data_long), 
                            error = function(e) {e$message} )

  ## records results of model, scaling up when we save them if needed
  if (class(lmer_mod)!="character"){ 
    #co_lmer = coef(lmer_mod)[[1]] # mean(co_lmer$"(Intercept)") & mean(co_lmer$"time_point") are the intercept and slope below
    su_lmer = summary(lmer_mod)
    model_intercept = su_lmer$coefficients["(Intercept)", "Estimate"] 
    model_slope = su_lmer$coefficients["time_point", "Estimate"]
    
    # confidence intervals http://stats.stackexchange.com/questions/117641/how-trustworthy-are-the-confidence-intervals-for-lmer-objects-through-effects-pa
    CI_mod = tryCatch(confint(lmer_mod, parm = parm_CI, method="Wald"), error = function(e) {e$message} ) 
    if (is.character(CI_mod)){ 
      CI_error[m] = CI_mod; 
      CI_mod = matrix(NA, ncol=2, nrow=length(parm_CI)); 
      rownames(CI_mod) = parm_CI
    }
    
    # pvalues
    anova_mod = anova(lmer_mod) # compares main model to null model for each fixed effect
    pvalue_interaction = anova_mod["time_point:Treatment", "Pr(>F)"]
    if (is.null(pvalue_interaction)){pvalue_interaction = NA}

    ## save fit results
    # control_intercept is the intercept
    # control_slope is the time_point value
    # difference in intervention and control group is the TreatmentIntBYtime_point value
    Fit_Results_lme4[m, ] <- c(Intercept = model_intercept/Scaled_by_all[metabolite], # intercept
                              CI_mod["(Intercept)", ]/Scaled_by_all[metabolite], #include CI
                              Slope = model_slope/Scaled_by_all[metabolite], # slope
                              CI_mod["time_point", ]/Scaled_by_all[metabolite], #include CI
                              TreatmentIntBYtime_point = su_lmer$coefficients["time_point:TreatmentIntervention", "Estimate"]/Scaled_by_all[metabolite], # estimate coefficient of Treatment*time_point + 
                              CI_mod["time_point:TreatmentIntervention", ]/Scaled_by_all[metabolite], #include CI
                              pvalue_interaction, # pvalue
                              Var_residuals = attr(VarCorr(lmer_mod),"sc")^2, # variance of residuals
                              Var_intercept = 1/dim(se.ranef(lmer_mod)$participant_id)[1] * sum((se.ranef(lmer_mod)$participant_id[, "(Intercept)"])^2) + var(ranef(lmer_mod)$participant_id[, "(Intercept)"])) # variance of the intercept???)
    
    #### summarise the change in metabolite from 16 to 36 weeks in the control arm
    # derive predicted level at 16 an 36 weeks for each outcome
    data_long_new = data_long[data_long$Treatment=="Control",] # limit to the control dataset
    #if (timeunit_4weeks){ # redefine the time points correctly spaced, making every woman the same (0, 16 and 36 weeks)
      data_long_new$time_point = rep(c(0, 12, 20)/4)
    # } else{
    #   data_long_new$time_point = rep(c(0, 12, 20)) # recall the time points are zero'd by removing 16
    # }
    p = predict(lmer_mod, newdata = data_long_new, allow.new.levels = T)
    data_predicted = rep(NA, times=dim(data_long)[1])
    data_predicted[as.numeric(names(p))] = p
    data_16 = data_predicted[seq(1, dim(data_long)[1], 3)]/Scaled_by_all[metabolite]
    data_36 = data_predicted[seq(3, dim(data_long)[1], 3)]/Scaled_by_all[metabolite]

    # summarise these predicted levels at each time point 
    summary_data_predicted[metabolite, c("Mean_16wk", "Median_16wk", "SD_16wk")] = c(mean(data_16, na.rm=TRUE), median(data_16, na.rm=TRUE), sd(data_16, na.rm=TRUE))
    summary_data_predicted[metabolite, c("Mean_36wk", "Median_36wk", "SD_36wk")] = c(mean(data_36, na.rm=TRUE), median(data_36, na.rm=TRUE), sd(data_36, na.rm=TRUE))

    # subtract the 16 from the 36 values and derive the mean (standard error of that mean) and 95%CI for it so we have the total mean difference between 16 and 36 weeks in control women 
    change =  data_36 - data_16
    MEAN_DIFF = mean(change, na.rm=TRUE)
    SE = sd(change, na.rm=TRUE) / sqrt(length(change[!is.na(change)]))
    summary_data_predicted[metabolite, c("Mean_diff", "SE_Mean_diff", "Lower_95CI_diff", "Upper_95CI_diff")] = c(MEAN_DIFF, SE, MEAN_DIFF - 1.96*SE, MEAN_DIFF + 1.96*SE)
    
    # mean diff in SD units of data_16
    MEAN_DIFF_16WeekSD_Units = mean(change / summary_data_predicted[metabolite, c("SD_16wk")], na.rm=TRUE) # same as summary_data_predicted[metabolite, c("Mean_diff")] / summary_data_predicted[metabolite, c("SD_16wk")] # 
    SE = sd(change / summary_data_predicted[metabolite, c("SD_16wk")], na.rm=TRUE) / sqrt(length(change[!is.na(change)]))
    summary_data_predicted[metabolite, c("Mean_diff_16weekSD_units", "SE_Mean_diff_16weekSD_units", "Lower_95CI_diff_16weekSD_units", "Upper_95CI_diff_16weekSD_units")] = 
      c(MEAN_DIFF_16WeekSD_Units, SE, MEAN_DIFF_16WeekSD_Units - 1.96*SE, MEAN_DIFF_16WeekSD_Units + 1.96*SE)

    # mean diff in SD units of change
    MEAN_DIFF_SD_Units = mean(change, na.rm=TRUE)/ sd(change, na.rm=TRUE) 
    SE = 1 / sqrt(length(change[!is.na(change)])) # SD cancels out
    summary_data_predicted[metabolite, c("Mean_diff_SD_units", "SE_Mean_diff_SD_units", "Lower_95CI_diff_SD_units", "Upper_95CI_diff_SD_units")] = 
      c(MEAN_DIFF_SD_Units, SE, MEAN_DIFF_SD_Units - 1.96*SE, MEAN_DIFF_SD_Units + 1.96*SE)
    
    
  } else{
    ERRORS_lme4[m] <- lmer_mod$error
  }
}

save(Fit_Results_lme4, ERRORS_lme4, CI_error, file=paste0(resultspath, "LME_Model_FullResults.Rdata"))
write.csv(summary_data_predicted, file=paste0(resultspath, "LME_Model_predicted_mean_change.csv"), row.names=TRUE)


          
          
##### code to run a lme model (lmer) and output results for every metabolite (where metabolites have been standardised)

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
data_long[, c("time_point")] = data_long[, c("time_point")] / 4

### name the metabolites by extracting the relevant column names
METABOLITES = colnames(data_long)[which(colnames(data_long)=="xxlvldlp"):dim(data_long)[2]]

### standardise at each time point to the control group
time_indicators = c(15, 28, 36)
w_control = which(data_long$Treatment%in%"Control")
for (tp in 1:3){
  time_totake = which(data_long$time_ind==time_indicators[tp])
  controlrows_totake = intersect(w_control, time_totake)
  for (m in 1:length(METABOLITES)){
    # find mean and sd just for the control at that time
    mean_control = mean(data_long[controlrows_totake, METABOLITES[m]], na.rm=TRUE)
    sd_control = sd(data_long[controlrows_totake, METABOLITES[m]], na.rm=TRUE)
    # adjust all at that time
    data_long[time_totake, METABOLITES[m]] = (data_long[time_totake, METABOLITES[m]] - mean_control) / sd_control
  }
}

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
summary_data_predicted = data.frame(matrix(NA, nrow=length(METABOLITES), ncol=14))
colnames(summary_data_predicted) = c("Mean_16wk", "Median_16wk", "SD_16wk", "Mean_36wk", "Median_36wk", "SD_36wk", 
                                     "Mean_diff", "SE_Mean_diff", "Lower_95CI_diff", "Upper_95CI_diff", 
                                     "Mean_diff_16weekSD_units", "SE_Mean_diff_16weekSD_units", "Lower_95CI_diff_16weekSD_units", "Upper_95CI_diff_16weekSD_units")
rownames(summary_data_predicted) = METABOLITES

for (m in 1:length(METABOLITES)){
  metabolite = METABOLITES[m]
  print(metabolite)
  lmer_mod <- tryCatch(lmer(as.formula(paste0(metabolite, " ~ BMI + Ethnicity + Parity + Age + Centre + time_point + Treatment:time_point + (time_point|participant_id)")), data=data_long), 
                            error = function(e) {e$message} )

  ## records results of model, scaling up when we save them if needed
  if (class(lmer_mod)!="character"){ 
    #co_lmer = coef(lmer_mod)[[1]]
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
    Fit_Results_lme4[m, ] <- c(Intercept = model_intercept, # intercept
                              CI_mod["(Intercept)", ], #include CI
                              Slope = model_slope, # slope
                              CI_mod["time_point", ], #include CI
                              TreatmentIntBYtime_point = su_lmer$coefficients["time_point:TreatmentIntervention", "Estimate"], # estimate coefficient of Treatment*time_point + 
                              CI_mod["time_point:TreatmentIntervention", ], #include CI
                              pvalue_interaction, # pvalue
                              Var_residuals = attr(VarCorr(lmer_mod),"sc")^2, # variance of residuals
                              Var_intercept = 1/dim(se.ranef(lmer_mod)$participant_id)[1] * sum((se.ranef(lmer_mod)$participant_id[, "(Intercept)"])^2) + var(ranef(lmer_mod)$participant_id[, "(Intercept)"])) # variance of the intercept???)

  } else{
    ERRORS_lme4[m] <- lmer_mod$error
  }
}

save(Fit_Results_lme4, ERRORS_lme4, CI_error, file=paste0(resultspath, "/LME_Model_Standardised_FullResults.Rdata"))


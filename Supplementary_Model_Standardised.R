##### code to run a lme model (lmer) and output results for every metabolite (where metabolites have been standardised)
## uses a spline with single knot point

### clear workspace
rm(list = ls())

### load packages
library(lmerTest)
library(arm)
library(splines)

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
Fit_Results_lme4 = data.frame(matrix(NA, nrow = length(METABOLITES), ncol = 19))
rownames(Fit_Results_lme4) = METABOLITES
colnames(Fit_Results_lme4) = c("Intercept", "Int 2.5CI", "Int 97.5CI", "Slope1", "Slope1 2.5CI", "Slope1 97.5CI", "Slope2", "Slope2 2.5CI", "Slope2 97.5CI",
                               "TI*gest1", "TI*gest1 2.5CI", "TI*gest1 97.5CI", "pvalue1", "TI*gest2", "TI*gest2 2.5CI", "TI*gest2 97.5CI", "pvalue2",  
                               "Var_residuals", "Var_intercept")
ERRORS_lme4 = CI_error = vector(len=length(METABOLITES))
parm_CI = c("(Intercept)", "gesttime1", "gesttime2", "gesttime1:TreatmentIntervention", "gesttime2:TreatmentIntervention") 

## spline

gesttime1 = rep(c(-1, 12, 12), times=1158)
gesttime2 = rep(c(0, 0, 8 ), times=1158)
#if (timeunit_4weeks){
  gesttime1 = gesttime1/4
  gesttime2 = gesttime2/4
#}


for (m in 1:length(METABOLITES)){
  metabolite = METABOLITES[m]
  print(metabolite)
  
  ## model
  lmer_mod <- tryCatch(lmer(as.formula(paste0(metabolite, " ~ BMI + Ethnicity + Parity + Age + Centre + gesttime1 + gesttime2 + Treatment:gesttime1 + Treatment:gesttime2 + (time_point|participant_id)")), data=data_long), 
                       error = function(e) {e$message} )
  
  ## records results of model, scaling up when we save them if needed
  if (class(lmer_mod)!="character"){ 
    #co_lmer = coef(lmer_mod)[[1]]
    su_lmer = summary(lmer_mod)
    model_intercept = su_lmer$coefficients["(Intercept)", "Estimate"] 
    model_slope1 = su_lmer$coefficients["gesttime1", "Estimate"]
    model_slope2 = su_lmer$coefficients["gesttime2", "Estimate"]
    model_interactionslope1 = su_lmer$coefficients["gesttime1:TreatmentIntervention", "Estimate"] 
    model_interactionslope2 = su_lmer$coefficients["gesttime2:TreatmentIntervention", "Estimate"] 
    
    # confidence intervals http://stats.stackexchange.com/questions/117641/how-trustworthy-are-the-confidence-intervals-for-lmer-objects-through-effects-pa
    CI_mod = tryCatch(confint(lmer_mod, parm = parm_CI, method="Wald"), error = function(e) {e$message} ) 
    if (is.character(CI_mod)){ # suggests maybe you could change dev.tol but didn't seem to work when I tried (more warnings) https://marc.info/?l=r-sig-mixed-models&m=140605544029403&w=4
      CI_error[m] = CI_mod; 
      CI_mod = matrix(NA, ncol=2, nrow=length(parm_CI)); 
      rownames(CI_mod) = parm_CI
    }

    # pvalues
    anova_mod = anova(lmer_mod) # compares main model to null model for each fixed effect
    pvalue_interaction_gesttime1 = anova_mod["gesttime1:Treatment", "Pr(>F)"]
    pvalue_interaction_gesttime2 = anova_mod["gesttime2:Treatment", "Pr(>F)"]
    if (is.null(pvalue_interaction_gesttime1)){pvalue_interaction_gesttime1 = NA}
    if (is.null(pvalue_interaction_gesttime2)){pvalue_interaction_gesttime2 = NA}
    
    ## save fit results
    Fit_Results_lme4[m, ] <- c(Intercept = model_intercept, # intercept
                               CI_mod["(Intercept)", ], #include CI
                               Slope1 = model_slope1, # slope
                               CI_mod["gesttime1", ], #include CI
                               Slope2 = model_slope2, # slope
                               CI_mod["gesttime2", ], #include CI
                               TreatmentIntBYgesttime1 = model_interactionslope1, 
                               CI_mod["gesttime1:TreatmentIntervention", ], #include CI
                               pvalue1 = pvalue_interaction_gesttime1, 
                               TreatmentIntBYgesttime2 = model_interactionslope2, 
                               CI_mod["gesttime2:TreatmentIntervention", ], #include CI
                               pvalue2 = pvalue_interaction_gesttime2, 
                               Var_residuals = attr(VarCorr(lmer_mod),"sc")^2, # variance of residuals
                               Var_intercept = 1/dim(se.ranef(lmer_mod)$participant_id)[1] * sum((se.ranef(lmer_mod)$participant_id[, "(Intercept)"])^2) + var(ranef(lmer_mod)$participant_id[, "(Intercept)"])) # variance of the intercept???)

  } else{
    ERRORS_lme4[m] <- lmer_mod$error
  }
}

save(Fit_Results_lme4, ERRORS_lme4, CI_error, file=paste0(resultspath, "/LME_SplineModel_Standardised_FullResults.Rdata"))


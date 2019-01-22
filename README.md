# UPBEAT-effects
R code for "The effect of a lifestyle intervention in obese pregnant women on change in gestational metabolic profiles: findings from the UK Pregnancies Better Eating and Activity Trial (UPBEAT) RCT"

Code for the model and analysis is described in the supplement for this paper on BioRxiv https://doi.org/10.1101/125740, and reproduced below. The code for the sensitivity analysis was added after review, and is relevant for the final published version of the paper.

## Note about the folder system the code expects
- "MainFolder
    - "RScripts"
    - "Data"
    - "Results"
        - "Raw"

Note that the data files used in the R code (data_long.rds and metabolome_data_dictionary.csv) are not provided. Researchers interested in accessing the UPBEAT data can do so via a data request form available from the study website http://www.medscinet.net/upbeat/default.aspx)

## Model scripts
- *Main_Model.R* Takes the data in long format, scales metabolite measures which are very small to improve the fits, fits the multilevel model to every metabolite and records the results and the predictions.
- *Main_Model_Standardised.R* As above, but standardises the metabolite measures first.
- *Supplementary_Model_Standardised.R* Takes the data in long format, standardises the metabolite measures, scales metabolite measures which are very small to improve the fits, fits the supplementary model (with knot point at 28 weeks) to every metabolite and records the results.

## Analysis scripts
- *UsefulFunctions_UPBEAT.R* Two useful functions written for analysis methods.
- *Analysis_Code.R* Analyse results from both models and produce tables and figures for the main text and supplement. 

## Sensitivity analysis
- *Functions_forSA.R* – script containing all the functions written for the Sensitivity Analyses. Separate functions for the MLM, bootstrap and Qreg methods, and for giving the proportion of outliers. In these functions the data is loaded, outliers are removed if needed, data is scaled accordingly before the analysis is run. Results are saved to the excel file directly from these functions.
- *Run_SA.R* - main R script for running the Sensitivity Analyses. Creates an excel file called “SensAnalysis_Month_Day_Year”, eg SensAnalysis_Jan_15_2019. It then runs each analysis and writes the results tables to separate sheets in the excel file, each sheet is clearly named for that analysis. The results tables have one row for every metabolite and columns including a p-value.



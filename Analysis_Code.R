################## analysis code to produce the tables and plots in the paper
## note that Table 1 and sTable1 are demographic data, and are not produced here
## note that Figure 1 is a flowchart, not produced here
## note that STable 2 is a list of metabolites and is not produced here
## note that sFigure 1 describes methods and stages of NMR spectroscropy platform and is not produced here

### clear workspace
rm(list = ls())

### load packages
library(scales)

### define path to main folder (contains subfolders called "RScripts", "Results" (with subfolder "Raw") and "Data")
#mainpath <- "MainFolder/"
mainpath <- "//ads.bris.ac.uk/filestore/BRMS/Research/Metabolic profiles in UPBEAT/"

### source useful functions
source(paste0(mainpath, "RScripts/UsefulFunctions_UPBEAT.R"))

########################################################################
## Tables and Figures for Main Model
########################################################################

### load results for unstandardised model
load(file=paste0(mainpath, "Paper_TablesFigures/Raw/LME_Model_FullResults.Rdata"))
# Fit_Results_lme4, ERRORS_lme4, CI_error, 

### sort out the groupings
dictionary <- read.csv(paste0(mainpath, "Data/metabolome_data_dictionary.csv"), row.names=1)
# unfactor dictionary columns
for (i in 1:ncol(dictionary)){
  dictionary[, i] = as.character(dictionary[, i])
}

# trim white space from dictionary entries
dictionary$Group <- trim(dictionary$Group)
dictionary$Biomarker <- trim(dictionary$Biomarker)

# assign Group and Biomarker to Fit_results
Fit_Results_lme4$Group = dictionary[rownames(Fit_Results_lme4), "Group"]
Fit_Results_lme4$Biomarker = dictionary[rownames(Fit_Results_lme4), "Biomarker"]

# add a gap before each new group and plot nothing, but add group name
groups = unique(dictionary$Group)
for (g in 1:length(groups)){
  tmp = rbind(rep(NA, dim(Fit_Results_lme4)[2]), Fit_Results_lme4[Fit_Results_lme4$Group==groups[g], ])
  rownames(tmp) = c(groups[g], Fit_Results_lme4$Biomarker[Fit_Results_lme4$Group==groups[g]])
  if (g==1){
    Fit_Results_withsections = tmp
  } else{
    Fit_Results_withsections = rbind(Fit_Results_withsections, tmp)
  }
}

########################
# sTable 3
########################

summary_data_predicted = read.csv(file=paste0(mainpath, "Paper_TablesFigures/Raw/LME_Model_predicted_mean_change.csv"), row.names=1)

summary_data_predicted_Group = dictionary[rownames(summary_data_predicted), "Group"]
summary_data_predicted_Biomarker = dictionary[rownames(summary_data_predicted), "Biomarker"]
summary_data_predicted$VarName = rownames(summary_data_predicted)
rownames(summary_data_predicted) = summary_data_predicted$Biomarker = summary_data_predicted_Biomarker # change rownames to helpful thing
summary_data_predicted$Group = summary_data_predicted_Group

# add a gap before each new group and add group name
groups = unique(dictionary$Group)
for (g in 1:length(groups)){
  tmp = rbind(rep(NA, dim(summary_data_predicted)[2]), summary_data_predicted[summary_data_predicted$Group==groups[g], ])
  rownames(tmp) = c(groups[g], summary_data_predicted$Biomarker[summary_data_predicted$Group==groups[g]])
  if (g==1){
    summary_data_predicted_withsections = tmp
  } else{
    summary_data_predicted_withsections = rbind(summary_data_predicted_withsections, tmp)
  }
}
NICE_summary_data_predicted = cbind(paste0(format_nicely(summary_data_predicted_withsections[, "Mean_diff"], 4), " (", format_nicely(summary_data_predicted_withsections[, "Lower_95CI_diff"], 4), ", ", format_nicely(summary_data_predicted_withsections[, "Upper_95CI_diff"], 4), ")"),
                                    paste0(format_nicely(summary_data_predicted_withsections[, "Mean_diff_16weekSD_units"], 4), " (", format_nicely(summary_data_predicted_withsections[, "Lower_95CI_diff_16weekSD_units"], 4), ", ", format_nicely(summary_data_predicted_withsections[, "Upper_95CI_diff_16weekSD_units"], 4), ")"))
rownames(NICE_summary_data_predicted) = rownames(summary_data_predicted_withsections) 
colnames(NICE_summary_data_predicted) = c("Mean absolute difference between 16 and 36 weeks of gestation in original units (95% CI)",
                                          "Mean absolute difference between 16 and 36 weeks of gestational age in SD units")  #c("Mean diff (95CI)", "Mean SD diff (95CI)")

NICE_summary_data_predicted[grep("NA", NICE_summary_data_predicted[, 1]), 1:2] = ""

write.csv(NICE_summary_data_predicted, file=paste0(mainpath, "Paper_TablesFigures/sTable3.csv"), row.names=TRUE)


########################
# sTable 4
########################

sTable4 = cbind(format_nicely(Fit_Results_withsections[, "Intercept"], 3) ,
                paste0("(", format_nicely(Fit_Results_withsections[, "Int 2.5CI"], 3), ", ", format_nicely(Fit_Results_withsections[, "Int 97.5CI"], 3), ")"),
                format_nicely(Fit_Results_withsections[, "Slope"], 3) ,
                paste0("(", format_nicely(Fit_Results_withsections[, "Slope 2.5CI"], 3), ", ", format_nicely(Fit_Results_withsections[, "Slope 97.5CI"], 3), ")"))
sTable4[grep("NA", sTable4[, 1]), ] = rep("", 4)
sTable4 = cbind(paste(sTable4[, 1], sTable4[, 2]), paste(sTable4[, 3], sTable4[, 4])) 
colnames(sTable4) = c("Mean concentration at 16-weeks (95%CI)", "Mean change in concentration per 4 weeks gestational age (95%CI)") #c("Intercept (95% CI)", "Slope (95% CI)")
rownames(sTable4) = rownames(Fit_Results_withsections) 
write.csv(sTable4, file=paste0(mainpath, "Paper_TablesFigures/sTable4.csv"))


########################
# sTable 5
########################

sTable5 = cbind(format_nicely(Fit_Results_withsections[, "TI*time"], 3) ,
                paste0("(", format_nicely(Fit_Results_withsections[, "TI*time 2.5CI"], 3), ", ", format_nicely(Fit_Results_withsections[, "TI*time 97.5CI"], 3), ")"),
                format_nicely(Fit_Results_withsections[, "pvalue"], 3))
sTable5[grep("NA", sTable5[, 1]), ] = rep("", 3)
sTable5 = cbind(paste(sTable5[, 1], sTable5[, 2]), sTable5[, 3]) 
colnames(sTable5) = c("Difference in mean rate of change in traits per 4 weeks of gestation between 16 and 36 weeks between women receiving intervention and control group (reference)",
                      "p-value") # c("Diff Intervention Slope (95% CI)", "p-value")
rownames(sTable5) = rownames(Fit_Results_withsections) 
write.csv(sTable5, file=paste0(mainpath, "Paper_TablesFigures/sTable5.csv"), row.names=TRUE)


########################
# Table 2 (a subset of sTable5) - FIRST COLUMN GENERATED HERE, ADDED TO LATER
########################

subset_groups = c("Extremely large VLDL", "Very large VLDL","Large VLDL", "Medium VLDL")
subset_otherLipids = c("Triglycerides in very large HDL", "Phospholipids in small HDL" ,"Mean diameter for VLDL particles","Serum total triglycerides",
"Triglycerides in VLDL" ,"Ratio of triglycerides to phosphoglycerides" ,"Estimated degree of unsaturation","Ratio of 18:2 linoleic acid to total fatty acids",
"Ratio of omega-6 fatty acids to total fatty acids","Ratio of polyunsaturated fatty acids to total fatty acids" ,"Ratio of saturated fatty acids to total fatty acids")
subset_otherTraits = c("Lactate", "Pyruvate", "Alanine" , "Acetate")

for (s in 1:length(subset_groups)){
  w = which(Fit_Results_withsections$Group%in%subset_groups[s])
  if (s==1){
    Table_Groups = Fit_Results_withsections[c(w[1]-1, w), ]
  } else{
    Table_Groups = rbind(Table_Groups, Fit_Results_withsections[c(w[1]-1, w), ])
  }
}

Table_otherLipids = Fit_Results_withsections[which(Fit_Results_withsections$Biomarker%in%subset_otherLipids), ]
Table_otherTraits = Fit_Results_withsections[which(Fit_Results_withsections$Biomarker%in%subset_otherTraits), ]

Whole_Table_unstandardised = rbind(Table_Groups,
                    rep(NA, dim(Fit_Results_lme4)[2]), Table_otherLipids,
                    rep(NA, dim(Fit_Results_lme4)[2]), Table_otherTraits)
rownames(Whole_Table_unstandardised)[33] = "Other lipids, lipoproteins affected by the intervention"
rownames(Whole_Table_unstandardised)[45] = "Other traits affected by the intervention"

Table2 = cbind(format_nicely(Whole_Table_unstandardised[, "TI*time"], 3) ,
      paste0("(", format_nicely(Whole_Table_unstandardised[, "TI*time 2.5CI"], 3), ", ", format_nicely(Whole_Table_unstandardised[, "TI*time 97.5CI"], 3), ")"),
      NA)
Table2[grep("NA", Table2[, 1]), ] = rep("", 3)
Table2 = cbind(paste(Table2[, 1], Table2[, 2]), Table2[, 3])
# Difference in mean rate of change in metabolic traits per 4-weeks of gestation between 16 and 36 weeks comparing women receiving intervention to control group
colnames(Table2) = c("In original units (see first column) per 4-weeks",
                      "In SD units per 4-weeks")
rownames(Table2) = rownames(Whole_Table_unstandardised)

write.csv(Table2, file=paste0(mainpath, "Paper_TablesFigures/Table2_col1.csv"), row.names=TRUE)


########################
# Figure 2
########################
# load simple labels
file=load(paste0(mainpath, "RScripts/metabolite_labels.Rdata"))
metabolite_labels$label.no.units[metabolite_labels$metabolite=="tgpg"] = "Triglycerides/Phosphoglycerides (%)"

## normal units
#diff = "Mean_diff"
#lowerCI = "Lower_95CI_diff"
#upperCI = "Upper_95CI_diff"

## SD units
diff = "Mean_diff_16weekSD_units"# "Mean_diff_SD_units"
lowerCI = "Lower_95CI_diff_16weekSD_units"# "Lower_95CI_diff_SD_units"
upperCI = "Upper_95CI_diff_16weekSD_units"# "Upper_95CI_diff_SD_units"

## load data
summary_data_predicted = read.csv(file=paste0(mainpath, "Paper_TablesFigures/Raw/LME_Model_predicted_mean_change.csv"), row.names=1)

summary_data_predicted_Group = dictionary[rownames(summary_data_predicted), "Group"]
summary_data_predicted_Biomarker = dictionary[rownames(summary_data_predicted), "Biomarker"]
summary_data_predicted$VarName = rownames(summary_data_predicted)
summary_data_predicted$Group = summary_data_predicted_Group

rownames(summary_data_predicted) = summary_data_predicted$Biomarker = summary_data_predicted_Biomarker # change rownames to helpful thing
rownames(summary_data_predicted)[which(rownames(summary_data_predicted)%in%"glycerol")] = "Glycerol"
#summary_data_predicted$Biomarker[which(rownames(summary_data_predicted)%in%"Glycerol")] = "Glycerol"
summary_data_predicted$label = metabolite_labels[rownames(summary_data_predicted), "label.no.units"]

## create pos and neg subsets
summary_data_predicted_pos = summary_data_predicted[summary_data_predicted[,diff]>0, ]
summary_data_predicted_neg = summary_data_predicted[summary_data_predicted[,diff]<0, ]

# add a gap before each new group and add group name
groups = unique(dictionary$Group)
started = "No"; first_g=1
for (g in 1:length(groups)){
  # for pos
  if (sum(summary_data_predicted_pos$Group==groups[g])>0){ # then that group exists
    started = "Yes"
    tmp = rbind(rep(NA, dim(summary_data_predicted_pos)[2]), summary_data_predicted_pos[summary_data_predicted_pos$Group==groups[g], ])
    rownames(tmp) = c(groups[g], summary_data_predicted_pos$Biomarker[summary_data_predicted_pos$Group==groups[g]])
    tmp$label[1] = groups[g]
    if (g==first_g){
      summary_data_predicted_pos_withsections = tmp
    } else{
      summary_data_predicted_pos_withsections = rbind(summary_data_predicted_pos_withsections, tmp)
    }    
  } else if (started == "No"){
    first_g = g + 1
  }
}
started = "No"; first_g=1
for (g in 1:length(groups)){
  # for neg
  if (sum(summary_data_predicted_neg$Group==groups[g])>0){ # then that group exists
    started = "Yes"
    tmp = rbind(rep(NA, dim(summary_data_predicted_neg)[2]), summary_data_predicted_neg[summary_data_predicted_neg$Group==groups[g], ])
    rownames(tmp) = c(groups[g], summary_data_predicted_neg$Biomarker[summary_data_predicted_neg$Group==groups[g]])
    tmp$label[1] = groups[g]
    if (g==first_g){
      summary_data_predicted_neg_withsections = tmp
    } else{
      summary_data_predicted_neg_withsections = rbind(summary_data_predicted_neg_withsections, tmp)
    }
  } else if (started == "No"){
    first_g = g + 1
  }
}

#rownames(summary_data_predicted_pos_withsections) = gsub("Concentration", "Conc.", rownames(summary_data_predicted_pos_withsections))
#rownames(summary_data_predicted_pos_withsections) = gsub("and extremely", "& ext.", rownames(summary_data_predicted_pos_withsections))
#rownames(summary_data_predicted_pos_withsections) = gsub("total fatty acids", "TFAs", rownames(summary_data_predicted_pos_withsections))

## normal units
#xlim_1.1 = xlim_1.2 = xlim_2.1 = xlim_2.2 = xlim_3.1 = c(-0.05, 0.5)
#space = c(1, 0)

## SD units
xlim_1.1 = xlim_1.2 = c(-0.1, 4.6)
xlim_2.1 = xlim_2.2 = c(-0.1, 4.6)
xlim_3.1 = c(-3, 0.5)
space1 = c(0.4, 0)
space2 = c(3.4, 0)
space3 = c(2, 0)

## plot positive sections
pdf(paste0(mainpath, "Paper_TablesFigures/Figure2a.pdf"), width=7.48031, height=9.84252)
par(mfrow=c(1, 2), mar=c(2.5, 0.5, 0.5, 0.5), oma = c(0, 0, 0, 1))
### column 1
y = 1:40 # dim(summary_data_predicted_pos_withsections)[1]
pch_list = rep(16, length(y)); #pch_list[rev(summary_data_predicted_pos_withsections[y, "pvalue"])<pvalue_sig] = 16; #pch_list[rev(summary_data_predicted_pos_withsections[y, "pvalue"])<pvalue_BFsig] = 17; 
plot(rev(summary_data_predicted_pos_withsections[y, diff]), y, pch=pch_list, col="red", axes=FALSE, ylab="", xlab="", xlim=xlim_1.1 - space1, ylim=c(min(y)-0.5, max(y)+0.5), xaxs="i", yaxs="i", cex=0.5) #xlim=c(round(min_x, 2), round(max_x, 2)))
# add details
gr = unique(summary_data_predicted_pos_withsections$Group[y]); gr=gr[!is.na(gr)]
for (g in 1:length(gr)){
  # add solid line at 0
  #where = y[grep(gr[g], rev(summary_data_predicted_pos_withsections$Group[y]))]; lines(x=c(0,0), y=c(min(where)-0.5, max(where)+0.5))
  # add dashed lines
  where = y[grep(gr[g], rev(summary_data_predicted_pos_withsections$Group[y]))]; lines(x=c(0,0), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2)
  #where = y[grep(gr[g], rev(summary_data_predicted_pos_withsections$Group[y]))]; lines(x=rep(xlim_1.1[1],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2); lines(x=rep(xlim_1.1[2],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2)
  # add grey sections
  where = y[grep(gr[g], rev(summary_data_predicted_pos_withsections$Group[y]))]
  greylines = where[seq(1, length(where), 2)]
  rect(rep((xlim_1.1 - space1)[1], length(greylines)), greylines-0.5, rep((xlim_1.1 - space1)[2], length(greylines)), greylines+0.5, col=alpha("grey", 0.2), border=NA)
}
# add points and CI lines (on top of grey and vertical lines)
for (i in 1:length(y)){ #dim(summary_data_predicted_pos_withsections)[1]){
  lines(summary_data_predicted_pos_withsections[rev(y)[i], c(lowerCI, upperCI)], rep(y[i], 2), col="red", lwd=1)
}
points(rev(summary_data_predicted_pos_withsections[y, diff]), y, pch=pch_list, col="red", cex=0.5) # re-add dots on top of line
#axis(side=1, at=seq(xlim_1.1[1], xlim_1.1[2], 2), padj=-2.5, cex.axis=0.6)
axis(side=1, at=c(xlim_1.1[1], xlim_1.1[2]), tck=0, labels=FALSE)
axis(side=1, padj=-2.5, cex.axis=0.6)
# add names
sizes = rep(0.6, times=length(y)); sizes[is.na(summary_data_predicted_pos_withsections$Biomarker[y])] = 0.8
text(x=xlim_1.1[1]-space1[1], y=y, labels = rev(summary_data_predicted_pos_withsections$label[y]), pos=4, cex=rev(sizes))
# add bottom x-axis
mtext(side=1, at=(xlim_1.1[2] - xlim_1.1[1])/2, "Mean absolute difference between\n16 and 36 weeks in the control group (SD units)", line=1.25, cex=0.5)


### column 2
y = 41:80 # dim(summary_data_predicted_pos_withsections)[1]
pch_list = rep(16, length(y)); #pch_list[rev(summary_data_predicted_pos_withsections[y, "pvalue"])<pvalue_sig] = 16; #pch_list[rev(summary_data_predicted_pos_withsections[y, "pvalue"])<pvalue_BFsig] = 17; 
plot(rev(summary_data_predicted_pos_withsections[y, diff]), y, pch=pch_list, col="red", axes=FALSE, ylab="", xlab="", xlim=xlim_1.2 - space1, ylim=c(min(y)-0.5, max(y)+0.5), xaxs="i", yaxs="i", cex=0.5) #xlim=c(round(min_x, 2), round(max_x, 2)))
# add details
gr = unique(summary_data_predicted_pos_withsections$Group[y]); gr=gr[!is.na(gr)]
for (g in 1:length(gr)){
  # add solid line at 0
  #where = y[grep(gr[g], rev(summary_data_predicted_pos_withsections$Group[y]))]; lines(x=c(0,0), y=c(min(where)-0.5, max(where)+0.5))
  # add dashed lines
  where = y[grep(gr[g], rev(summary_data_predicted_pos_withsections$Group[y]))]; lines(x=c(0,0), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2)
  #where = y[grep(gr[g], rev(summary_data_predicted_pos_withsections$Group[y]))]; lines(x=rep(xlim_1.2[1],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2); lines(x=rep(xlim_1.2[2],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2)
  # add grey sections
  where = y[grep(gr[g], rev(summary_data_predicted_pos_withsections$Group[y]))]
  greylines = where[seq(1, length(where), 2)]
  rect(rep((xlim_1.2 - space1)[1], length(greylines)), greylines-0.5, rep((xlim_1.2 - space1)[2], length(greylines)), greylines+0.5, col=alpha("grey", 0.2), border=NA)
}
# add points and CI lines (on top of grey and vertical lines)
for (i in 1:length(y)){ #dim(summary_data_predicted_pos_withsections)[1]){
  lines(summary_data_predicted_pos_withsections[rev(y)[i], c(lowerCI, upperCI)], rep(y[i], 2), col="red", lwd=1)
}
points(rev(summary_data_predicted_pos_withsections[y, diff]), y, pch=pch_list, col="red", cex=0.5) # re-add dots on top of line
#axis(side=1, at=seq(xlim_1.2[1], xlim_1.2[2], 2), padj=-2.5, cex.axis=0.6)
axis(side=1, at=c(xlim_1.2[1], xlim_1.2[2]), tck=0, labels=FALSE)
axis(side=1, padj=-2.5, cex.axis=0.6)
# add names
sizes = rep(0.6, times=length(y)); sizes[is.na(summary_data_predicted_pos_withsections$Biomarker[y])] = 0.8
text(x=xlim_1.2[1]-space1[1], y=y, labels = rev(summary_data_predicted_pos_withsections$label[y]), pos=4, cex=rev(sizes))
# add bottom x-axis
mtext(side=1, at=(xlim_1.2[2] - abs(xlim_1.2[1]))/2, "Mean absolute difference between\n16 and 36 weeks in the control group (SD units)", line=1.25, cex=0.5)

dev.off()

pdf(paste0(mainpath, "Paper_TablesFigures/Figure2b.pdf"), width=7.48031, height=9.84252)
par(mfrow=c(1, 2), mar=c(2.5, 0.5, 0.5, 0.5), oma = c(0, 0, 0, 1))
### column 1
y = 81:118 # dim(summary_data_predicted_pos_withsections)[1]
pch_list = rep(16, length(y)); #pch_list[rev(summary_data_predicted_pos_withsections[y, "pvalue"])<pvalue_sig] = 16; #pch_list[rev(summary_data_predicted_pos_withsections[y, "pvalue"])<pvalue_BFsig] = 17; 
plot(rev(summary_data_predicted_pos_withsections[y, diff]), y, pch=pch_list, col="red", axes=FALSE, ylab="", xlab="", xlim=xlim_2.1 - space2, ylim=c(min(y)-0.5, max(y)+0.5), xaxs="i", yaxs="i", cex=0.5) #xlim=c(round(min_x, 2), round(max_x, 2)))
# add details
gr = unique(summary_data_predicted_pos_withsections$Group[y]); gr=gr[!is.na(gr)]
for (g in 1:length(gr)){
  # add solid line at 0
  #where = y[grep(gr[g], rev(summary_data_predicted_pos_withsections$Group[y]))]; lines(x=c(0,0), y=c(min(where)-0.5, max(where)+0.5))
  # add dashed lines
  where = y[grep(gr[g], rev(summary_data_predicted_pos_withsections$Group[y]))]; lines(x=c(0,0), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2)
  #where = y[grep(gr[g], rev(summary_data_predicted_pos_withsections$Group[y]))]; lines(x=rep(xlim_2.1[1],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2); lines(x=rep(xlim_2.1[2],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2)
  # add grey sections
  where = y[grep(gr[g], rev(summary_data_predicted_pos_withsections$Group[y]))]
  greylines = where[seq(1, length(where), 2)]
  rect(rep((xlim_2.1 - space2)[1], length(greylines)), greylines-0.5, rep((xlim_2.1 - space2)[2], length(greylines)), greylines+0.5, col=alpha("grey", 0.2), border=NA)
}
# add points and CI lines (on top of grey and vertical lines)
for (i in 1:length(y)){ #dim(summary_data_predicted_pos_withsections)[1]){
  lines(summary_data_predicted_pos_withsections[rev(y)[i], c(lowerCI, upperCI)], rep(y[i], 2), col="red", lwd=1)
}
points(rev(summary_data_predicted_pos_withsections[y, diff]), y, pch=pch_list, col="red", cex=0.5) # re-add dots on top of line
#axis(side=1, at=seq(xlim_2.1[1], xlim_2.1[2], 2), padj=-2.5, cex.axis=0.6)
axis(side=1, at=c(xlim_2.1[1], xlim_2.1[2]), tck=0, labels=FALSE)
axis(side=1, padj=-2.5, cex.axis=0.6)
# add names
sizes = rep(0.6, times=length(y)); sizes[is.na(summary_data_predicted_pos_withsections$Biomarker[y])] = 0.8
text(x=xlim_2.1[1]-space2[1], y=y, labels = rev(summary_data_predicted_pos_withsections$label[y]), pos=4, cex=rev(sizes))
# add bottom x-axis
mtext(side=1, at=(xlim_2.1[2] - abs(xlim_2.1[1]))/2, "Mean absolute difference between\n16 and 36 weeks in the control group (SD units)", line=1.25, cex=0.5)


### column 2
y = 119:161 # dim(summary_data_predicted_pos_withsections)[1]
pch_list = rep(16, length(y)); #pch_list[rev(summary_data_predicted_pos_withsections[y, "pvalue"])<pvalue_sig] = 16; #pch_list[rev(summary_data_predicted_pos_withsections[y, "pvalue"])<pvalue_BFsig] = 17; 
plot(rev(summary_data_predicted_pos_withsections[y, diff]), y, pch=pch_list, col="red", axes=FALSE, ylab="", xlab="", xlim=xlim_2.2 - space2, ylim=c(min(y)-0.5, max(y)+0.5), xaxs="i", yaxs="i", cex=0.5) #xlim=c(round(min_x, 2), round(max_x, 2)))
# add details
gr = unique(summary_data_predicted_pos_withsections$Group[y]); gr=gr[!is.na(gr)]
for (g in 1:length(gr)){
  # add solid line at 0
  #where = y[grep(gr[g], rev(summary_data_predicted_pos_withsections$Group[y]))]; lines(x=c(0,0), y=c(min(where)-0.5, max(where)+0.5))
  # add dashed lines
  where = y[grep(gr[g], rev(summary_data_predicted_pos_withsections$Group[y]))]; lines(x=c(0,0), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2)
  #where = y[grep(gr[g], rev(summary_data_predicted_pos_withsections$Group[y]))]; lines(x=rep(xlim_2.2[1],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2); lines(x=rep(xlim_2.2[2],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2)
  # add grey sections
  where = y[grep(gr[g], rev(summary_data_predicted_pos_withsections$Group[y]))]
  greylines = where[seq(1, length(where), 2)]
  rect(rep((xlim_2.2-space2)[1], length(greylines)), greylines-0.5, rep((xlim_2.2-space2)[2], length(greylines)), greylines+0.5, col=alpha("grey", 0.2), border=NA)
}
# add points and CI lines (on top of grey and vertical lines)
for (i in 1:length(y)){ #dim(summary_data_predicted_pos_withsections)[1]){
  lines(summary_data_predicted_pos_withsections[rev(y)[i], c(lowerCI, upperCI)], rep(y[i], 2), col="red", lwd=1)
}
points(rev(summary_data_predicted_pos_withsections[y, diff]), y, pch=pch_list, col="red", cex=0.5) # re-add dots on top of line
#axis(side=1, at=seq(xlim_2.2[1], xlim_2.2[2], 2), padj=-2.5, cex.axis=0.6)
axis(side=1, at=c(xlim_2.2[1], xlim_2.2[2]), tck=0, labels=FALSE)
axis(side=1, padj=-2.5, cex.axis=0.6)
# add names
sizes = rep(0.6, times=length(y)); sizes[is.na(summary_data_predicted_pos_withsections$Biomarker[y])] = 0.8
text(x=xlim_2.2[1]-space2[1], y=y, labels = rev(summary_data_predicted_pos_withsections$label[y]), pos=4, cex=rev(sizes))
# add bottom label
mtext(side=1, at=(xlim_2.2[2] - abs(xlim_2.2[1]))/2, "Mean absolute difference between\n16 and 36 weeks in the control group (SD units)", line=1.25, cex=0.5)

dev.off()

## negatives
pdf(paste0(mainpath, "Paper_TablesFigures/Figure2c.pdf"), width=7.48031/2, height=9.84252)
par(mar=c(2.5, 0.5, 0.5, 0.5), oma = c(0, 0, 0, 1))
### column 1
y = 1:41 # dim(summary_data_predicted_neg_withsections)[1]
pch_list = rep(16, length(y)); #pch_list[rev(summary_data_predicted_neg_withsections[y, "pvalue"])<pvalue_sig] = 16; #pch_list[rev(summary_data_predicted_neg_withsections[y, "pvalue"])<pvalue_BFsig] = 17; 
plot(rev(summary_data_predicted_neg_withsections[y, diff]), y, pch=pch_list, col="red", axes=FALSE, ylab="", xlab="", xlim=xlim_3.1 - space3, ylim=c(min(y)-0.5, max(y)+0.5), xaxs="i", yaxs="i", cex=0.5) #xlim=c(round(min_x, 2), round(max_x, 2)))
# add details
gr = unique(summary_data_predicted_neg_withsections$Group[y]); gr=gr[!is.na(gr)]
for (g in 1:length(gr)){
  # add solid line at 0
  #where = y[grep(gr[g], rev(summary_data_predicted_neg_withsections$Group[y]))]; lines(x=c(0,0), y=c(min(where)-0.5, max(where)+0.5))
  # add dashed lines
  where = y[grep(gr[g], rev(summary_data_predicted_neg_withsections$Group[y]))]; lines(x=c(0,0), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2)
  #where = y[grep(gr[g], rev(summary_data_predicted_neg_withsections$Group[y]))]; lines(x=rep(xlim_3.1[1],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2); lines(x=rep(xlim_3.1[2],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2)
  # add grey sections
  where = y[grep(gr[g], rev(summary_data_predicted_neg_withsections$Group[y]))]
  greylines = where[seq(1, length(where), 2)]
  rect(rep((xlim_3.1 - space3)[1], length(greylines)), greylines-0.5, rep((xlim_3.1 - space3)[2], length(greylines)), greylines+0.5, col=alpha("grey", 0.2), border=NA)
}
# add points and CI lines (on top of grey and vertical lines)
for (i in 1:length(y)){ #dim(summary_data_predicted_neg_withsections)[1]){
  lines(summary_data_predicted_neg_withsections[rev(y)[i], c(lowerCI, upperCI)], rep(y[i], 2), col="red", lwd=1)
}
points(rev(summary_data_predicted_neg_withsections[y, diff]), y, pch=pch_list, col="red", cex=0.5) # re-add dots on top of line
#axis(side=1, at=c(xlim_3.1[1], xlim_3.1[2]), labels=c("", ""), padj=-2.5, cex.axis=0.6)
#axis(side=1, at=seq(xlim_3.1[2], xlim_3.1[1], -2), padj=-2.5, cex.axis=0.6)
axis(side=1, at=c(xlim_3.1[1], xlim_3.1[2]), tck=0, labels=FALSE)
axis(side=1, padj=-2.5, cex.axis=0.6)
# add names
sizes = rep(0.6, times=length(y)); sizes[is.na(summary_data_predicted_neg_withsections$Biomarker[y])] = 0.8
text(x=xlim_3.1[1]-space3[1], y=y, labels = rev(summary_data_predicted_neg_withsections$label[y]), pos=4, cex=rev(sizes))
# add bottom label
mtext(side=1, at=(xlim_3.1[2] - abs(xlim_3.1[1]))/2, "Mean absolute difference between\n16 and 36 weeks in the control group (SD units)", line=1.25, cex=0.5)

dev.off()

########################
# Figure 3
########################
### for the Standardised model - therefore we must start analysis again

### clear workspace
rm(list = ls())

### load packages
library(scales)

### define path to main folder (contains subfolders called "RScripts", "Results" (with subfolder "Raw") and "Data")
mainpath <- "//ads.bris.ac.uk/filestore/BRMS/Research/Metabolic profiles in UPBEAT/"

### source useful functions
source(paste0(mainpath, "RScripts/UsefulFunctions_UPBEAT.R"))

### load results for standardised model
load(file=paste0(mainpath, "Paper_TablesFigures/Raw/LME_Model_Standardised_FullResults.Rdata"))
# Fit_Results_lme4, ERRORS_lme4, CI_error, 

### sort out the groupings
dictionary <- read.csv(paste0(mainpath, "Data/metabolome_data_dictionary.csv"), row.names=1)
# unfactor dictionary columns
for (i in 1:ncol(dictionary)){
  dictionary[, i] = as.character(dictionary[, i])
}

# trim white space from dictionary entries
dictionary$Group <- trim(dictionary$Group)
dictionary$Biomarker <- trim(dictionary$Biomarker)

# assign Group and Biomarker to Fit_results
Fit_Results_lme4$Group = dictionary[rownames(Fit_Results_lme4), "Group"]
Fit_Results_lme4$Biomarker = dictionary[rownames(Fit_Results_lme4), "Biomarker"]

METABOLITES = rownames(Fit_Results_lme4)

## find min and max 
min_x = min(c(Fit_Results_lme4[, "TI*time 2.5CI"], Fit_Results_lme4[, "TI*time 97.5CI"]), na.rm=TRUE)
max_x = max(c(Fit_Results_lme4[, "TI*time 2.5CI"], Fit_Results_lme4[, "TI*time 97.5CI"]), na.rm=TRUE)

file=load(paste0(mainpath, "RScripts/metabolite_labels.Rdata"))
metabolite_labels$label.no.units[metabolite_labels$metabolite=="tgpg"] = "Triglycerides/Phosphoglycerides (%)"

Fit_Results_lme4$label = Fit_Results_lme4$Biomarker
Fit_Results_lme4$label[Fit_Results_lme4$label%in%"glycerol"] = "Glycerol"
Fit_Results_lme4$label = metabolite_labels[Fit_Results_lme4$label, "label.no.units"]

## add FDR 
p_values = Fit_Results_lme4[, "pvalue"]
pvalues_order = order(p_values, decreasing =TRUE)
m_tests = length(p_values)
tmp = data.frame(p_value = p_values[pvalues_order], orig_order = (1:m_tests)[pvalues_order], sorted_order = (1:m_tests))
p_values_FDR = (tmp$p_value*m_tests)/(m_tests - (tmp$sorted_order - 1))
tmp = cbind(tmp, p_values_FDR=p_values_FDR)
Fit_Results_lme4$p_values_FDR = p_values_FDR[order(tmp$orig_order)]

# add * to label if FDR<0.05
FDR_005 = which(Fit_Results_lme4$p_values_FDR<0.05)
Fit_Results_lme4$label[FDR_005] = paste0(Fit_Results_lme4$label[FDR_005], "*")


## add a gap before each new group and plot nothing, but add group name
groups = unique(dictionary$Group)
for (g in 1:length(groups)){
  tmp = rbind(rep(NA, dim(Fit_Results_lme4)[2]), Fit_Results_lme4[Fit_Results_lme4$Group==groups[g], ])
  rownames(tmp) = c(groups[g], Fit_Results_lme4$Biomarker[Fit_Results_lme4$Group==groups[g]])
  tmp$label[1] = groups[g]
  if (g==1){
    Fit_Results_toplot = tmp
  } else{
    Fit_Results_toplot = rbind(Fit_Results_toplot, tmp)
  }
}

xlim_1.1 = xlim_1.2 = xlim_2.1 = xlim_2.2 = xlim_3.1 = c(-0.075, 0.075)
space = c(0.1, 0)

pvalue_sig = 0.05
#pvalue_BFsig = 0.05 / length(METABOLITES)

rownames(Fit_Results_toplot) = gsub("Concentration", "Conc.", rownames(Fit_Results_toplot))
rownames(Fit_Results_toplot) = gsub("and extremely", "& ext.", rownames(Fit_Results_toplot))

pdf(paste0(mainpath, "Paper_TablesFigures/Figure3a.pdf"), width=7.48031, height=9.84252)
par(mfrow=c(1, 2), mar=c(2.5, 0.5, 0.5, 0.5), oma = c(0, 0, 0, 1))
### column 1
y = 1:40 # dim(Fit_Results_toplot)[1]
pch_list = rep(16, length(y)); #pch_list[rev(Fit_Results_toplot[y, "pvalue"])<pvalue_sig] = 16; #pch_list[rev(Fit_Results_toplot[y, "pvalue"])<pvalue_BFsig] = 17; 
plot(rev(Fit_Results_toplot[y, "TI*time"]), y, pch=pch_list, col="red", axes=FALSE, ylab="", xlab="", xlim=xlim_1.1 - space, ylim=c(min(y)-0.5, max(y)+0.5), xaxs="i", yaxs="i") #xlim=c(round(min_x, 2), round(max_x, 2)))
# add details
gr = unique(Fit_Results_toplot$Group[y]); gr=gr[!is.na(gr)]
for (g in 1:length(gr)){
  # add solid line at 0
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]; lines(x=c(0,0), y=c(min(where)-0.5, max(where)+0.5))
  # add dashed lines
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]; lines(x=rep(xlim_1.1[1],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2); lines(x=rep(xlim_1.1[2],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2)
  # add grey sections
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]
  greylines = where[seq(1, length(where), 2)]
  rect(rep(-5, length(greylines)), greylines-0.5, rep(5, length(greylines)), greylines+0.5, col=alpha("grey", 0.2), border=NA)
}
# add points and CI lines (on top of grey and vertical lines)
for (i in 1:length(y)){ #dim(Fit_Results_toplot)[1]){
  lines(Fit_Results_toplot[rev(y)[i], c("TI*time 2.5CI", "TI*time 97.5CI")], rep(y[i], 2), col="red", lwd=2)
}
points(rev(Fit_Results_toplot[y, "TI*time"]), y, pch=pch_list, col="red") # re-add dots on top of line
axis(side=1, at=c(xlim_1.1[1], 0, xlim_1.1[2]), padj=-2.5, cex.axis=0.6)
# add names
sizes = rep(0.6, times=length(y)); sizes[is.na(Fit_Results_toplot$Biomarker[y])] = 0.8
text(x=xlim_1.1[1]-space[1], y=y, labels = rev(Fit_Results_toplot$label[y]), pos=4, cex=rev(sizes))
# add bottom x-axis
mtext(side=1, at=0, "Difference in mean change in outcome\nintervention vs control (SD units)", line=1.25, cex=0.5)


### column 2
y = 41:80 # dim(Fit_Results_toplot)[1]
pch_list = rep(16, length(y)); #pch_list[rev(Fit_Results_toplot[y, "pvalue"])<pvalue_sig] = 16; #pch_list[rev(Fit_Results_toplot[y, "pvalue"])<pvalue_BFsig] = 17; 
plot(rev(Fit_Results_toplot[y, "TI*time"]), y, pch=pch_list, col="red", axes=FALSE, ylab="", xlab="", xlim=xlim_1.2 - space, ylim=c(min(y)-0.5, max(y)+0.5), xaxs="i", yaxs="i") #xlim=c(round(min_x, 2), round(max_x, 2)))
# add details
gr = unique(Fit_Results_toplot$Group[y]); gr=gr[!is.na(gr)]
for (g in 1:length(gr)){
  # add solid line at 0
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]; lines(x=c(0,0), y=c(min(where)-0.5, max(where)+0.5))
  # add dashed lines
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]; lines(x=rep(xlim_1.2[1],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2); lines(x=rep(xlim_1.2[2],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2)
  # add grey sections
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]
  greylines = where[seq(1, length(where), 2)]
  rect(rep(-5, length(greylines)), greylines-0.5, rep(5, length(greylines)), greylines+0.5, col=alpha("grey", 0.2), border=NA)
}
# add points and CI lines (on top of grey and vertical lines)
for (i in 1:length(y)){ #dim(Fit_Results_toplot)[1]){
  lines(Fit_Results_toplot[rev(y)[i], c("TI*time 2.5CI", "TI*time 97.5CI")], rep(y[i], 2), col="red", lwd=2)
}
points(rev(Fit_Results_toplot[y, "TI*time"]), y, pch=pch_list, col="red") # re-add dots on top of line
axis(side=1, at=c(xlim_1.2[1], 0, xlim_1.2[2]), padj=-2.5, cex.axis=0.6)
# add names
sizes = rep(0.6, times=length(y)); sizes[is.na(Fit_Results_toplot$Biomarker[y])] = 0.8
text(x=xlim_1.2[1]-space[1], y=y, labels = rev(Fit_Results_toplot$label[y]), pos=4, cex=rev(sizes))
# add bottom x-axis
mtext(side=1, at=0, "Difference in mean change in outcome\nintervention vs control (SD units)", line=1.25, cex=0.5)

dev.off()



pdf(paste0(mainpath, "Paper_TablesFigures/Figure3b.pdf"), width=7.48031, height=9.84252)
par(mfrow=c(1, 2), mar=c(2.5, 0.5, 0.5, 0.5), oma = c(0, 0, 0, 1))
### column 1
y = 81:116 # dim(Fit_Results_toplot)[1]
pch_list = rep(16, length(y)); #pch_list[rev(Fit_Results_toplot[y, "pvalue"])<pvalue_sig] = 16; #pch_list[rev(Fit_Results_toplot[y, "pvalue"])<pvalue_BFsig] = 17; 
plot(rev(Fit_Results_toplot[y, "TI*time"]), y, pch=pch_list, col="red", axes=FALSE, ylab="", xlab="", xlim=xlim_2.1 - space, ylim=c(min(y)-0.5, max(y)+0.5), xaxs="i", yaxs="i") #xlim=c(round(min_x, 2), round(max_x, 2)))
# add details
gr = unique(Fit_Results_toplot$Group[y]); gr=gr[!is.na(gr)]
for (g in 1:length(gr)){
  # add solid line at 0
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]; lines(x=c(0,0), y=c(min(where)-0.5, max(where)+0.5))
  # add dashed lines
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]; lines(x=rep(xlim_2.1[1],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2); lines(x=rep(xlim_2.1[2],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2)
  # add grey sections
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]
  greylines = where[seq(1, length(where), 2)]
  rect(rep(-5, length(greylines)), greylines-0.5, rep(5, length(greylines)), greylines+0.5, col=alpha("grey", 0.2), border=NA)
}
# add points and CI lines (on top of grey and vertical lines)
for (i in 1:length(y)){ #dim(Fit_Results_toplot)[1]){
  lines(Fit_Results_toplot[rev(y)[i], c("TI*time 2.5CI", "TI*time 97.5CI")], rep(y[i], 2), col="red", lwd=2)
}
points(rev(Fit_Results_toplot[y, "TI*time"]), y, pch=pch_list, col="red") # re-add dots on top of line
axis(side=1, at=c(xlim_2.1[1], 0, xlim_2.1[2]), padj=-2.5, cex.axis=0.6)
# add names
sizes = rep(0.6, times=length(y)); sizes[is.na(Fit_Results_toplot$Biomarker[y])] = 0.8
text(x=xlim_2.1[1]-space[1], y=y, labels = rev(Fit_Results_toplot$label[y]), pos=4, cex=rev(sizes))
# add bottom x-axis
mtext(side=1, at=0, "Difference in mean change in outcome\nintervention vs control (SD units)", line=1.25, cex=0.5)


### column 2
y = 117:159 # dim(Fit_Results_toplot)[1]
pch_list = rep(16, length(y)); #pch_list[rev(Fit_Results_toplot[y, "pvalue"])<pvalue_sig] = 16; #pch_list[rev(Fit_Results_toplot[y, "pvalue"])<pvalue_BFsig] = 17; 
plot(rev(Fit_Results_toplot[y, "TI*time"]), y, pch=pch_list, col="red", axes=FALSE, ylab="", xlab="", xlim=xlim_2.2 - space, ylim=c(min(y)-0.5, max(y)+0.5), xaxs="i", yaxs="i") #xlim=c(round(min_x, 2), round(max_x, 2)))
# add details
gr = unique(Fit_Results_toplot$Group[y]); gr=gr[!is.na(gr)]
for (g in 1:length(gr)){
  # add solid line at 0
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]; lines(x=c(0,0), y=c(min(where)-0.5, max(where)+0.5))
  # add dashed lines
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]; lines(x=rep(xlim_2.2[1],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2); lines(x=rep(xlim_2.2[2],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2)
  # add grey sections
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]
  greylines = where[seq(1, length(where), 2)]
  rect(rep(-5, length(greylines)), greylines-0.5, rep(5, length(greylines)), greylines+0.5, col=alpha("grey", 0.2), border=NA)
}
# add points and CI lines (on top of grey and vertical lines)
for (i in 1:length(y)){ #dim(Fit_Results_toplot)[1]){
  lines(Fit_Results_toplot[rev(y)[i], c("TI*time 2.5CI", "TI*time 97.5CI")], rep(y[i], 2), col="red", lwd=2)
}
points(rev(Fit_Results_toplot[y, "TI*time"]), y, pch=pch_list, col="red") # re-add dots on top of line
axis(side=1, at=c(xlim_2.2[1], 0, xlim_2.2[2]), padj=-2.5, cex.axis=0.6)
# add names
sizes = rep(0.6, times=length(y)); sizes[is.na(Fit_Results_toplot$Biomarker[y])] = 0.8
text(x=xlim_2.2[1]-space[1], y=y, labels = rev(Fit_Results_toplot$label[y]), pos=4, cex=rev(sizes))
# add bottom label
mtext(side=1, at=0, "Difference in mean change in outcome\nintervention vs control (SD units)", line=1.25, cex=0.5)

dev.off()


pdf(paste0(mainpath, "Paper_TablesFigures/Figure3c.pdf"), width=7.48031/2, height=9.84252*(2/3))
par(mar=c(2.5, 0.5, 0.5, 0.5), oma = c(0, 0, 0, 1))
### column 1
y = 160:186 # dim(Fit_Results_toplot)[1]
pch_list = rep(16, length(y)); #pch_list[rev(Fit_Results_toplot[y, "pvalue"])<pvalue_sig] = 16; #pch_list[rev(Fit_Results_toplot[y, "pvalue"])<pvalue_BFsig] = 17; 
plot(rev(Fit_Results_toplot[y, "TI*time"]), y, pch=pch_list, col="red", axes=FALSE, ylab="", xlab="", xlim=xlim_3.1 - space, ylim=c(min(y)-0.5, max(y)+0.5), xaxs="i", yaxs="i") #xlim=c(round(min_x, 2), round(max_x, 2)))
# add details
gr = unique(Fit_Results_toplot$Group[y]); gr=gr[!is.na(gr)]
for (g in 1:length(gr)){
  # add solid line at 0
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]; lines(x=c(0,0), y=c(min(where)-0.5, max(where)+0.5))
  # add dashed lines
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]; lines(x=rep(xlim_3.1[1],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2); lines(x=rep(xlim_3.1[2],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2)
  # add grey sections
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]
  greylines = where[seq(1, length(where), 2)]
  rect(rep(-5, length(greylines)), greylines-0.5, rep(5, length(greylines)), greylines+0.5, col=alpha("grey", 0.2), border=NA)
}
# add points and CI lines (on top of grey and vertical lines)
for (i in 1:length(y)){ #dim(Fit_Results_toplot)[1]){
  lines(Fit_Results_toplot[rev(y)[i], c("TI*time 2.5CI", "TI*time 97.5CI")], rep(y[i], 2), col="red", lwd=2)
}
points(rev(Fit_Results_toplot[y, "TI*time"]), y, pch=pch_list, col="red") # re-add dots on top of line
axis(side=1, at=c(xlim_3.1[1], 0, xlim_3.1[2]), padj=-2.5, cex.axis=0.6)
# add names
sizes = rep(0.6, times=length(y)); sizes[is.na(Fit_Results_toplot$Biomarker[y])] = 0.8
text(x=xlim_3.1[1]-space[1], y=y, labels = rev(Fit_Results_toplot$label[y]), pos=4, cex=rev(sizes))
# add bottom label
mtext(side=1, at=0, "Difference in mean change in outcome\nintervention vs control (SD units)", line=1.25, cex=0.5)

dev.off()


########################
# Table 2 (a subset of sTable5) - SECOND COLUMN (FOR STANDARDISED MODEL) GENERATED HERE
########################

subset_groups = c("Extremely large VLDL", "Very large VLDL","Large VLDL", "Medium VLDL")
subset_otherLipids = c("Triglycerides in very large HDL", "Phospholipids in small HDL" ,"Mean diameter for VLDL particles","Serum total triglycerides",
                       "Triglycerides in VLDL" ,"Ratio of triglycerides to phosphoglycerides" ,"Estimated degree of unsaturation","Ratio of 18:2 linoleic acid to total fatty acids",
                       "Ratio of omega-6 fatty acids to total fatty acids","Ratio of polyunsaturated fatty acids to total fatty acids" ,"Ratio of saturated fatty acids to total fatty acids")
subset_otherTraits = c("Lactate", "Pyruvate", "Alanine" , "Acetate")

for (s in 1:length(subset_groups)){
  w = which(Fit_Results_toplot$Group%in%subset_groups[s])
  if (s==1){
    Table_Groups = Fit_Results_toplot[c(w[1]-1, w), ]
  } else{
    Table_Groups = rbind(Table_Groups, Fit_Results_toplot[c(w[1]-1, w), ])
  }
}

Table_otherLipids = Fit_Results_toplot[which(Fit_Results_toplot$Biomarker%in%subset_otherLipids), ]
Table_otherTraits = Fit_Results_toplot[which(Fit_Results_toplot$Biomarker%in%subset_otherTraits), ]

Whole_Table_standardised = rbind(Table_Groups,
                                 rep(NA, dim(Fit_Results_lme4)[2]), Table_otherLipids,
                                 rep(NA, dim(Fit_Results_lme4)[2]), Table_otherTraits)
rownames(Whole_Table_standardised)[33] = "Other lipids, lipoproteins affected by the intervention"
rownames(Whole_Table_standardised)[45] = "Other traits affected by the intervention"

Table2_col2 = cbind(format_nicely(Whole_Table_standardised[, "TI*time"], 3) ,
                    paste0("(", format_nicely(Whole_Table_standardised[, "TI*time 2.5CI"], 3), ", ", format_nicely(Whole_Table_standardised[, "TI*time 97.5CI"], 3), ")"),
                    NA)
Table2_col2[grep("NA", Table2_col2[, 1]), ] = rep("", 3)
Table2_col2 = cbind(paste(Table2_col2[, 1], Table2_col2[, 2]), Table2_col2[, 3])
Table2_col2 = Table2_col2[, 2:1]
# Difference in mean rate of change in metabolic traits per 4-weeks of gestation between 16 and 36 weeks comparing women receiving intervention to control group
colnames(Table2_col2) = c("In original units (see first column) per 4-weeks",
                          "In SD units per 4-weeks")
rownames(Table2_col2) = rownames(Whole_Table_standardised)

write.csv(Table2_col2, file=paste0(mainpath, "Paper_TablesFigures/Table2_col2.csv"), row.names=TRUE)

# produce whole Table
Table2_col1 <- read.csv(file=paste0(mainpath, "Paper_TablesFigures/Table2_col1.csv"), row.names=1, header=TRUE)

Table2 <- data.frame(as.character(Table2_col1[, 1]), as.character(Table2_col2[, 2]))
colnames(Table2) = c("In original units (see first column) per 4-weeks",
                          "In SD units per 4-weeks")
rownames(Table2) = rownames(Whole_Table_standardised)

write.csv(Table2, file=paste0(mainpath, "Paper_TablesFigures/Table2.csv"), row.names=TRUE)


########################################################################
## Figures for Supplementary model (spline, with knot point at 28 weeks)
########################################################################

### clear workspace
rm(list = ls())

### load packages
library(scales)

### define path to main folder (contains subfolders called "RScripts", "Results" (with subfolder "Raw") and "Data")
mainpath <- "//ads.bris.ac.uk/filestore/BRMS/Research/Metabolic profiles in UPBEAT/"

### source useful functions
source(paste0(mainpath, "RScripts/UsefulFunctions_UPBEAT.R"))

### load spline fit, standardised results
load(file=paste0(mainpath, "Paper_TablesFigures/Raw/LME_SplineModel_Standardised_FullResults.Rdata"))
# Fit_Results_lme4, ERRORS_lme4, CI_error, 

### sort out the groupings
dictionary <- read.csv(paste0(mainpath, "Data/metabolome_data_dictionary.csv"), row.names=1)
for (i in 1:ncol(dictionary)){
  dictionary[, i] = as.character(dictionary[, i])
}

# trim white space from dictionary entries
dictionary$Group <- trim(dictionary$Group)
dictionary$Biomarker <- trim(dictionary$Biomarker)

# assign Group and Biomarker to Fit_results
Fit_Results_lme4$Group = dictionary[rownames(Fit_Results_lme4), "Group"]
Fit_Results_lme4$Biomarker = dictionary[rownames(Fit_Results_lme4), "Biomarker"]

########################
# sFigure 2
########################

METABOLITES = rownames(Fit_Results_lme4)

## find min and max - these are really big though, in comparison to the rest
min_x = min(c(Fit_Results_lme4[, "TI*gest1 2.5CI"], Fit_Results_lme4[, "TI*gest1 97.5CI"],Fit_Results_lme4[, "TI*gest2 2.5CI"], Fit_Results_lme4[, "TI*gest2 97.5CI"]), na.rm=TRUE)
max_x = max(c(Fit_Results_lme4[, "TI*gest1 2.5CI"], Fit_Results_lme4[, "TI*gest1 97.5CI"],Fit_Results_lme4[, "TI*gest2 2.5CI"], Fit_Results_lme4[, "TI*gest2 97.5CI"]), na.rm=TRUE)

file=load(paste0(mainpath, "RScripts/metabolite_labels.Rdata"))
metabolite_labels$label.no.units[metabolite_labels$metabolite=="tgpg"] = "Triglycerides/Phosphoglycerides (%)"

Fit_Results_lme4$label = Fit_Results_lme4$Biomarker
Fit_Results_lme4$label[Fit_Results_lme4$label%in%"glycerol"] = "Glycerol"
Fit_Results_lme4$label = metabolite_labels[Fit_Results_lme4$label, "label.no.units"]

# add a gap before each new group and plot nothing, but add group name
groups = unique(dictionary$Group)
for (g in 1:length(groups)){
  tmp = rbind(rep(NA, dim(Fit_Results_lme4)[2]), Fit_Results_lme4[Fit_Results_lme4$Group==groups[g], ])
  rownames(tmp) = c(groups[g], Fit_Results_lme4$Biomarker[Fit_Results_lme4$Group==groups[g]])
  tmp$label[1] = groups[g]
  if (g==1){
    Fit_Results_toplot = tmp
  } else{
    Fit_Results_toplot = rbind(Fit_Results_toplot, tmp)
  }
}

xlim_1.1 = xlim_1.2 = xlim_2.1 = xlim_2.2 = xlim_3.1 = c(-0.1, 0.1)
space = c(0.1, 0)

pvalue_sig = 0.05
#pvalue_BFsig = 0.05 / length(METABOLITES)

rownames(Fit_Results_toplot) = gsub("Concentration", "Conc.", rownames(Fit_Results_toplot))
rownames(Fit_Results_toplot) = gsub("and extremely", "& ext.", rownames(Fit_Results_toplot))

pdf(paste0(mainpath, "Paper_TablesFigures/sFigure2a.pdf"), width=7.48031, height=9.84252)
par(mfrow=c(1, 2), mar=c(2.5, 0.5, 0.5, 0.5), oma = c(0, 0, 0, 1))
### column 1
y = 1:40 # dim(Fit_Results_toplot)[1]
pch_list1 = rep(16, length(y)); #pch_list1[rev(Fit_Results_toplot[y, "pvalue1"])<pvalue_sig] = 16; #pch_list1[rev(Fit_Results_toplot[y, "pvalue1"])<pvalue_BFsig] = 17; 
pch_list2 = rep(16, length(y)); #pch_list2[rev(Fit_Results_toplot[y, "pvalue2"])<pvalue_sig] = 16; #pch_list2[rev(Fit_Results_toplot[y, "pvalue2"])<pvalue_BFsig] = 17;
# diff slope 1
plot(rev(Fit_Results_toplot[y, "TI*gest1"]), y+0.25, pch=pch_list1, col="red", axes=FALSE, ylab="", xlab="", xlim=xlim_1.1 - space, ylim=c(min(y)-0.5, max(y)+0.5), xaxs="i", yaxs="i") #xlim=c(round(min_x, 2), round(max_x, 2)))
# diff slope 2
points(rev(Fit_Results_toplot[y, "TI*gest2"]), y-0.25, pch=pch_list2, col="blue")
# add details
gr = unique(Fit_Results_toplot$Group[y]); gr=gr[!is.na(gr)]
for (g in 1:length(gr)){
  # add solid line at 0
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]; lines(x=c(0,0), y=c(min(where)-0.5, max(where)+0.5))
  # add dashed lines
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]; lines(x=rep(xlim_1.1[1],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2); lines(x=rep(xlim_1.1[2],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2)
  # add grey sections
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]
  greylines = where[seq(1, length(where), 2)]
  rect(rep(-5, length(greylines)), greylines-0.5, rep(5, length(greylines)), greylines+0.5, col=alpha("grey", 0.2), border=NA)
}
# add points and CI lines (on top of grey and vertical lines)
points(rev(Fit_Results_toplot[y, "TI*gest1"]), y+0.25, pch=pch_list1, col="red") # re-add dots on top of line
points(rev(Fit_Results_toplot[y, "TI*gest2"]), y-0.25, pch=pch_list2, col="blue")
for (i in 1:length(y)){ #dim(Fit_Results_toplot)[1]){
  lines(Fit_Results_toplot[rev(y)[i], c("TI*gest1 2.5CI", "TI*gest1 97.5CI")], rep(y[i]+0.25, 2), col="red", lwd=2)
  lines(Fit_Results_toplot[rev(y)[i], c("TI*gest2 2.5CI", "TI*gest2 97.5CI")], rep(y[i]-0.25, 2), col="blue", lwd=2)
}
axis(side=1, at=c(xlim_1.1[1], 0, xlim_1.1[2]), padj=-2.5, cex.axis=0.6)
# add names
sizes = rep(0.6, times=length(y)); sizes[is.na(Fit_Results_toplot$Biomarker[y])] = 0.8
text(x=xlim_1.1[1]-space[1], y=y, labels = rev(Fit_Results_toplot$label[y]), pos=4, cex=rev(sizes))
# add bottom label
mtext(side=1, at=(xlim_1.1[2] - abs(xlim_1.1[1]))/2, "Difference in mean change in outcome\nintervention vs control (SD units)", line=1.25, cex=0.5)


### column 2
y = 41:80 # dim(Fit_Results_toplot)[1]
pch_list1 = rep(16, length(y)); #pch_list1[rev(Fit_Results_toplot[y, "pvalue1"])<pvalue_sig] = 16; #pch_list1[rev(Fit_Results_toplot[y, "pvalue1"])<pvalue_BFsig] = 17; 
pch_list2 = rep(16, length(y)); #pch_list2[rev(Fit_Results_toplot[y, "pvalue2"])<pvalue_sig] = 16; #pch_list2[rev(Fit_Results_toplot[y, "pvalue2"])<pvalue_BFsig] = 17;
# diff slope 1
plot(rev(Fit_Results_toplot[y, "TI*gest1"]), y+0.25, pch=pch_list1, col="red", axes=FALSE, ylab="", xlab="", xlim=xlim_1.1 - space, ylim=c(min(y)-0.5, max(y)+0.5), xaxs="i", yaxs="i") #xlim=c(round(min_x, 2), round(max_x, 2)))
# diff slope 2
points(rev(Fit_Results_toplot[y, "TI*gest2"]), y-0.25, pch=pch_list2, col="blue")
# add details
gr = unique(Fit_Results_toplot$Group[y]); gr=gr[!is.na(gr)]
for (g in 1:length(gr)){
  # add solid line at 0
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]; lines(x=c(0,0), y=c(min(where)-0.5, max(where)+0.5))
  # add dashed lines
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]; lines(x=rep(xlim_1.2[1],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2); lines(x=rep(xlim_1.2[2],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2)
  # add grey sections
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]
  greylines = where[seq(1, length(where), 2)]
  rect(rep(-5, length(greylines)), greylines-0.5, rep(5, length(greylines)), greylines+0.5, col=alpha("grey", 0.2), border=NA)
}
# add points and CI lines (on top of grey and vertical lines)
points(rev(Fit_Results_toplot[y, "TI*gest1"]), y+0.25, pch=pch_list1, col="red") # re-add dots on top of line
points(rev(Fit_Results_toplot[y, "TI*gest2"]), y-0.25, pch=pch_list2, col="blue")# re-add dots on top of line
for (i in 1:length(y)){ #dim(Fit_Results_toplot)[1]){
  lines(Fit_Results_toplot[rev(y)[i], c("TI*gest1 2.5CI", "TI*gest1 97.5CI")], rep(y[i]+0.25, 2), col="red", lwd=2)
  lines(Fit_Results_toplot[rev(y)[i], c("TI*gest2 2.5CI", "TI*gest2 97.5CI")], rep(y[i]-0.25, 2), col="blue", lwd=2)
}
axis(side=1, at=c(xlim_1.2[1], 0, xlim_1.2[2]), padj=-2.5, cex.axis=0.6)
# add names
sizes = rep(0.6, times=length(y)); sizes[is.na(Fit_Results_toplot$Biomarker[y])] = 0.8
text(x=xlim_1.2[1]-space[1], y=y, labels = rev(Fit_Results_toplot$label[y]), pos=4, cex=rev(sizes))
# add bottom label
mtext(side=1, at=(xlim_1.2[2] - abs(xlim_1.2[1]))/2, "Difference in mean change in outcome\nintervention vs control (SD units)", line=1.25, cex=0.5)


dev.off()


pdf(paste0(mainpath, "Paper_TablesFigures/sFigure2b.pdf"), width=7.48031, height=9.84252)
par(mfrow=c(1, 2), mar=c(2.5, 0.5, 0.5, 0.5), oma = c(0, 0, 0, 1))
### column 1
y = 81:116 # dim(Fit_Results_toplot)[1]
pch_list1 = rep(16, length(y)); #pch_list1[rev(Fit_Results_toplot[y, "pvalue1"])<pvalue_sig] = 16; #pch_list1[rev(Fit_Results_toplot[y, "pvalue1"])<pvalue_BFsig] = 17; 
pch_list2 = rep(16, length(y)); #pch_list2[rev(Fit_Results_toplot[y, "pvalue2"])<pvalue_sig] = 16; #pch_list2[rev(Fit_Results_toplot[y, "pvalue2"])<pvalue_BFsig] = 17;
# diff slope 1
plot(rev(Fit_Results_toplot[y, "TI*gest1"]), y+0.25, pch=pch_list1, col="red", axes=FALSE, ylab="", xlab="", xlim=xlim_2.1 - space, ylim=c(min(y)-0.5, max(y)+0.5), xaxs="i", yaxs="i") #xlim=c(round(min_x, 2), round(max_x, 2)))
# diff slope 2
points(rev(Fit_Results_toplot[y, "TI*gest2"]), y-0.25, pch=pch_list2, col="blue")
# add details
gr = unique(Fit_Results_toplot$Group[y]); gr=gr[!is.na(gr)]
for (g in 1:length(gr)){
  # add solid line at 0
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]; lines(x=c(0,0), y=c(min(where)-0.5, max(where)+0.5))
  # add dashed lines
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]; lines(x=rep(xlim_2.1[1],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2); lines(x=rep(xlim_2.1[2],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2)
  # add grey sections
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]
  greylines = where[seq(1, length(where), 2)]
  rect(rep(-5, length(greylines)), greylines-0.5, rep(5, length(greylines)), greylines+0.5, col=alpha("grey", 0.2), border=NA)
}
# add points and CI lines (on top of grey and vertical lines)
points(rev(Fit_Results_toplot[y, "TI*gest1"]), y+0.25, pch=pch_list1, col="red") # re-add dots on top of line
points(rev(Fit_Results_toplot[y, "TI*gest2"]), y-0.25, pch=pch_list2, col="blue") # re-add dots on top of line
for (i in 1:length(y)){ #dim(Fit_Results_toplot)[1]){
  lines(Fit_Results_toplot[rev(y)[i], c("TI*gest1 2.5CI", "TI*gest1 97.5CI")], rep(y[i]+0.25, 2), col="red", lwd=2)
  lines(Fit_Results_toplot[rev(y)[i], c("TI*gest2 2.5CI", "TI*gest2 97.5CI")], rep(y[i]-0.25, 2), col="blue", lwd=2)
}
axis(side=1, at=c(xlim_2.1[1], 0, xlim_2.1[2]), padj=-2.5, cex.axis=0.6)
# add names
sizes = rep(0.6, times=length(y)); sizes[is.na(Fit_Results_toplot$Biomarker[y])] = 0.8
text(x=xlim_2.1[1]-space[1], y=y, labels = rev(Fit_Results_toplot$label[y]), pos=4, cex=rev(sizes))
# add bottom label
mtext(side=1, at=(xlim_2.1[2] - abs(xlim_2.1[1]))/2, "Difference in mean change in outcome\nintervention vs control (SD units)", line=1.25, cex=0.5)


### column 2
y = 117:159 # dim(Fit_Results_toplot)[1]
pch_list1 = rep(16, length(y)); #pch_list1[rev(Fit_Results_toplot[y, "pvalue1"])<pvalue_sig] = 16; #pch_list1[rev(Fit_Results_toplot[y, "pvalue1"])<pvalue_BFsig] = 17; 
pch_list2 = rep(16, length(y)); #pch_list2[rev(Fit_Results_toplot[y, "pvalue2"])<pvalue_sig] = 16; #pch_list2[rev(Fit_Results_toplot[y, "pvalue2"])<pvalue_BFsig] = 17;
# diff slope 1
plot(rev(Fit_Results_toplot[y, "TI*gest1"]), y+0.25, pch=pch_list1, col="red", axes=FALSE, ylab="", xlab="", xlim=xlim_2.1 - space, ylim=c(min(y)-0.5, max(y)+0.5), xaxs="i", yaxs="i") #xlim=c(round(min_x, 2), round(max_x, 2)))
# diff slope 2
points(rev(Fit_Results_toplot[y, "TI*gest2"]), y-0.25, pch=pch_list2, col="blue")
# add details
gr = unique(Fit_Results_toplot$Group[y]); gr=gr[!is.na(gr)]
for (g in 1:length(gr)){
  # add solid line at 0
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]; lines(x=c(0,0), y=c(min(where)-0.5, max(where)+0.5))
  # add dashed lines
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]; lines(x=rep(xlim_2.2[1],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2); lines(x=rep(xlim_2.2[2],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2)
  # add grey sections
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]
  greylines = where[seq(1, length(where), 2)]
  rect(rep(-5, length(greylines)), greylines-0.5, rep(5, length(greylines)), greylines+0.5, col=alpha("grey", 0.2), border=NA)
}
# add points and CI lines (on top of grey and vertical lines)
points(rev(Fit_Results_toplot[y, "TI*gest1"]), y+0.25, pch=pch_list1, col="red") # re-add dots on top of line
points(rev(Fit_Results_toplot[y, "TI*gest2"]), y-0.25, pch=pch_list2, col="blue") # re-add dots on top of line
for (i in 1:length(y)){ #dim(Fit_Results_toplot)[1]){
  lines(Fit_Results_toplot[rev(y)[i], c("TI*gest1 2.5CI", "TI*gest1 97.5CI")], rep(y[i]+0.25, 2), col="red", lwd=2)
  lines(Fit_Results_toplot[rev(y)[i], c("TI*gest2 2.5CI", "TI*gest2 97.5CI")], rep(y[i]-0.25, 2), col="blue", lwd=2)
}
axis(side=1, at=c(xlim_2.2[1], 0, xlim_2.2[2]), padj=-2.5, cex.axis=0.6)
# add names
sizes = rep(0.6, times=length(y)); sizes[is.na(Fit_Results_toplot$Biomarker[y])] = 0.8
text(x=xlim_2.2[1]-space[1], y=y, labels = rev(Fit_Results_toplot$label[y]), pos=4, cex=rev(sizes))
# add bottom label
mtext(side=1, at=(xlim_2.2[2] - abs(xlim_2.2[1]))/2, "Difference in mean change in outcome\nintervention vs control (SD units)", line=1.25, cex=0.5)

dev.off()

## for part 3 only
xlim_1.1 = xlim_1.2 = xlim_2.1 = xlim_2.2 = xlim_3.1 = c(-0.2, 0.2)
space = c(0.1, 0)

pdf(paste0(mainpath, "Paper_TablesFigures/sFigure2c.pdf"), width=7.48031/2, height=9.84252*(2/3))
par(mar=c(2.5, 0.5, 0.5, 0.5), oma = c(0, 0, 0, 1))
### column 1
y = 160:186 # dim(Fit_Results_toplot)[1]
pch_list1 = rep(16, length(y)); #pch_list1[rev(Fit_Results_toplot[y, "pvalue1"])<pvalue_sig] = 16; #pch_list1[rev(Fit_Results_toplot[y, "pvalue1"])<pvalue_BFsig] = 17; 
pch_list2 = rep(16, length(y)); #pch_list2[rev(Fit_Results_toplot[y, "pvalue2"])<pvalue_sig] = 16; #pch_list2[rev(Fit_Results_toplot[y, "pvalue2"])<pvalue_BFsig] = 17;
# diff slope 1
plot(rev(Fit_Results_toplot[y, "TI*gest1"]), y+0.25, pch=pch_list1, col="red", axes=FALSE, ylab="", xlab="", xlim=xlim_1.1 - space, ylim=c(min(y)-0.5, max(y)+0.5), xaxs="i", yaxs="i") #xlim=c(round(min_x, 2), round(max_x, 2)))
# diff slope 2
points(rev(Fit_Results_toplot[y, "TI*gest2"]), y-0.25, pch=pch_list2, col="blue")
# add details
gr = unique(Fit_Results_toplot$Group[y]); gr=gr[!is.na(gr)]
for (g in 1:length(gr)){
  # add solid line at 0
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]; lines(x=c(0,0), y=c(min(where)-0.5, max(where)+0.5))
  # add dashed lines
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]; lines(x=rep(xlim_3.1[1],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2); lines(x=rep(xlim_3.1[2],2), y=c(min(where)-0.5, max(where)+0.5), col="grey", lty=2)
  # add grey sections
  where = y[grep(gr[g], rev(Fit_Results_toplot$Group[y]))]
  greylines = where[seq(1, length(where), 2)]
  rect(rep(-5, length(greylines)), greylines-0.5, rep(5, length(greylines)), greylines+0.5, col=alpha("grey", 0.2), border=NA)
}
# add points and CI lines (on top of grey and vertical lines)
points(rev(Fit_Results_toplot[y, "TI*gest1"]), y+0.25, pch=pch_list1, col="red") # re-add dots on top of line
points(rev(Fit_Results_toplot[y, "TI*gest2"]), y-0.25, pch=pch_list2, col="blue") # re-add dots on top of line
for (i in 1:length(y)){ #dim(Fit_Results_toplot)[1]){
  lines(Fit_Results_toplot[rev(y)[i], c("TI*gest1 2.5CI", "TI*gest1 97.5CI")], rep(y[i]+0.25, 2), col="red", lwd=2)
  lines(Fit_Results_toplot[rev(y)[i], c("TI*gest2 2.5CI", "TI*gest2 97.5CI")], rep(y[i]-0.25, 2), col="blue", lwd=2)
}
axis(side=1, at=c(xlim_3.1[1], 0, xlim_3.1[2]), padj=-2.5, cex.axis=0.6)
# add names
sizes = rep(0.6, times=length(y)); sizes[is.na(Fit_Results_toplot$Biomarker[y])] = 0.8
text(x=xlim_3.1[1]-space[1], y=y, labels = rev(Fit_Results_toplot$label[y]), pos=4, cex=rev(sizes))
# add bottom label
mtext(side=1, at=(xlim_3.1[2] - abs(xlim_3.1[1]))/2, "Difference in mean change in outcome\nintervention vs control (SD units)", line=1.25, cex=0.5)

dev.off()


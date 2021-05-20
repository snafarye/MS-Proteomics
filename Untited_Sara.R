
library(MSstats)
?MSstats

# read in the data file
PP <- protein_peptides
Proteins <- read.csv('proteins.csv')
Peptide <- read.csv('peptide.csv')

dim(PP);colnames(PP)
dim(Proteins);colnames(Proteins)
dim(Peptide);colnames(Peptide)


length(unique(data$Peptide))

PP.proposed <- dataProcess(raw = PP,
                                normalization = 'equalizeMedians',
                                summaryMethod = 'TMP',
                                censoredInt = "NA", cutoffCensored = "minFeature",
                                MBimpute = TRUE,
                                maxQuantileforCensored=0.999)




## https://github.com/MeenaChoi/MSstats_material/tree/master/example%20dataset/SRM_CPTAC
#########################################
## Part 1: Install MSstats (only need to install once)
#########################################
# 1. Please install all dependancy packages first.
install.packages(c("gplots","lme4","ggplot2","reshape","reshape2","data.table","Rcpp"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("limma","marray","preprocessCore","MSnbase"))

# 2. Install MSstats
install.packages(pkgs = "MSstats_3.0.6.tar.gz", repos = NULL, type ="source") 

# MSstats installation is complete ~~~

# Please exit R and re-enter R again. 
# Once re-enter R, please start with the below commands to exercute MSstats.


#########################################
## Part 2: prepare Input data for MSstats
#########################################
# Output report from Skyline (.sky) using report format "MSstats"


#########################################
## Part 3: Run MSstats (start here if you already install MSstats)
#########################################
# load the library
library(MSstats)

# Help file
?MSstats

# Input data
# Installing
install.packages("readr")
# Loading
library("readr")
raw <- read.table(file = "forMSstats_Study7.txt", header = TRUE, sep = ",")

#raw<-read.csv("forMSstats_Study7.csv")
head(raw)
dim(raw)
#find if any missing values
table(is.na(raw))

#=====================
# Function: dataProcess
# pre-processing data: quality control of MS runs
? dataProcess
quantData<-dataProcess(raw)
# output of dataProcess includes several data types.
names(quantData)
# the data after reformatting and normalization
head(quantData$ProcessedData)


# output QuantData
write.csv(quantData$ProcessedData,"Study7_QuantData.csv")

#=====================
# Function: dataProcessPlots
# visualization 
? dataProcessPlots
dataProcessPlots(data=quantData,type="ProfilePlot",address="Study7_")
dataProcessPlots(data=quantData,type="QCPlot",address="Study7_")
dataProcessPlots(data=quantData,type="ConditionPlot",address="Study7_")

dataProcessPlots(data = quantData, type="QCplot", ylimUp=35,
                 width=5, height=5)

#=====================
# Function: modelBasedQCPlots
# visualization for model-based quality control in fitting model.
?modelBasedQCPlots

modelBasedQCPlots(data=quantData$ModelQC,type="ResidualPlots",featureName=FALSE,address="Study7_")
modelBasedQCPlots(data=quantData$ModelQC,type="QQPlots",feature.QQPlot="byFeature",address="Study7_")


#=====================
# Function: groupComparison
# generate testing results of protein inferences across concentrations
?groupComparison
levels(quantData$ProcessedData$GROUP_ORIGINAL)
comparison<-matrix(c(0,0,0,0,0,1,-1),nrow=1)
row.names(comparison)<-"I-J"

resultOneComparison<-groupComparison(contrast.matrix=comparison,data=quantData)
resultOneComparison$ComparisonResult

#=====================
# Function: groupComparisonPlots
# visualization for testing results
?groupComparisonPlots

# testing with more than one comparisons
comparison1<-matrix(c(0,0,0,0,0,1,-1),nrow=1)
comparison2<-matrix(c(0,0,0,0,1,0,-1),nrow=1)
comparison3<-matrix(c(0,0,0,1,0,0,-1),nrow=1)
comparison<-rbind(comparison1,comparison2,comparison3)
row.names(comparison)<-c("I-J","H-J","G-J")

resultMultiComparisons<-groupComparison(contrast.matrix=comparison,data=quantData)
resultMultiComparisons$ComparisonResult

write.csv(resultMultiComparisons$ComparisonResult, file="Study7_SignificanceTestingResult.csv")


#=====================
# Function: groupComparisonPlots
# visualization for testing results
?groupComparisonPlots

# Visualization 1: Volcano plot
# (1) default setup: FDR cutoff = 0.05; fold change cutoff = NA
groupComparisonPlots(data=resultMultiComparisons$ComparisonResult,type="VolcanoPlot", ylimUp=70,address="Study7_Ex1_")

# (2) Both FDR cutoff = 0.05; fold change cutoff = 2
groupComparisonPlots(data=resultMultiComparisons$ComparisonResult,type="VolcanoPlot", ylimUp=70,FCcutoff=2,ProteinName=FALSE,address="Study7_Ex2_")

# Visualization 2: Heatmap (required more than one comparisons)
# (1) default setup: FDR cutoff = 0.05; fold change cutoff = NA
groupComparisonPlots(data=resultMultiComparisons$ComparisonResult,type="Heatmap",address="Study7_Ex1_")

# (2) Both FDR cutoff = 0.05; fold change cutoff = 2
groupComparisonPlots(data=resultMultiComparisons$ComparisonResult,type="Heatmap",FCcutoff=2,address="Study7_Ex2_")

# Visualization 3: Comparison plot
groupComparisonPlots(data=resultMultiComparisons$ComparisonResult,type="ComparisonPlot",address="Study7_")



#=====================
# Function: designSampleSize
# calulate number of biological replicates per group for your next experiment
?designSampleSize

result.sample<-designSampleSize(data=resultMultiComparisons$fittedmodel,numSample=TRUE,desiredFC=c(1.05,1.25),FDR=0.01,power=0.95)
result.sample


#=====================
# Function: designSampleSizePlots
# visualization for sample size calculation
?designSampleSizePlots

pdf("Study7_SampleSizePlot.pdf")
designSampleSizePlots(data=result.sample)
dev.off()


#=====================
# Function: quantification
# model-based quantification
?quantification

# (1) group quantification
groupQuant<-quantification(data=quantData$ProcessedData,type="Group")
write.csv(groupQuant, file="Study7_GroupQuantification.csv")

# (2) sample quantification
sampleQuant<-quantification(data=quantData$ProcessedData,type="Sample")
write.csv(sampleQuant, file="Study7_SampleQuantification.csv")








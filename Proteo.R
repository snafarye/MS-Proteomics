#Libraries and Install

library(devtools)
#For creating and useing other packages

library(BiocManager)
#For doing bioinformatics in R. Contains MSstats, topGO, and limma packages

#Installation of most current packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MSstats")
BiocManager::install("topGO")
BiocManager::install("limma")
BiocManager::install("viper")
BiocManager::install("pcaMethods")
BiocManager::install("impute")
BiocManager::install("UniProt.ws")
BiocManager::install("affy")
devtools::install_github("bartongroup/proteusLabelFree")
devtools::install_github("bartongroup/proteusTMT")
devtools::install_github("bartongroup/proteusSILAC")
devtools::install_github("bartongroup/Proteus", build_opts= c("--no-resave-data", "--no-manual"), build_vignettes=TRUE)
install_github("eahrne/SafeQuant")

library(MSstats)
#For comparison analysis of experimental and control samples

library(topGO)
#For significance testing of abudance and structural changes

library(limma)
#For linear analysis of microarrays data

library(proteus)
#For differential analysis

library(SafeQuant)
#For statistical validation of differential expression

library(pheatmap)
#For heat map generation

library(tidyverse)
#General Data Science and visualization Tools

library(reticulate)
#Allows for Python integration with R

py_install("Pandas")
pd=import("pandas")
#For data processing and containment

py_install("Seaborn")
Sb=import("seaborn")
#for data visualization

py_install("SciPy")
sc=import("scipy")
#For statistical comparison and regression

py_install("genewalk", pip=TRUE)
gw=import("genewalk")
#For quantifying relevant, context specific gene functions

#External Applications

#PyMOL
#For molecular visualization

#Spectronaut
#For spectral library generation

#Progenesis QI
#For small molecule discovery and analysis

PHP=read_csv("https://raw.githubusercontent.com/CharlesColgan/MS-Proteomics/main/Software%20Testing/phosphoPeptide.csv")
proteins=read_csv("https://raw.githubusercontent.com/CharlesColgan/MS-Proteomics/main/Software%20Testing/proteins.csv")
peptide=read_csv("https://raw.githubusercontent.com/CharlesColgan/MS-Proteomics/main/Software%20Testing/SampleDDAdata/SampleDDAdata/proteins.csv")
protein_peptides=read_csv("https://raw.githubusercontent.com/CharlesColgan/MS-Proteomics/main/Software%20Testing/SampleDDAdata/SampleDDAdata/protein-peptides.csv")

PHP$PTM=factor(PHP$PTM)
PHP$AScore=factor(PHP$AScore)
PHP$Accession=factor(PHP$Accession)
PHP$"Source File"=factor(PHP$"Source File")
summary(PHP)

#https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf

#Install
BiocManager::install("limma")

#Libraries
library(limma)
library(affy)

#Read in experiment design
targets=read.delim("~/Documents/R/limma_TEST/targets.txt")

#
eset=justRMA(filenames=targets$FileName)
colnames(eset)=row.names(targets)
head(exprs(eset))
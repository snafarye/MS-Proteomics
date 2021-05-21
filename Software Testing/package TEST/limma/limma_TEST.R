#https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf

#Install
BiocManager::install("limma")

#Libraries
library(limma)
library(affy)

#Read in experiment design
targets=read.delim("~/Documents/R/targets.txt")

#
esetLI=justRMA(filenames=targets$FileName)
colnames(esetLI)=row.names(targets)
head(exprs(esetLI))

#Gene annotation
Alias=sub("_.*","",row.names(esetLI))
fData(esetLI)=alias2SymbolUsingNCBI(Alias,
                                  "Escherichia_coli_str._K-12_substr._MG1655.gene_info.gz",
                                  required=c("GeneID","Symbol","type_of_gene"))
HasSymbol=!is.na(fData(esetLI)$Symbol)
esetLI=esetLI[HasSymbol,]

#Differential expression

#Create design matrix
Exp=factor(targets$Experiment)
Strain=factor(targets$Strain,levels=c("lrp+","lrp-"))
design=model.matrix(~Exp+Strain)

#Apply linear modeling
fit=lmFit(esetLI,design)

#Statistical analysis of linear fit
Bayes_fit=eBayes(fit,trend=TRUE,robust=TRUE)
results=decideTests(Bayes_fit)
summary(results)

#Statistical analysis for min. FC threshold
treatment=treat(fit,lfc=log2(2),trend=TRUE,robust=TRUE)
results=decideTests(treatment)
summary(results)

#Show top ranked results from stat analyses
topTable(Bayes_fit,n=30)
topTreat(treatment,n=30)
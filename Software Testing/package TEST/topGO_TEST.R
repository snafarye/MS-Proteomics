#https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf

#nessecary packages can be installed via BiocManager

#Libraries
library(topGO)
library(ALL)
library(genefilter)
library(lattice)

#Data
data(ALL)
data(geneList)

#GO annotations
affyLib=paste(annotation(ALL), "db", sep = ".")
library(package=affyLib,character.only=TRUE)

#see top expressed genes for a given cutoff
topDiffGenes=function(allScore){
  return(allScore<0.01)
}

#see total number of top expressed genes
sum(topDiffGenes(geneList))

#Create topGo object for enrichment analysis
sampleGOdata=new("topGOdata",
                 description="Simple session",
                 ontology="BP",
                 allGenes=geneList,
                 geneSel=topDiffGenes,
                 nodeSize=10,
                 annot=annFUN.db,
                 affyLib=affyLib)
sampleGOdata

#run Fisher Exact Test on gene counts
#Uses topGO object, GO graph struture, and test statistic
#two methods, classic and elim
#elim is more conservative in general
resultFisher=runTest(sampleGOdata,algorithm="classic",statistic="fisher")
resultFisher.elim=runTest(sampleGOdata,algorithm="elim",statistic="fisher")

resultFisher#classic method
resultFisher.elim#elim method

#Run Kolmogorov-Smirnov Test on gene scores
resultKS=runTest(sampleGOdata,algorithm="classic",statistic="ks")
resultKS.elim=runTest(sampleGOdata,algorithm="elim",statistic="ks")

resultKS#classic method
resultKS.elim#elim method

#Tabulate data summary of enrichment analysis
#Comparing KS vs Fisher and Classic vs Elim
#Shows the n top nodes based on orderBy
#ranksOf shows ranks in selected result
allRes=GenTable(sampleGOdata,
                classicFisher=resultFisher,
                elimFisher=resultFisher.elim,
                classicKS=resultKS,
                elimKS=resultKS.elim,
                orderBy="classicKS",
                ranksOf="classicFisher",
                topNodes=0)
allRes

#P-Values comparison of classic and elim methods for Fisher
pValue.classic=score(resultKS)
pValue.elim=score(resultKS.elim)[names(pValue.classic)]

#Plot elements
gstat=termStat(sampleGOdata, names(pValue.classic))
gSize=gstat$Annotated/max(gstat$Annotated)*4
colMap=function(x){
  .col=rep(rev(heat.colors(length(unique(x)))),time=table(x))
  return(.col[match(1:length(x), order(x))])
}
gCol=colMap(gstat$Significant)

#Scatter Plot
plot(pValue.classic,
     pValue.elim,
     xlab="p-value classic",
     ylab="p-value elim",
     pch=19,
     cex=gSize,col=gCol)

#GO terms where elim is less conservative than classic
sel.go=names(pValue.classic)[pValue.elim < pValue.classic]
cbind(termStat(sampleGOdata, sel.go),
      elim=pValue.elim[sel.go],
      classic=pValue.classic[sel.go])

#GO subgraph
#based on p-values and the n most significant terms
showSigOfNodes(sampleGOdata,
               termsP.value=score(resultKS.elim),
               firstSigNodes=5,
               useInfo='all')

#Filter out low expression values and small variability across samples

filter=filterfun(pOverA(0.20,log2(100)),function(x)(IQR(x)>0.25))
#filters out gene if the proportion of expression values
#larger than a given cutoff exceeds a given threshold
#filter out genes under the 1st quartile

selProbes=genefilter(ALL,filter)
esetGO=ALL[selProbes,]
esetGO

#Get P-values from chosen statistics for interresting genes
y=as.integer(sapply(esetGO$BT,function(x) return(substr(x,1,1)=='T')))
#One-hot encode between two class labels

geneList=getPvalues(exprs(esetGO),
                    classlabel=y,
                    test="t",
                    alternative="greater",
                    correction="BH")
#can now be used to make a topGOdata object
#Also possible to use pre-made GO-gene mappings

#Creating topGOdata object can include ontology types
#possible 
#geneSel is the criteria for selecting the genes
GOdata=new("topGOdata",
              description="GO analysis of ALL data; B-cell vs T-cell",
              ontology="BP",
              allGenes=geneList,
              geneSel=topDiffGenes,
              annot=annFUN.db,
              nodeSize=5,
              affyLib=affyLib)

#Adjusting number of probes to use for best test results
#based on filter we created earlier
allProb=featureNames(ALL)#All Probe names
groupProb=integer(length(allProb))+1
groupProb[allProb %in% genes(GOdata)]=0
groupProb[!selProbes]=2
groupProb=factor(groupProb,labels=c("Used","Not annotated","Filtered"))
table(groupProb)

#Visualize groups of probes
pValue=getPvalues(exprs(ALL),classlabel=y,alternative="greater")
geneVar=apply(exprs(ALL),1,var)
dd=data.frame(x=geneVar[allProb],y=log10(pValue[allProb]),groups=groupProb)
xyplot(y~x|groups,data=dd,groups=groups)

#topGOdata object uses
description(GOdata)#description, can be modified
genes(GOdata)#names of genes
numGenes(GOdata)#number of genes
geneScore(GOdata,whichGenes=sample(genes(GOdata),10))#get scores for a selection of genes
sigGenes(GOdata)#list of the significant genes
numSigGenes(GOdata)#number of significant genes
updateGenes(GOdata,geneScore(GOdata,use.names=TRUE),topDiffGenes)#update gene list
graph(GOdata)#show GO graph structure
usedGO(GOdata)#shows GO terms used
sel.terms=sample(usedGO(GOdata),10)#generic selection of terms
genesInTerm(GOdata,sel.terms)#genes annotated to a selection of terms
countGenesInTerm(GOdata,sel.terms)#number of annotated genes
scoresInTerm(GOdata,sel.terms)#scores of annotated genes
termStat(GOdata,sel.terms)#stats for selected terms
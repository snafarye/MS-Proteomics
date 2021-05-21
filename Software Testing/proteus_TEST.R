#http://www.compbio.dundee.ac.uk/user/mgierlinski/proteus/proteus.html

#Install
devtools::install_github("bartongroup/proteusLabelFree")

#Library
library(proteusLabelFree)
library(proteus)
data(proteusLabelFree)

#Input Evidence and Metadata
evidenceFile=system.file("extdata","evidence.txt.gz",package="proteusLabelFree")
metadataFile=system.file("extdata","metadata.txt",package="proteusLabelFree")
evi=readEvidenceFile(evidenceFile)
meta=read.delim(metadataFile,header=TRUE,sep="\t")

#base columns needed in evidence, can be changed and more added
#additionally, min. one column of measurements needed
str(evidenceColumns)
#metadata requires columns labeled experiment, measure, sample, and condition
#experiment columns in evidence and metadata should match
#metadata measure should contain name of measure column in evidence
#sample is a unique name
#condition is condition names

#Create peptide database
pepdat=makePeptideTable(evi,meta)
summary(pepdat)

#Create protein database
prodat=makeProteinTable(pepdat)
summary(prodat)

#Counts
plotCount(pepdat)
plotCount(prodat)

#Distance
plotDistanceMatrix(pepdat)
plotDistanceMatrix(prodat)

#PCA
plotPCA(pepdat)
plotPCA(prodat)

#Normalize data, can use other norm functions
normalizeData(prodat,norm.fun=normalizeMedian)
normalizeData(pepdat,norm.fun=normalizeMedian)

#Create protein group database
proteinGroupsFile=system.file("extdata","proteinGroups.txt.gz",package="proteusLabelFree")
prot.MQ=readProteinGroups(proteinGroupsFile,meta)

#Differential expression, passes through limma package
#defualts B-H adjust
res=limmaDE(prodat,conditions=c("A","B"),transform.fun=log2,sig.level=0.05)






#!/usr/bin/env Rscript

##author: Luke Hebert
##date begun: 11/28/2019
##description:
  ##takes a flat file of absolute RNASeq counts (including genes in first column, samplesIDs in 1st row, and condition status in 2nd row)
  ##performs differential expression analysis comparing conditions using DESeq2 package
  ##outputs plots and results to local directory

##useful resources
  ##somewhat useful for installing packages: https://www.youtube.com/watch?v=7UKMU5HK380
  ##widely useful, esp. for plots: https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf
  ##useful for PCA labels: https://support.bioconductor.org/p/90791/

#load DESeq2 package
library(DESeq2)
library(ggplot2)

#get arguments from command line
args = commandArgs(trailingOnly=TRUE)

#get the counts and the conditions input file names
countsFilename <- args[1]
conditionsFilename <- args[2]

###countsFilename <- "/Users/lhebert2-l/Desktop/Kathryn\ Kuehn/DifferentialExpressionAnalysis_Results2/cellM/countMatrix_cellM.csv"
###conditionsFilename <- "/Users/lhebert2-l/Desktop/Kathryn\ Kuehn/DifferentialExpressionAnalysis_Results2/cellM/conditionMatrix_cellM.csv"

#read counts file as a matrix
countsDataframe <- read.csv(countsFilename, sep=",", row.names=1, header = TRUE)
countsMatrix <- as.matrix(countsDataframe)

#read conditions file as a matrix
conditionsDataframe <- read.csv(conditionsFilename, sep=",", row.names=1, header = TRUE)
conditionsMatrix <- as.matrix(conditionsDataframe)

#convert matrices to a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData= countsMatrix, colData=conditionsMatrix, design=~condition)

#perform DE analysis for each gene using
dds <- DESeq(dds)

#create a character string for labeling all output which identifies what kind of input was used; this includes location/pathway characters
outLabel_temp <- gsub(".csv","_",countsFilename,fixed=TRUE)
outLabel <- gsub("countMatrix_","",outLabel_temp,fixed=TRUE)

#plot sparsity
dds_esf <- estimateSizeFactors(dds)
png(paste(outLabel,"plotSparsity.png",sep=""))
plotSparsity(dds_esf)
dev.off()
#plot dispersion estimate
dds_de <- estimateDispersions(dds)
png(paste(outLabel,"plotDispEst.png",sep=""))
plotDispEsts(dds_de)
dev.off()
#plot principal component analysis of rlog transformed data
rld <- rlog(dds)
z <- plotPCA(rld)
nudge <- position_nudge(y = 1)
png(paste(outLabel,"plotPCA.png",sep=""))
z + geom_text(aes(label = name), position = nudge)
# + xlim(-10,12.5)
dev.off()

if (grepl("pt", countsFilename, fixed=TRUE)) 
{
  #apply the DESeq2 pipeline so that it compares Spina Bifida samples to Control samples
  resultsDF <- results(dds, contrast = c("condition","SB","CTRL"))
  #save the results as a flat file
  write.csv(resultsDF, paste(outLabel,"resultsSB.csv",sep=""))
  #plot MA
  png(paste(outLabel,"plotMA.png",sep=""))
  plotMA(resultsDF)
  dev.off()
  
} else if (grepl("cell", countsFilename, fixed=TRUE))
{
  #apply the DESeq2 pipeline so that it compares the NO folic acid cells to NORMAL folic acid cells
  resultsDF1 <- results(dds, contrast = c("condition","NO","NORMAL"))
  #save the results as a flat file
  write.csv(resultsDF1, paste(outLabel,"resultsNoFA.csv",sep=""))
  #plot MA
  png(paste(outLabel,"plotNoMA.png",sep=""))
  plotMA(resultsDF1)
  dev.off()
  
  #apply the DESeq2 pipeline so that it compares the HIGH folic acid cells to NORMAL folic acid cells
  resultsDF2 <- results(dds, contrast = c("condition","HIGH","NORMAL"))
  #save the results as a flat file
  write.csv(resultsDF2, paste(outLabel,"resultsHighFA.csv",sep=""))
  #plot MA
  png(paste(outLabel,"plotHighMA.png",sep=""))
  plotMA(resultsDF2)
  dev.off()
}


# R for Beginners Vignette: DiffBind

This is the vignette that will accompany the R script, where I go through line by line how the package works.

This script requires: Bam files for input and chip replicates, and the Peak files (MACS2).

First, lines with the "#" sign in front of the script are comments and will not run in the terminal.

# Install packages
```
#Install Diffbind
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DiffBind")
```

This will install the DiffBind package. There will likely be other packages that DiffBind will auto-install that are required. To check if the installation happened correctly, try loading the package. If there was a problem with the installation, the package will not load and you will get an error. The error messages are usually pretty self explanatory, or can be googled pretty easily.

```
#Load Diffbind
library(DiffBind)
```
Next, we need to load a few other packages for this vignette. You will likely need to install these, and packages that are required by these packages. 

```
#Load other packages: 
library(ChIPseeker)
library(clusterProfiler)
library(ggplot2)
library(VennDiagram)
```

Next, we need to load the genome that we are using. I use hg19 normally, but this can be changed to hg38 or mm10. The txdb is a UCSC genome package that we need to install and load. 

```
#Load the genome
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#for mouse, the package to install and load:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
#https://bioconductor.org/packages/release/data/annotation/html/TxDb.Mmusculus.UCSC.mm10.knownGene.html

library("TxDb.Mmusculus.UCSC.mm10.knownGene")
```
# Begin the analysis
## Import your information

```
#Import mm10 information
txdb<- TxDb.Mmusculus.UCSC.mm10.knownGene
```

```
#Import your data
peaks<- dba(sampleSheet="metadata.csv")
samples<- read.csv("metadata.csv")
```

Your data should be formatted as a csv file and should be located in the working directory. I have provided a template here:
[a relative link](DiffBind_example_metadata.csv)

# Analyse the data

A first look at the data, uncorrected and normalized. This will pull up a clustered heatmap of all the replicates.
```
plot(peaks)
```
Next, we are going to count the reads and calculate a binding matrix with scores based on read counts for every sample (affinity scores). For narrowpeaks, it is advisable to use the "summits" option to re-center each peak aroundthe point of greatest enrichment.

```
peaks<- dba.count(peaks, summits = 250)
plot(peaks)
```
Next, we are going to establish a contrast (minmembers = the lowest number of replicates)

```
peaks<- dba.contrast(peaks, categories = DBA_TISSUE, minMembers = 2)
```

#Perform the analysis:
peaks<- dba.analyze(peaks)
plot(peaks, contrast = 1)

peaks.db<- dba.report(peaks)
dba.plotPCA(peaks,DBA_TISSUE, label=DBA_TISSUE)
dba.plotPCA(peaks,contrast=1, label=DBA_TISSUE)
dba.plotMA(peaks)
```




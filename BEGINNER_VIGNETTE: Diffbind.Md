<b>R for Beginners Vignette: DiffBind</b><br>

This is the vignette that will accompany the R script, where I go through line by line how the package works.

This script requires: Bam files for input and chip replicates, and the Peak files (MACS2).

First, lines with the "#" sign in front of the script are comments and will not run in the terminal.

#Install Diffbind

`if (!requireNamespace("BiocManager", quietly = TRUE)) \
  install.packages("BiocManager") \
BiocManager::install("DiffBind")`

This will install the DiffBind package. There will likely be other packages that DiffBind will auto-install that are required. 

`#Load Diffbind
library(DiffBind)`

If there was a problem with the installation, the package will not load and you will get an error. The error messages are usually pretty self explanatory, or can be googled pretty easily.

`#Load other packages that I have already installed. You will likely need to install these, and packages that are required by these packages. 
library(ChIPseeker)
library(clusterProfiler)
library(ggplot2)
library(VennDiagram)`

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#for mouse, the package to install and load:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
#https://bioconductor.org/packages/release/data/annotation/html/TxDb.Mmusculus.UCSC.mm10.knownGene.html

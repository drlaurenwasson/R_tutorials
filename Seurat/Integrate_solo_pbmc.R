#!/usr/bin/env Rscript
#This script was written by Lauren Wasson on 2021_09_21. It will load the 10X object, create a seurat object, add the solo score, and save the output as both an h5ad and a Seurat object
##Usage: Rscript Integrate_solo_pmbc.R <working_directory> <keyfile>
args = commandArgs(trailingOnly=TRUE)


# test to make sure the things are installed
if("remotes" %in% rownames(installed.packages()) == FALSE) {install.packages("remotes")}
library(remotes)
if("SeuratDisk" %in% rownames(installed.packages()) == FALSE) {remotes::install_github("mojaveazure/seurat-disk")}
if("SeuratData" %in% rownames(installed.packages()) == FALSE) {devtools::install_github('satijalab/seurat-data')}

#Load libraries
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(tidyverse)

args<- commandArgs(trailingOnly=TRUE)
setwd(args[1])
samples<- read.table(args[2])
identifiers = as.character(samples$V1)


for (i in 1:length(identifiers)) {
filename = paste(identifiers[i], ".pbmc", sep="")
filename = gsub("-","_",filename)
matrix = paste("/n/data1/hms/genetics/seidman/DeLaughter/pipeline/",identifiers[i],"/outs/filtered_feature_bc_matrix", sep="")
pbmc.data<- Read10X(data.dir=matrix)
pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells = 0, min.features = 10, project = identifiers[i])
solofile = paste("/n/data1/hms/genetics/seidman/DeLaughter/pipeline/",identifiers[i],"/outs/", identifiers[i], "_pbmc_withsolo_METADATA.csv", sep="")
solo <- read.delim(solofile, sep = ",")
pbmc[["solo_score"]] <- solo$solo_score
pbmc[["predicted_doublets_solo"]]<- solo$predicted_doublets_solo
assign(noquote(filename), pbmc)
remove(pbmc)
}


pbmcs<- mget(grep(pattern = ".pbmc", x = ls(),
                       value = TRUE))

pbmcstomerge<- pbmcs[-1]

pbmc<-merge(x = pbmcs[[1]], y = pbmcstomerge)

objectname<- str_remove(basename(args[2]), "_keyfile.txt")
finame<- paste(objectname,"object_",ncol(pbmc),"cells.RData", sep = "")
save(pbmc, file = finame)


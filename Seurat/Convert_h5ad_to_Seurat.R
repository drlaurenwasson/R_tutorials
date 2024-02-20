#Conversion of h5ad single cell RNA_seq atlases to Seurat Objects
#20230126 LWasson

#Install packages and load libraries

install.packages("Seurat")
install.packages("devtools")
devtools::install_github('satijalab/seurat-data')

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(anndata)

library(reticulate)
ad <- import("anndata", convert = FALSE)

#Adipocytes
Convert("scRNA-seq/hca_heart_adipocytes_raw.h5ad", dest = "h5seurat", overwrite = TRUE)
adipo <- LoadH5Seurat("scRNA-seq/hca_heart_adipocytes_raw.h5seurat")
adipo
save(adipo, file = "HCA_Atlas_Adipocytes.RData")

#Fibroblasts
Convert("scRNA-seq/hca_heart_fibroblasts_raw.h5ad", dest = "h5seurat", overwrite = TRUE)
fibro <- LoadH5Seurat("scRNA-seq/hca_heart_fibroblasts_raw.h5seurat")
fibro
save(fibro, file = "HCA_Atlas_Fibroblasts.RData")
fibrosp<- subset(fibro, subset = region == "SP")
save(fibrosp, file = "HCA_Atlas_Fibroblasts_SP.RData")

# Immune cells
Convert("hca_heart_immune_download.h5ad", dest = "h5seurat", overwrite = TRUE)
immune <- LoadH5Seurat("hca_heart_immune_download.h5seurat")
immune
save(immune, file = "HCA_Atlas_ImmuneCells.RData")
immunesp<- subset(immune, subset = region == "SP")
save(immunesp, file = "HCA_Atlas_ImmuneCells_SP.RData")

# Neuronal cells
Convert("hca_heart_neuronal_raw.h5ad", dest = "h5seurat", overwrite = TRUE)
neuro <- LoadH5Seurat("hca_heart_neuronal_raw.h5seurat")
neuro
save(neuro, file = "HCA_Atlas_NeuroCells.RData")
neurosp<- subset(neuro, subset = region == "SP")
save(neurosp, file = "HCA_Atlas_NeuroCells_SP.RData")

# Vascular cells
Convert("hca_heart_vascular_raw.h5ad", dest = "h5seurat", overwrite = TRUE)
vascular <- LoadH5Seurat("hca_heart_vascular_raw.h5seurat")
vascular
save(vascular, file = "HCA_Atlas_VascularCells.RData")
vascularsp<- subset(vascular, subset = region == "SP")
save(vascularsp, file = "HCA_Atlas_VascularCells_SP.RData")

# Ventricular CM
Convert("hca_heart_ventricular_CM_raw.h5ad", dest = "h5seurat", overwrite = TRUE)
cm <- LoadH5Seurat("hca_heart_ventricular_CM_raw.h5seurat")
cm
save(cm, file = "HCA_Atlas_Ventricular_CM.RData")
cmsp<- subset(cm, subset = region == "SP")
save(cmsp, file = "HCA_Atlas_Ventricular_CM_SP.RData")



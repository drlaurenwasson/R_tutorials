#Install (this only has to be done once)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

#Load (this must be done every time)
library(DESeq2)

#Load additional libraries that might be useful (these may need to be installed)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(VennDiagram)
library(gProfileR)

counts_file<- read.table("CHD7_data_1.txt")

#remove non-protein coding genes from analysis
protein_genes<- read.table("protein_genes.txt")
protein_genes<- as.character(protein_genes$V1)

counts_protein_coding <- counts[rownames(counts) %in% protein_genes,]

metadata<- read.csv("metadata.csv")
all(names(counts_protein_coding) %in% rownames(metadata))

#DDS by "exp"
dds<- DESeqDataSetFromMatrix(countData = counts_protein_coding, colData = metadata, design = ~ exp )

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds<- estimateSizeFactors(dds)
sizeFactors(dds)
combined_normalized_counts<- counts(dds, normalized =TRUE)
write.table(combined_normalized_counts, file= "dds_combined_normalized_counts_proteingenes.txt", sep = "\t", quote = FALSE)

dds$exp<- relevel(dds$exp, ref = "wt" 

rld <- rlog(dds, blind = TRUE)

#Plot PCA
tiff(filename = "PCA_analysis.tiff",
     width = 480 ,height = 480, units = "px", pointsize = 12)
plotPCA(rld, intgroup="exp")
dev.off()

pcaData <- plotPCA(rld, intgroup="exp", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
mycolors<- c("#209324","#ef0000","#1166d6","#e27106", "#e810aa", "#000000")
pcaplot<- ggplot(pcaData, aes(PC1,PC2, color = exp), scale_color_manual(values = mycolors)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
pcaplot+ scale_color_manual(values = mycolors)

cortest<- cor(counts_protein_coding, use = "all.obs", method = "pearson")
cortest
write.csv(cortest, file= "corrtest.csv", row.names= TRUE, col.names = TRUE)


dds <- DESeq(dds)
plotDispEsts(dds)

rld <- rlog(dds, blind=F)
topVarianceGenes <- head(order(rowVars(assay(rld)), decreasing=T),1000)
matrix <- assay(rld)[ topVarianceGenes, ]
matrix <- matrix - rowMeans(matrix)
#df<- as.data.frame(metadata(dds))
annotation_col<- as.data.frame(metadata$exp)
row.names(annotation_col)<- row.names(metadata)
colnames(annotation_col)<- "exp"
pheatmap(matrix, clustering_distance_rows = "correlation", annotation_col = annotation_col, show_rownames  = FALSE, show_colnames = FALSE, filename = "2019_01_10_top1000misexpressedheatmap_norows_nocols.pdf")
dev.off()
pheatmap(matrix, clustering_distance_rows = "correlation", annotation_col = annotation_col, cellwidth = 15, cellheight = 12, fontsize = 8, filename = "2019_01_10_top1000misexpressedheatmap_genenames.pdf")
dev.off()

#Top 250 genes
top250VarianceGenes <- head(order(rowVars(assay(rld)), decreasing=T),250)
matrix2 <- assay(rld)[ top250VarianceGenes, ]
matrix2 <- matrix2 - rowMeans(matrix2)
#df<- as.data.frame(metadata(dds))
annotation_col<- as.data.frame(metadata$exp)
row.names(annotation_col)<- row.names(metadata)
colnames(annotation_col)<- "exp"
pheatmap(matrix2, clustering_distance_rows = "correlation", annotation_col = annotation_col, show_rownames  = FALSE, show_colnames = FALSE, filename = "2019_01_10_top250misexpressedheatmap_norows_nocols.pdf")
dev.off()
pheatmap(matrix2, clustering_distance_rows = "correlation", annotation_col = annotation_col, cellwidth = 15, cellheight = 12, fontsize = 8, filename = "2019_01_10_top250misexpressedheatmap_genenames.pdf")
dev.off()

pca<- princomp(combined_normalized_counts)
#grab the top 50 genes from the top 4 components
write.table(row.names(as.data.frame(sort(abs(pca$scores[,"Comp.1"]),decreasing=TRUE)[1:50])), file = "PCA_component1", quote = F, row.names = FALSE)
pca.1<- as.character(read.table("PCA_component1")$V1)
write.table(row.names(as.data.frame(sort(abs(pca$scores[,"Comp.2"]),decreasing=TRUE)[1:50])), file = "PCA_component2", quote = F, row.names = FALSE)
pca.2<- as.character(read.table("PCA_component2")$V1)
write.table(row.names(as.data.frame(sort(abs(pca$scores[,"Comp.3"]),decreasing=TRUE)[1:50])), file = "PCA_component3", quote = F, row.names = FALSE)
pca.3<- as.character(read.table("PCA_component3")$V1)
write.table(row.names(as.data.frame(sort(abs(pca$scores[,"Comp.4"]),decreasing=TRUE)[1:50])), file = "PCA_component4", quote = F, row.names = FALSE)
pca.4<- as.character(read.table("PCA_component4")$V1)


res_table<- results(dds, contrast = c("group", "CHD7_het", "wt"))
summary(res_table)

#Set thresholds
padj.cutoff <- 0.01
lfc.cutoff <- 0.58

threshold<- res_table$padj < padj.cutoff & abs(res_table$log2FoldChange) > lfc.cutoff
length(which(threshold == TRUE))
res_table$threshold<- threshold

plotCounts(dds, gene="CHD4", intgroup = "exp" )

#Volcano plot
df_results<- data.frame(res_table)
ggplot(df_results) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold))  +
  ggtitle('INSERT TITLE HERE') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) 
ggsave("Volcano_plot.tiff", plot = last_plot(), device = "tiff")
dev.off()

#Distinguish which genes are up and which are down
res_table_up<- subset(res_table, threshold == TRUE & res_table$log2FoldChange > 0)
sig_genes_up<- row.names(res_table_up)

write.table(res_table_up, file = "INSERT FILE NAME HERE", row.names = T, quote = F)
write.table(sig_genes_up, file = "INSERT FILE NAME HERE", row.names = F, quote = F)

res_table_down<- subset(res_table, threshold == TRUE & res_table$log2FoldChange < 0)
sig_genes_down<- row.names(res_table_down)

write.table(res_table_down, file = "INSERT FILE NAME HERE", row.names = T, quote = F)
write.table(sig_genes_down, file = "INSERT FILE NAME HERE", row.names = F, quote = F)

#Making gene tables and gene lists for sig genes that fit criteria (no matter up or down)
res_table_sig<- subset(res_table, threshold == TRUE)
#write.table(res_table_sig, file = "INSERT FILE NAME.txt", row.names = T, quote = F)
sig_genes<- unique(row.names(res_table_sig))
#write.table(sig_genes, file = "INSERT FILE NAME.txt", row.names = F, quote = F)

#Extract normalized counts for significant genes
norm_counts_sig<- combined_normalized_counts[sig_genes, ]

#Annotate heatmap
annotation<- metadata[,"exp", drop=F]
#Set a color palette
heat.colors <- brewer.pal(6,"YlOrRd")
#Run pheatmap
pheatmap(norm_counts_sig, color = heat.colors, cluster_rows = T, show_rownames=F,
         annotation= annotation, border_color=NA, fontsize = 10, scale="row",
         fontsize_row = 10, height=20, show_colnames = FALSE,filename = "heatmap.pdf")
dev.off()

gprofiler_results<- gprofiler(query = sig_genes, 
                                   organism = "hsapiens",
                                   ordered_query = F, 
                                   exclude_iea = F, 
                                   max_p_value = 0.05, 
                                   max_set_size = 0,
                                   correction_method = "fdr",
                                   hier_filtering = "none", 
                                   domain_size = "annotated",
                                   custom_bg = "")

write.table(gprofiler_results,"GProfiler_results.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms <- cbind(gprofiler_results$term.id, gprofiler_results$p.value)
GOs <- allterms[grep('GO:', allterms),]
write.table(GOs, "GO_terms.txt", sep="\t", quote=F, row.names=F, col.names=F)


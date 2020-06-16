
# DESeq2 in R for beginners

More detailed information can be found here: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# Install and load Packages
```
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
```

# What you need to start
1. Counts file

## Of Note- Input data
DESEQ2 is an R package that takes an annotated counts file and performs differential expression analysis. There are two main inputs for DESEQ2, and it is your choice as to which you might use: 
1) an annotated counts matrix file (developed using a program like FeatureCounts), which requires bam files as input (alignment to STAR). This is what we will be using in this markdown.
2) an annotated counts file that was generated using transcript abundance quantifiers (Salmon, Sailfish, kallisto) from fastq files. This method is trending nowadays and is useful for quantifying alternate transcripts.
The commands are almost the same for the two files, but I will mark the difference in the script.

## Of Note- Normalized versus un-normalized reads
Also note that we are using un-normalized counts here (not RPKM), which are not normalized to gene LENGTH. DESEQ2 is designed to compare Gene A in Sample 1 versus Sample 2. Therefore, the length of the gene is irrelevant because we are comparing the same gene in two different samples. DESEQ2 does have a count normalization step within the analysis , so you could analyze data from separate batches (I would add a batch column to my metadata to test for potential batch effect) If you wanted to compare Gene A abundance to Gene B abundance in Sample 1, you would do this using RPKM, which DESEQ does have a command to calculate.

2. Metadata file

A metadata file contains all of the relevant information to an experiment. I have provided an example [here](https://github.com/drlaurenwasson/R_tutorials/blob/master/DESeq2/metadata.csv). 

My example experiment compares different mutant iPSC lines at Day 0 of differentiation. In my example I have the sample name, a description of the sample (to me, its the same, but to others it might not be), the cell line, the genotype, and what I call "exp", which is the variable that I choose to categorize the data by for analysis. "Exp" is a column in which I have no special characters (such as "/", which messes with the analysis but I use to describe the genotype), and is short enough to use as a label for things like PCA and heatmaps and other graphs. Other columns of metadata that might be useful include: time point, drug concentration, mouse strain background (if you have different strains), batch, date of sequencing, the person who did the experiment (if you're combining data).

# Data Analysis

First, load in the counts file.
```
counts_file<- read.table("CHD7_data_1.txt")
```
You may want to only include protein coding genes in your analysis. 

```
#remove non-protein coding genes from analysis
protein_genes<- read.table("protein_genes.txt")
protein_genes<- as.character(protein_genes$V1)

counts_protein_coding <- counts[rownames(counts) %in% protein_genes,]
```
Next, read in your metadata file
```
metadata<- read.csv("metadata.csv")
```

The column names in the counts file MUST match the row names in the metadata file. 
```
all(names(counts_protein_coding) %in% rownames(metadata))
```
This must return "TRUE" to continue. If not, I find it's easiest to alter the rownames in the metadata file to match the column names of the counts file. Resave the metadata file and re-upload it into R.

## Generate the dds (the backbone file of the analysis)
Here I choose to perform the analysis using the groupings specified in the "exp" column of my metadata file. You can specify primary and secondary groupings if you wish (time point, genotype). Information for that is found in the detailed vignette.
This is the command for a counts matrix file ("DESeqDataSetFromMatrix"). The commands would be different for txtimport from Sailfish or Salmon, and can be found on the page referenced above. 

```
#DDS by "exp"
dds<- DESeqDataSetFromMatrix(countData = counts_protein_coding, colData = metadata, design = ~ exp )
```
## Normalize counts/filter
Here you can do things like only keep rows that have more than 10 total reads (between all the genotypes) to filter out very lowly expressed transcripts.

```
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
```

You can also normalize for read depth here

```
dds<- estimateSizeFactors(dds)
sizeFactors(dds)
combined_normalized_counts<- counts(dds, normalized =TRUE)
write.table(combined_normalized_counts, file= "dds_combined_normalized_counts_proteingenes.txt", sep = "\t", quote = FALSE)
```
The last thing you need to do is assign a "reference" sample to do the comparisons. "By default, R will choose a reference level for factors based on alphabetical order. Then, if you never tell the DESeq2 functions which level you want to compare against (e.g. which level represents the control group), the comparisons will be based on the alphabetical order of the levels." To set your reference, type in the reference from your metadata file (In my example, it's "wt").
```
dds$exp<- relevel(dds$exp, ref = "wt" )
```

## Transform counts for PCA analysis
Here we are log-transforming the dds file above to use in a PCA analysis. The PCA will plot the log-transformed data, and group by a variable that you provide (here it's "exp"). It will save the PCA in your working directory as a tiff file named "PCA_analysis.tiff".

```
rld <- rlog(dds, blind = TRUE)

#Plot PCA
tiff(filename = "PCA_analysis.tiff",
     width = 480 ,height = 480, units = "px", pointsize = 12)
plotPCA(rld, intgroup="exp")
dev.off()
```

Now, to make the PCA fancy! Here, I have used ggplot2 to plot the same graph with a few modifications

- chosen my own color scheme (a variable called my_colors)
- chosen which components to actually graph (by default, a PCA graph will graph the first two principal components, but sometimes you might want to plot different components). This is represented by aes(PC1,PC2). I also chose here which variable to color by (exp).
- changed the size of each dot (geom_point)
- added X and Y labels to show the % variance
```
pcaData <- plotPCA(rld, intgroup="exp", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
mycolors<- c("#209324","#ef0000","#1166d6","#e27106", "#e810aa", "#000000")
pcaplot<- ggplot(pcaData, aes(PC1,PC2, color = exp), scale_color_manual(values = mycolors)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
pcaplot+ scale_color_manual(values = mycolors)
```
You can then export this PCA and overwrite the first one, if you'd like.

## Correlation test
This is a quick R function that runs a correlation test to say how "alike" is each sample to every other sample. Values closest to 1 are highly correlative. Simetimes, if I get a funky-looking PCA, this can help me determine the weird sample in the dataset to remove it.

```
cortest<- cor(counts_protein_coding, use = "all.obs", method = "pearson")
cortest
write.csv(cortest, file= "corrtest.csv", row.names= TRUE, col.names = TRUE)
```
## Run DESEQ
```
dds <- DESeq(dds)
plotDispEsts(dds)
```

## Generate some heatmaps of most DE genes (1000 most DE genes, and 250 most DEGs)
```
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
```

## Another way to do PCA, with the added bonus of taking the genes that make up the principle components
```
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
```
# Generate Results tables
Here we generate results tables and actually do the comparison between wild type and mutant. In this example we are comparing "CHD7_het" to "wt"

```
res_table<- results(dds, contrast = c("group", "CHD7_het", "wt"))
summary(res_table)

```
## Determine significant genes by fold change and p-value
This will subset your results table using an adjusted p-value and a log2 fold change cutoff (0.58 translates to 1.5 fold). The absolute log2 fold change is used so that it will take results that are differential in both directions. The last command will add a column to your table called "threshold", which is "True" if the gene is significant, and "False" if the gene is not significantly different.
```
#Set thresholds
padj.cutoff <- 0.01
lfc.cutoff <- 0.58

threshold<- res_table$padj < padj.cutoff & abs(res_table$log2FoldChange) > lfc.cutoff
length(which(threshold == TRUE))
res_table$threshold<- threshold
```

## Plot expression for a single gene
```
plotCounts(dds, gene="CHD4", intgroup = "exp" )
```
## Genetate a volcano plot of all genes, with significant genes colored
This uses the R package ggplot2 to generate a volcano plot of our data
```
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
```
## Distinguish genes that are up versus down
```
#Distinguish which genes are up and which are down
res_table_up<- subset(res_table, threshold == TRUE & res_table$log2FoldChange > 0)
sig_genes_up<- row.names(res_table_up)

write.table(res_table_up, file = "INSERT FILE NAME HERE", row.names = T, quote = F)
write.table(sig_genes_up, file = "INSERT FILE NAME HERE", row.names = F, quote = F)

res_table_down<- subset(res_table, threshold == TRUE & res_table$log2FoldChange < 0)
sig_genes_down<- row.names(res_table_down)

write.table(res_table_down, file = "INSERT FILE NAME HERE", row.names = T, quote = F)
write.table(sig_genes_down, file = "INSERT FILE NAME HERE", row.names = F, quote = F)
```

## Combine the two to make a total DEG list
```
#Making gene tables and gene lists for sig genes that fit criteria (no matter up or down)
res_table_sig<- subset(res_table, threshold == TRUE)
#write.table(res_table_sig, file = "INSERT FILE NAME.txt", row.names = T, quote = F)
sig_genes<- unique(row.names(res_table_sig))
#write.table(sig_genes, file = "INSERT FILE NAME.txt", row.names = F, quote = F)
```

## Heatmap
```
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
```


#biocLite("VennDiagram")
library(VennDiagram)
A.cds.sig.new <- c(sig1D4_combined_up, sig1D4_combined_down)
write.table(A.cds.sig.new, file = "A.cds.sig.new_15fold.txt", row.names = F, quote = F, col.names=F)
A.cds.sig.up.new <- sig1D4_combined_up
A.cds.sig.down.new <- sig1D4_combined_down
C.cds.sig.new <- c(sigG1_combined_up, sigG1_combined_down)
write.table(C.cds.sig.new, file = "C.cds.sig.new_15fold.txt", row.names = F, quote = F, col.names=F)
C.cds.sig.up.new <- sigG1_combined_up
C.cds.sig.down.new <- sigG1_combined_down
D.cds.sig <- c(sigCHD414_combined_up, sigCHD414_combined_down)
D.cds.sig.up <- sigCHD414_combined_up
D.cds.sig.down <- sigCHD414_combined_down
E.cds.sig <- c(sigCHD427_combined_up, sigCHD427_combined_down)
E.cds.sig.up <- sigCHD427_combined_up
E.cds.sig.down <- sigCHD427_combined_down

AC.cds.sig.up<- A.cds.sig.up[A.cds.sig.up %in% C.cds.sig.up]
write.table(AC.cds.sig.up, file = "2019_01_11_1D41G1_genes.sig_UP_15fold.txt", row.names = F, quote = F, col.names=F)
AC.cds.sig.down<- A.cds.sig.down[A.cds.sig.down %in% C.cds.sig.down]
write.table(AC.cds.sig.down, file = "2019_01_11_1D41G1_genes.sig_DOWN_15fold.txt", row.names = F, quote = F, col.names=F)
AC.cds.sig<- c(AC.cds.sig.up,AC.cds.sig.down)
write.table(AC.cds.sig, file = "2019_01_11_1D41G1_genes.sig_15fold.txt", row.names = F, quote = F, col.names=F)

AC.cds.sig.up.new<- A.cds.sig.up.new[A.cds.sig.up.new %in% C.cds.sig.up.new]
#write.table(AC.cds.sig.up.new, file = "2019_01_11_1D41G1_genes.sig_UP_15fold.txt", row.names = F, quote = F, col.names=F)
AC.cds.sig.down.new<- A.cds.sig.down.new[A.cds.sig.down.new %in% C.cds.sig.down.new]
#write.table(AC.cds.sig.down, file = "2019_01_11_1D41G1_genes.sig_DOWN_15fold.txt", row.names = F, quote = F, col.names=F)

#Compare the two CHD4 hets to each other
DE.cds.sig.up<- D.cds.sig.up[D.cds.sig.up %in% E.cds.sig.up]
write.table(DE.cds.sig.up, file = "CHD414_CHD427_UP_15fold.txt", row.names = F, quote = F, col.names=F)
DE.cds.sig.down<- D.cds.sig.down[D.cds.sig.down %in% E.cds.sig.down]
write.table(DE.cds.sig.down, file = "CHD414_CHD427_DOWN_15fold.txt", row.names = F, quote = F, col.names=F)
DE.cds.sig<- c(DE.cds.sig.up,DE.cds.sig.down)
write.table(AD.cds.sig, file = "CHD414_CHD7_genes.sig_15fold.txt", row.names = F, quote = F, col.names=F)

#Compare the CHD4 het to CHD7 (old and new)
AD.cds.sig.up<- A.cds.sig.up[A.cds.sig.up %in% D.cds.sig.up]
write.table(AD.cds.sig.up, file = "CHD414_CHD7_UP_15fold.txt", row.names = F, quote = F, col.names=F)
AD.cds.sig.down<- A.cds.sig.down[A.cds.sig.down %in% D.cds.sig.down]
write.table(AD.cds.sig.down, file = "CHD414_CHD7_DOWN_15fold.txt", row.names = F, quote = F, col.names=F)
AD.cds.sig<- c(AD.cds.sig.up,AD.cds.sig.down)
write.table(AD.cds.sig, file = "CHD414_CHD7_genes.sig_15fold.txt", row.names = F, quote = F, col.names=F)

AD.cds.sig.up.new<- A.cds.sig.up.new[A.cds.sig.up.new %in% D.cds.sig.up]
write.table(AD.cds.sig.up.new, file = "CHD414_CHD7NEW_UP_15fold.txt", row.names = F, quote = F, col.names=F)
AD.cds.sig.down.new<- A.cds.sig.down.new[A.cds.sig.down.new %in% D.cds.sig.down]
write.table(AD.cds.sig.down, file = "CHD414_CHD7NEW_genes.sig_DOWN_15fold.txt", row.names = F, quote = F, col.names=F)
AD.cds.sig.new<- c(AD.cds.sig.up.new,AD.cds.sig.down.new)
write.table(AD.cds.sig, file = "CHD414_CHD7NEW_genes.sig_15fold.txt", row.names = F, quote = F, col.names=F)
AD.cds.sig[!AD.cds.sig %in% AD.cds.sig.new]

#AE.cds.sig.up<- A.cds.sig.up[A.cds.sig.up %in% E.cds.sig.up]
#write.table(AE.cds.sig.up, file = "2019_01_11_1D41G1_genes.sig_UP_15fold.txt", row.names = F, quote = F, col.names=F)
#AE.cds.sig.down<- A.cds.sig.down[A.cds.sig.down %in% D.cds.sig.down]
#write.table(AE.cds.sig.down, file = "2019_01_11_1D41G1_genes.sig_DOWN_15fold.txt", row.names = F, quote = F, col.names=F)
#AE.cds.sig<- c(AE.cds.sig.up,AE.cds.sig.down)
#write.table(AE.cds.sig, file = "2019_01_11_1D41G1_genes.sig_15fold.txt", row.names = F, quote = F, col.names=F)

AE.cds.sig.up<- A.cds.sig.up.new[A.cds.sig.up.new %in% E.cds.sig.up]
write.table(AE.cds.sig.up, file = "CHD427_CHD7NEW_genes.sig_UP_15fold.txt", row.names = F, quote = F, col.names=F)
AE.cds.sig.down<- A.cds.sig.down.new[A.cds.sig.down.new %in% E.cds.sig.down]
write.table(AE.cds.sig.down, file = "2019_01_11_1D41G1_genes.sig_DOWN_15fold.txt", row.names = F, quote = F, col.names=F)
AE.cds.sig.new<- c(AE.cds.sig.up,AE.cds.sig.down)
write.table(AE.cds.sig, file = "2019_01_11_1D41G1_genes.sig_15fold.txt", row.names = F, quote = F, col.names=F)
AE.cds.sig[!AE.cds.sig %in% AE.cds.sig.new]

ADE.cds.sig.up<- AD.cds.sig.up.new[AD.cds.sig.up.new %in% E.cds.sig.up]
ADE.cds.sig.down<- AD.cds.sig.down.new[AD.cds.sig.down.new %in% E.cds.sig.down]

##STOP HERE
#Venn Diagrams
pdf("2019_01_10_CHD7_D3U_het_UP.pdf")
venn.plot <- venn.diagram(list(A.cds.sig.up, B.cds.sig.up), NULL, fill=c("blue", "red"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("CHD7 Q1599X/+ UP", "CHD7 K1597X/+ UP"))
grid.draw(venn.plot)
dev.off()

pdf("2019_01_10_CHD7_D3U_het_DOWN.pdf")
venn.plot <- venn.diagram(list(A.cds.sig.down, B.cds.sig.down), NULL, fill=c("blue", "red"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("CHD7 Q1599X/+ DOWN", "CHD7 K1597X/+ DOWN"))
grid.draw(venn.plot)
dev.off()

pdf("2019_01_10_CHD7_D3U_1D41G1_UP.pdf")
venn.plot <- venn.diagram(list(A.cds.sig.up, C.cds.sig.up), NULL, fill=c("blue", "orange"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("CHD7 Q1599X/+ UP", "CHD7 Q1599X/Q1599X UP"))
grid.draw(venn.plot)
dev.off()

pdf("2019_01_10_CHD7_D3U_1D41G1_DOWN.pdf")
venn.plot <- venn.diagram(list(A.cds.sig.down, C.cds.sig.down), NULL, fill=c("blue", "orange"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("CHD7 Q1599X/+ DOWN", "CHD7 Q1599X/Q1599X DOWN"))
grid.draw(venn.plot)
dev.off()

pdf("2019_01_10_CHD7_D3U_1G61G1_UP.pdf")
venn.plot <- venn.diagram(list(B.cds.sig.up, C.cds.sig.up), NULL, fill=c("red", "orange"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("CHD7 K1597X/+ UP", "CHD7 Q1599X/Q1599X UP"))
grid.draw(venn.plot)
dev.off()

pdf("2019_01_10_CHD7_D3U_1G61G1_DOWN.pdf")
venn.plot <- venn.diagram(list(B.cds.sig.down, C.cds.sig.down), NULL, fill=c("red", "orange"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("CHD7 K1597X/+ DOWN", "CHD7 Q1599X/Q1599X DOWN"))
grid.draw(venn.plot)
dev.off()


pdf("2019_01_10_CHD7_D3U_het_DOWN.pdf")
venn.plot <- venn.diagram(list(A.cds.sig.down, B.cds.sig.down), NULL, fill=c("blue", "red"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("CHD7 Q1599X/+ UP", "CHD7 K1597X/+ UP"))
grid.draw(venn.plot)
dev.off()


pdf("2019_01_11_LOF.pdf")
venn.plot<- draw.triple.venn(1905, 2356, 1694, 886, 935, 549, 401, category =
                               c("", "", ""), rotation = 1, reverse = FALSE, euler.d =
                               TRUE, scaled = TRUE, lwd = rep(2, 3), lty =
                               rep("solid", 3), col = rep("black", 3), fill = c("blue", "red", "orange"),
                             alpha = rep(0.5, 3), label.col = rep("black", 7), cex
                             = 5, fontface = rep("plain", 7), fontfamily =
                               rep("serif", 7), cat.pos = c(-40, 40, 180), cat.dist =
                               c(0.05, 0.05, 0.025), cat.col = rep("black", 3),
                             cat.cex = 1, cat.fontface = 4,
                             cat.fontfamily = rep("serif", 3), cat.just =
                               list(c(0.5, 1), c(0.5, 1), c(0.5, 0)), cat.default.pos
                             = "outer", cat.prompts = FALSE, rotation.degree = 0,
                             rotation.centre = c(0.5, 0.5), ind = TRUE, sep.dist =
                               0.05, offset = 0, cex.prop = NULL, print.mode = "raw",
                             sigdigs = 3, direct.area = FALSE)
grid.draw(venn.plot)
dev.off()

#Turn the results tables into data frames
res_table_1D4_combined_df<- as.data.frame(res_table_1D4_combined)
res_table_1G6_combined_df<- as.data.frame(res_table_1G6_combined)
res_table_1G1_combined_df<- as.data.frame(res_table_G1_combined)
res_table_CHD414_combined_df<- as.data.frame(res_table_CHD414_combined)
#res_table_CHD427_combined_df<- as.data.frame(res_table_CHD427_combined)

#scatter plots between the lines
pdf("2019_01_11_1D4vs1G6_scatter.pdf")
plot(res_table_1D4_combined_df$log2FoldChange, res_table_1G6_combined_df$log2FoldChange, main = "CHD7 Q1599X/+ vs CHD7 K1597X/+", xlab= "CHD7 Q1599X/+ log2foldchange", ylab = "CHD7 K1597X/+ log2foldchange") 
abline(fit<-lm(res_table_1D4_combined_df$log2FoldChange~res_table_1G6_combined_df$log2FoldChange), col = "red")
legend("bottomright", bty = "n", legend=paste("R2=",format(summary(fit)$adj.r.squared, digits=4)))
dev.off()

pdf("2019_01_11_1D4vs1G1_scatter.pdf")
plot(res_table_1D4_combined_df$log2FoldChange, res_table_1G1_combined_df$log2FoldChange, main = "CHD7 Q1599X/+ vs CHD7 Q1599X/Q1599X", xlab= "CHD7 Q1599X/+ log2foldchange", ylab = "CHD7 Q1599X/Q1599X log2foldchange") 
abline(fit<-lm(res_table_1D4_combined_df$log2FoldChange~res_table_1G1_combined_df$log2FoldChange), col = "red")
legend("bottomright", bty = "n", legend=paste("R2=",format(summary(fit)$adj.r.squared, digits=4)))
dev.off()

pdf("2019_01_11_1G6vs1G1_scatter.pdf")
plot(res_table_1G6_combined_df$log2FoldChange, res_table_1G1_combined_df$log2FoldChange, main = "CHD7 K1597X/+ vs CHD7 Q1599X/Q1599X", xlab= "CHD7 K1597X/+ log2foldchange", ylab = "CHD7 Q1599X/Q1599X log2foldchange") 
abline(fit<-lm(res_table_1G6_combined_df$log2FoldChange~res_table_1G1_combined_df$log2FoldChange), col = "red")
legend("bottomright", bty = "n", legend=paste("R2=",format(summary(fit)$adj.r.squared, digits=4)))
dev.off()


library(gProfileR)

# Running gprofiler to identify enriched processes among significant genes
gprofiler_results_1D4 <- gprofiler(query = A.cds.sig, 
                                   organism = "hsapiens",
                                   ordered_query = F, 
                                   exclude_iea = F, 
                                   max_p_value = 0.05, 
                                   max_set_size = 0,
                                   correction_method = "fdr",
                                   hier_filtering = "none", 
                                   domain_size = "annotated",
                                   custom_bg = "")

write.table(gprofiler_results_1D4,"2019_01_11_gProfileR_1D4.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_1D4 <- cbind(gprofiler_results_1D4$term.id, gprofiler_results_1D4$p.value)
GOs_1D4 <- allterms_1D4[grep('GO:', allterms_1D4),]
write.table(GOs_1D4, "2019_01_11_GOS_1D4.txt", sep="\t", quote=F, row.names=F, col.names=F)

gprofiler_results_1D4_unique <- gprofiler(query = A.unique.cds.sig, 
                                          organism = "hsapiens",
                                          ordered_query = F, 
                                          exclude_iea = F, 
                                          max_p_value = 0.05, 
                                          max_set_size = 0,
                                          correction_method = "fdr",
                                          hier_filtering = "none", 
                                          domain_size = "annotated",
                                          custom_bg = "")

write.table(gprofiler_results_1D4_unique,"2019_01_11_gProfileR_1D4_unique.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_1D4_unique <- cbind(gprofiler_results_1D4_unique$term.id, gprofiler_results_1D4_unique$p.value)
GOs_1D4_unique <- allterms_1D4_unique[grep('GO:', allterms_1D4_unique),]
write.table(GOs_1D4_unique, "2019_01_11_GOS_1D4_unique.txt", sep="\t", quote=F, row.names=F, col.names=F)

gprofiler_results_1G6 <- gprofiler(query = B.cds.sig, 
                                   organism = "hsapiens",
                                   ordered_query = F, 
                                   exclude_iea = F, 
                                   max_p_value = 0.05, 
                                   max_set_size = 0,
                                   correction_method = "fdr",
                                   hier_filtering = "none", 
                                   domain_size = "annotated",
                                   custom_bg = "")

write.table(gprofiler_results_1G6,"2019_01_11_gProfileR_1G6.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_1G6 <- cbind(gprofiler_results_1G6$term.id, gprofiler_results_1G6$p.value)
GOs_1G6 <- allterms_1G6[grep('GO:', allterms_1G6),]
write.table(GOs_1G6, "2019_01_11_GOS_1G6.txt", sep="\t", quote=F, row.names=F, col.names=F)

# Running gprofiler to identify enriched processes among significant genes
gprofiler_results_1G6_unique <- gprofiler(query = B.unique.cds.sig, 
                                          organism = "hsapiens",
                                          ordered_query = F, 
                                          exclude_iea = F, 
                                          max_p_value = 0.05, 
                                          max_set_size = 0,
                                          correction_method = "fdr",
                                          hier_filtering = "none", 
                                          domain_size = "annotated",
                                          custom_bg = "")

write.table(gprofiler_results_1G6_unique,"2019_01_11_gProfileR_1G6_unique.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_1G6_unique <- cbind(gprofiler_results_1G6_unique$term.id, gprofiler_results_1G6_unique$p.value)
GOs_1G6_unique <- allterms_1G6_unique[grep('GO:', allterms_1G6_unique),]
write.table(GOs_1G6_unique, "2019_01_11_GOS_1G6_unique.txt", sep="\t", quote=F, row.names=F, col.names=F)

gprofiler_results_lof <- gprofiler(query = lof.cds.sig, 
                                   organism = "hsapiens",
                                   ordered_query = F, 
                                   exclude_iea = F, 
                                   max_p_value = 0.05, 
                                   max_set_size = 0,
                                   correction_method = "fdr",
                                   hier_filtering = "none", 
                                   domain_size = "annotated",
                                   custom_bg = "")

write.table(gprofiler_results_lof,"2019_01_11_gProfileR_LOF.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_lof <- cbind(gprofiler_results_lof$term.id, gprofiler_results_lof$p.value)
GOs_lof <- allterms_lof[grep('GO:', allterms_lof),]
write.table(GOs_lof, "2019_01_11_GOS_LOF.txt", sep="\t", quote=F, row.names=F, col.names=F)

gprofiler_results_abcde <- gprofiler(query = ABCDE.cds.sig, 
                                     organism = "hsapiens",
                                     ordered_query = F, 
                                     exclude_iea = F, 
                                     max_p_value = 0.05, 
                                     max_set_size = 0,
                                     correction_method = "fdr",
                                     hier_filtering = "none", 
                                     domain_size = "annotated",
                                     custom_bg = "")

write.table(gprofiler_results_abde,"2019_01_11_gProfileR_LOF.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_lof <- cbind(gprofiler_results_lof$term.id, gprofiler_results_lof$p.value)
GOs_lof <- allterms_lof[grep('GO:', allterms_lof),]
write.table(GOs_lof, "2019_01_11_GOS_LOF.txt", sep="\t", quote=F, row.names=F, col.names=F)




#STOP HERE _ MISSENSE

pdf("2018_08_02_CHD7_D3U_combined_bcbio_hetvsCHD414_UP.pdf")
venn.plot <- venn.diagram(list(A.cds.sig.up, D.cds.sig.up), NULL, fill=c("purple", "green"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("CHD7 +/- UP", "CHD7 CHD414/CHD414 UP"))
grid.draw(venn.plot)
dev.off()

pdf("2018_08_02_CHD7_D3U_combined_bcbio_hetsvsCHD414_DOWN.pdf")
venn.plot <- venn.diagram(list(A.cds.sig.down, D.cds.sig.down), NULL, fill=c("purple", "green"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("CHD7 +/- DOWN", "CHD7 CHD414/CHD414 DOWN"))
grid.draw(venn.plot)
dev.off()

pdf("2018_08_02_CHD7_D3U_combined_hetsvsCHD414.pdf")
venn.plot <- draw.pairwise.venn(area1 = 3404, area2 = 4069, cross.area=1976, fill=c("purple", "green"), alpha=c(0.5,0.5), cex = 5, cat.fontface=4, category=c("", ""), cat.pos = c(200, 160))
grid.draw(venn.plot)
dev.off()
phyper(1976,3404,(57373-3404),4069, lower.tail = F, log.p = TRUE)


pdf("2018_08_02_CHD7_D3U_combined_bcbio_1G1vsCHD414_UP.pdf")
venn.plot <- venn.diagram(list(C.cds.sig.up, D.cds.sig.up), NULL, fill=c("orange", "green"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("CHD7 -/- UP", "CHD7 CHD414/CHD414 UP"))
grid.draw(venn.plot)
dev.off()

pdf("2018_08_02_CHD7_D3U_combined_bcbio_1G1vsCHD414_DOWN.pdf")
venn.plot <- venn.diagram(list(C.cds.sig.down, D.cds.sig.down), NULL, fill=c("orange", "green"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("CHD7 -/- DOWN", "CHD7 CHD414/CHD414 DOWN"))
grid.draw(venn.plot)
dev.off()

pdf("2018_08_02_CHD7_D3U_combined_1G1vsCHD414.pdf")
venn.plot <- draw.pairwise.venn(area1 = 1790, area2 = 4069, cross.area=752, fill=c("orange", "green"), alpha=c(0.5,0.5), cex = 5, cat.fontface=4, category=c("", ""), cat.pos = c(200, 160))
grid.draw(venn.plot)
dev.off()
phyper(752,1790,(57373-1790),4069, lower.tail = F, log.p = TRUE)


pdf("2018_08_02_CHD7_D3U_combined_bcbio_hetvsCHD427_UP.pdf")
venn.plot <- venn.diagram(list(A.cds.sig.up, E.cds.sig.up), NULL, fill=c("purple", "pink"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("CHD7 +/- UP", "CHD7 CHD427/CHD427 UP"))
grid.draw(venn.plot)
dev.off()

pdf("2018_08_02_CHD7_D3U_combined_bcbio_hetsvsCHD427_DOWN.pdf")
venn.plot <- venn.diagram(list(A.cds.sig.down, E.cds.sig.down), NULL, fill=c("purple", "pink"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("CHD7 +/- DOWN", "CHD7 CHD427/CHD427 DOWN"))
grid.draw(venn.plot)
dev.off()

pdf("2018_08_02_CHD7_D3U_combined_hetsvsCHD427.pdf")
venn.plot <- draw.pairwise.venn(area1 = 3404, area2 = 1079, cross.area=684, fill=c("purple", "pink"), alpha=c(0.5,0.5), cex = 5, cat.fontface=4, category=c("", ""), cat.pos = c(200, 160))
grid.draw(venn.plot)
dev.off()
phyper(684,3404,(57373-1304),1079, lower.tail = F, log.p = TRUE)


pdf("2018_08_02_CHD7_D3U_combined_bcbio_1G1vsCHD427_UP.pdf")
venn.plot <- venn.diagram(list(C.cds.sig.up, E.cds.sig.up), NULL, fill=c("orange", "pink"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("CHD7 -/- UP", "CHD7 CHD427/CHD427 UP"))
grid.draw(venn.plot)
dev.off()

pdf("2018_08_02_CHD7_D3U_combined_bcbio_1G1vsCHD427_DOWN.pdf")
venn.plot <- venn.diagram(list(C.cds.sig.down, E.cds.sig.down), NULL, fill=c("orange", "pink"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("CHD7 -/- DOWN", "CHD7 CHD427/CHD427 DOWN"))
grid.draw(venn.plot)
dev.off()

pdf("2018_08_02_CHD7_D3U_combined_1G1vsCHD427.pdf")
venn.plot <- draw.pairwise.venn(area1 = 1790, area2 = 1079, cross.area=450, fill=c("orange", "pink"), alpha=c(0.5,0.5), cex = 5, cat.fontface=4, category=c("", ""), cat.pos = c(200, 160))
grid.draw(venn.plot)
dev.off()
phyper(450,1790,(57373-1790),1079, lower.tail = F, log.p = TRUE)

pdf("2018_08_02_CHD7_D3U_combined_bcbio_1G1vsCHD427_UP.pdf")
venn.plot <- venn.diagram(list(C.cds.sig.up, E.cds.sig.up), NULL, fill=c("orange", "pink"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("CHD7 -/- UP", "CHD7 CHD427/CHD427 UP"))
grid.draw(venn.plot)
dev.off()

pdf("2018_08_02_CHD7_D3U_combined_bcbio_1G1vsCHD427_DOWN.pdf")
venn.plot <- venn.diagram(list(C.cds.sig.down, E.cds.sig.down), NULL, fill=c("orange", "pink"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("CHD7 -/- DOWN", "CHD7 CHD427/CHD427 DOWN"))
grid.draw(venn.plot)
dev.off()

pdf("2018_08_02_CHD7_D3U_combined_1G1vsCHD427.pdf")
venn.plot <- draw.pairwise.venn(area1 = 1790, area2 = 1079, cross.area=450, fill=c("orange", "pink"), alpha=c(0.5,0.5), cex = 5, cat.fontface=4, category=c("", ""), cat.pos = c(200, 160))
grid.draw(venn.plot)
dev.off()
phyper(450,1790,(57373-1790),1079, lower.tail = F, log.p = TRUE)

pdf("2018_08_02_CHD7_D3U_combined_bcbio_CHD414vsCHD427_UP.pdf")
venn.plot <- venn.diagram(list(D.cds.sig.up, E.cds.sig.up), NULL, fill=c("green", "pink"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("", ""))
grid.draw(venn.plot)
dev.off()

pdf("2018_08_02_CHD7_D3U_combined_bcbio_CHD414vsCHD427_DOWN.pdf")
venn.plot <- venn.diagram(list(D.cds.sig.down, E.cds.sig.down), NULL, fill=c("green", "pink"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("", ""))
grid.draw(venn.plot)
dev.off()

pdf("2018_08_02_CHD7_D3U_combined_CHD414vsCHD427.pdf")
venn.plot <- draw.pairwise.venn(area1 = 4069, area2 = 1079, cross.area=736, fill=c("green", "pink"), alpha=c(0.5,0.5), cex = 5, cat.fontface=4, category=c("", ""), cat.pos = c(200, 160))
grid.draw(venn.plot)
dev.off()
phyper(736,4069,(57373-4069),1079, lower.tail = F, log.p = TRUE)

###heatmap data

library (ggplot2)
install.packages("reshape")
library(reshape)
library(RColorBrewer)
#------------------
# CREATE DATA FRAME
#------------------
df.team_data <- expand.grid(teams = c("CHD7+/-", "CHD7-/-", "CHD414/CHD414", "CHD427/CHD427")
                            ,metrics = c("CHD7+/-", "CHD7-/-", "CHD414/CHD414", "CHD427/CHD427")
)

# add variable: performance
percentages_homo<- c(NA, 30.55, 58.0, 20.09, 59, NA, 42, 25.1, 48.5, 18.48, NA, 18.08, 63.3, 41.7, 68.2, NA )
df.team_data$performance <- percentages_homo
df.team_data$performance<- as.numeric(df.team_data$performance)

df.team_data$teams<- factor(df.team_data$teams, levels = c("CHD427/CHD427", "CHD414/CHD414", "CHD7-/-", "CHD7+/-"))

hm.palette <- colorRampPalette(brewer.pal(9, 'Blues'), space='Lab')


#---------------------------
# PLOT: heatmap
# - here, we use geom_tile()
#---------------------------
tiff("CHD7_Differential_RNA-seq.tiff")
ggplot(data = df.team_data, aes(x = metrics, y = teams, fill = performance)) +
  geom_tile()+
  coord_equal()+
  scale_fill_gradientn(colours = hm.palette(100))+
  scale_x_discrete(position = "top")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+
  labs(fill = "Percent Overlap")
dev.off()












#STOP HERE

areahet<- 3404
areahomo<- 1790
areaCHD414<- 4069
areaCHD427<- 1079
areahet_homo<- 1040
areahet_CHD414<- 1976
areahomo_D2355<- 752
areahet_CHD427<- 684
areahomo_CHD427<- 450
areaCHD414_CHD427<- 736
area_het_CHD414_CHD427<- 556

pdf("2018_08_03_hets_CHD414_CHD427.pdf")
venn.plot<- draw.triple.venn(areahet, areaCHD414, areaCHD427, areahet_CHD414, areaCHD414_CHD427, areahet_CHD427, area_het_CHD414_CHD427, category =
                               c("", "", ""), rotation = 1, reverse = FALSE, euler.d =
                               TRUE, scaled = TRUE, lwd = rep(2, 3), lty =
                               rep("solid", 3), col = rep("black", 3), fill = c("purple", "green", "pink"),
                             alpha = rep(0.5, 3), label.col = rep("black", 7), cex
                             = 5, fontface = rep("plain", 7), fontfamily =
                               rep("serif", 7), cat.pos = c(-40, 40, 180), cat.dist =
                               c(0.05, 0.05, 0.025), cat.col = rep("black", 3),
                             cat.cex = 1, cat.fontface = 4,
                             cat.fontfamily = rep("serif", 3), cat.just =
                               list(c(0.5, 1), c(0.5, 1), c(0.5, 0)), cat.default.pos
                             = "outer", cat.prompts = FALSE, rotation.degree = 0,
                             rotation.centre = c(0.5, 0.5), ind = TRUE, sep.dist =
                               0.05, offset = 0, cex.prop = NULL, print.mode = "raw",
                             sigdigs = 3, direct.area = FALSE)
grid.draw(venn.plot)
dev.off()

#Turn the results tables into data frames
res_table_hets_combined_df<- as.data.frame(res_table_het_combined_test)
res_table_G1_combined_df<- as.data.frame(res_table_G1_combined)
res_table_CHD414_combined_df<- as.data.frame(res_table_CHD414_combined)
res_table_CHD427_combined_df<- as.data.frame(res_table_CHD427_combined)

#2017_10_04_Try to make scatter plots between the lines
pdf("2018_08_02_hetsvshomo_scatter.pdf")
plot(res_table_hets_combined_df$log2FoldChange, res_table_G1_combined_df$log2FoldChange, main = "CHD7 +/- vs CHD7 -/-", xlab= "CHD7 +/- log2foldchange", ylab = "CHD7 -/- log2foldchange") 
abline(fit<-lm(res_table_hets_combined_df$log2FoldChange~res_table_G1_combined_df$log2FoldChange), col = "red")
legend("bottomright", bty = "n", legend=paste("R2=",format(summary(fit)$adj.r.squared, digits=4)))
dev.off()


#hets vs 1G1
#commongeneslof_table_up<- cbind(res_table_hets_combined_df[rownames(res_table_hets_combined_df) %in% commongeneslof_up,],res_table_G1_combined_df[rownames(res_table_G1_combined_df) %in% commongeneslof_up,])
#commongeneslof_table_down<- cbind(res_table_hets_combined_df[rownames(res_table_hets_combined_df) %in% commongeneslof_down,],res_table_G1_combined_df[rownames(res_table_G1_combined_df) %in% commongeneslof_down,])
#commongeneslof_table<- rbind(commongeneslof_table_up, commongeneslof_table_down)
#write.table(commongeneslof, file= "2018_08_02_CombinedHetsGenesNamesCommonLOF.txt", row.names = F, quote = F)
#write.csv(commongeneslof_table, file= "2018_08_02_GenesNamesCommonLOF_Table.csv", row.names = T, quote = F)


#MISSENSE
#CHD414 to both hets
commongeneshetCHD414_up<- A.cds.sig.up[which(A.cds.sig.up%in% D.cds.sig.up)]
commongeneshetCHD414_table_up<-cbind(res_table_hets_combined_df[rownames(res_table_hets_combined_df) %in% commongeneshetCHD414_up,],res_table_CHD414_combined_df[rownames(res_table_CHD414_combined_df) %in% commongeneshetCHD414_up,])
colnames(commongeneshetCHD414_table_up)<- c("baseMean_hets", "log2FoldChange_hets", "lfcSE_hets", "stat_hets", "pvalue_hets","padj_hets","threshold_hets", "threshold_UP_hets","threshold_DOWN_hets","baseMean_CHD414", "log2FoldChange_CHD414", "lfcSE_CHD414", "stat_CHD414", "pvalue_CHD414","padj_CHD414","threshold_CHD414", "threshold_UP_CHD414","threshold_DOWN_CHD414")  
write.table(commongeneshetCHD414_up, file= "2018_08_03_GenesNamesCommonHets_CHD414Missense_UP.txt", row.names = F, quote = F)
write.csv(commongeneshetCHD414_table_up, file= "2018_08_03_GenesNamesCommonHet_CHD414Missense_Table_UP.csv", row.names = T, quote = F)

commongeneshetCHD414_down<- A.cds.sig.down[which(A.cds.sig.down%in% D.cds.sig.down)]
commongeneshetCHD414_table_down<-cbind(res_table_hets_combined_df[rownames(res_table_hets_combined_df) %in% commongeneshetCHD414_down,],res_table_CHD414_combined_df[rownames(res_table_CHD414_combined_df) %in% commongeneshetCHD414_down,])
colnames(commongeneshetCHD414_table_down)<- c("baseMean_hets", "log2FoldChange_hets", "lfcSE_hets", "stat_hets", "pvalue_hets","padj_hets","threshold_hets", "threshold_down_hets","threshold_DOWN_hets","baseMean_CHD414", "log2FoldChange_CHD414", "lfcSE_CHD414", "stat_CHD414", "pvalue_CHD414","padj_CHD414","threshold_CHD414", "threshold_down_CHD414","threshold_DOWN_CHD414")  
write.table(commongeneshetCHD414_down, file= "2018_08_03_GenesNamesCommonHets_CHD414Missense_down.txt", row.names = F, quote = F)
write.csv(commongeneshetCHD414_table_down, file= "2018_08_03_GenesNamesCommonHet_CHD414Missense_Table_down.csv", row.names = T, quote = F)

commongeneshetCHD414<- c(commongeneshetCHD414_up,commongeneshetCHD414_down)
#commongeneshetCHD414_table<- rbind(commongeneshetCHD414_table_up, commongeneshetCHD414_table_down)
write.table(commongeneshetCHD414, file= "2018_08_08_GenesNamesCommonHets_CHD414Missense.txt", row.names = F, quote = F)
#write.csv(commongeneshetCHD414_table, file= "2018_08_08_GenesNamesCommonHets_CHD414Missense_TABLE.csv", row.names = T, quote = F)

#CHD414 to 1G1
commongenes_CHD414_1G1_up<- C.cds.sig.up[which(C.cds.sig.up %in% D.cds.sig.up)]
commongenes_CHD414_1G1_table_up<- cbind(res_table_G1_combined_df[rownames(res_table_G1_combined_df) %in% commongenes_CHD414_1G1_up,],res_table_CHD414_combined_df[rownames(res_table_CHD414_combined_df) %in% commongenes_CHD414_1G1_up,])
colnames(commongenes_CHD414_1G1_table_up)<- c("baseMean_1G1", "log2FoldChange_1G1", "lfcSE_1G1", "stat_1G1", "pvalue_1G1","padj_1G1","threshold_1G1", "threshold_UP_1G1","threshold_DOWN_1G1","baseMean_CHD414", "log2FoldChange_CHD414", "lfcSE_CHD414", "stat_CHD414", "pvalue_CHD414","padj_CHD414","threshold_CHD414", "threshold_UP_CHD414","threshold_DOWN_CHD414") 
write.table(commongenes_CHD414_1G1_up, file= "2018_08_07_GenesNamesCommon1G1_CHD414_UP.txt", row.names = F, quote = F)
write.csv(commongenes_CHD414_1G1_table_up, file= "2018_08_07_GenesNamesCommon1G1_CHD414_UP_TABLE.csv", row.names = T, quote = F)

commongenes_CHD414_1G1_down<- C.cds.sig.down[which(C.cds.sig.down %in% D.cds.sig.down)]
commongenes_CHD414_1G1_table_down<- cbind(res_table_G1_combined_df[rownames(res_table_G1_combined_df) %in% commongenes_CHD414_1G1_down,],res_table_CHD414_combined_df[rownames(res_table_CHD414_combined_df) %in% commongenes_CHD414_1G1_down,])
colnames(commongenes_CHD414_1G1_table_down)<- c("baseMean_1G1", "log2FoldChange_1G1", "lfcSE_1G1", "stat_1G1", "pvalue_1G1","padj_1G1","threshold_1G1", "threshold_UP_1G1","threshold_DOWN_1G1","baseMean_CHD414", "log2FoldChange_CHD414", "lfcSE_CHD414", "stat_CHD414", "pvalue_CHD414","padj_CHD414","threshold_CHD414", "threshold_UP_CHD414","threshold_DOWN_CHD414") 
write.table(commongenes_CHD414_1G1_down, file= "2018_08_07_GenesNamesCommon1G1_CHD414_DOWN.txt", row.names = F, quote = F)
write.csv(commongenes_CHD414_1G1_table_down, file= "2018_08_07_GenesNamesCommon1G1_CHD414_DOWN_TABLE.csv", row.names = T, quote = F)

commongenes_CHD414_1G1<- c(commongenes_CHD414_1G1_up,commongenes_CHD414_1G1_down)
commongenes_CHD414_1G1_table<- rbind(commongenes_CHD414_1G1_table_up, commongenes_CHD414_1G1_table_down)
write.table(commongenes_CHD414_1G1, file= "2018_08_07_GenesNamesCommon1G1_CHD414.txt", row.names = F, quote = F)
write.csv(commongenes_CHD414_1G1_table, file= "2018_08_07_GenesNamesCommon1G1_CHD414_TABLE.csv", row.names = T, quote = F)

#STOP HERE
#CHD427 to both hets
commongeneshetCHD427missense_up<- A.cds.sig.up[which(A.cds.sig.up%in% E.cds.sig.up)]
commongeneshetCHD427missense_table_up<-cbind(res_table_hets_combined_df[rownames(res_table_hets_combined_df) %in% commongeneshetCHD427missense_up,],res_table_CHD427_combined_df[rownames(res_table_CHD427_combined_df) %in% commongeneshetCHD427missense_up,])
colnames(commongeneshetCHD427missense_table_up)<- c("baseMean_hets", "log2FoldChange_hets", "lfcSE_hets", "stat_hets", "pvalue_hets","padj_hets","threshold_hets", "threshold_UP_hets","threshold_DOWN_hets","baseMean_CHD427", "log2FoldChange_CHD427", "lfcSE_CHD427", "stat_CHD427", "pvalue_CHD427","padj_CHD427","threshold_CHD427", "threshold_UP_CHD427","threshold_DOWN_CHD427")  
write.table(commongeneshetCHD427missense_up, file= "2018_08_03_GenesNamesCommonHets_CHD427Missense_UP.txt", row.names = F, quote = F)
write.csv(commongeneshetCHD427missense_table_up, file= "2018_08_03_GenesNamesCommonHet_CHD427Missense_Table_UP.csv", row.names = T, quote = F)

commongeneshetCHD427missense_down<- A.cds.sig.down[which(A.cds.sig.down%in% E.cds.sig.down)]
commongeneshetCHD427missense_table_down<-cbind(res_table_hets_combined_df[rownames(res_table_hets_combined_df) %in% commongeneshetCHD427missense_down,],res_table_CHD427_combined_df[rownames(res_table_CHD427_combined_df) %in% commongeneshetCHD427missense_down,])
colnames(commongeneshetCHD427missense_table_down)<- c("baseMean_hets", "log2FoldChange_hets", "lfcSE_hets", "stat_hets", "pvalue_hets","padj_hets","threshold_hets", "threshold_UP_hets","threshold_DOWN_hets","baseMean_CHD427", "log2FoldChange_CHD427", "lfcSE_CHD427", "stat_CHD427", "pvalue_CHD427","padj_CHD427","threshold_CHD427", "threshold_UP_CHD427","threshold_DOWN_CHD427")
write.table(commongeneshetCHD427missense_down, file= "2018_08_03_GenesNamesCommonHets_CHD427Missense_DOWN.txt", row.names = F, quote = F)
write.csv(commongeneshetCHD427missense_table_down, file= "2018_08_03_GenesNamesCommonHet_CHD427Missense_Table_DOWN.csv", row.names = T, quote = F)

commongenes_CHD427_hets<- c(commongeneshetCHD427missense_up,commongeneshetCHD427missense_down)
commongenes_CHD427_hets_table<- rbind(commongeneshetCHD427missense_table_up, commongeneshetCHD427missense_table_down)
write.table(commongenes_CHD427_hets, file= "2018_01_25_GenesNamesCommonhets_CHD427.txt", row.names = F, quote = F)
write.csv(commongenes_CHD427_hets_table, file= "2018_01_25_GenesNamesCommonhets_CHD427_TABLE.csv", row.names = T, quote = F)


#CHD427 to CHD414
commongenesmissense_up<- D.cds.sig.up[which(D.cds.sig.up %in% E.cds.sig.up)]
commongenesmissense_table_up<- cbind(res_table_CHD414_combined_df[rownames(res_table_CHD414_combined_df) %in% commongenesmissense_up,],res_table_CHD427_combined_df[rownames(res_table_CHD427_combined_df) %in% commongenesmissense_up,])
colnames(commongenesmissense_table_up)<- c("baseMean_CHD414", "log2FoldChange_CHD414", "lfcSE_CHD414", "stat_CHD414", "pvalue_CHD414","padj_CHD414","threshold_CHD414", "threshold_UP_CHD414","threshold_DOWN_CHD414","baseMean_CHD427/CHD427", "log2FoldChange_CHD427/CHD427", "lfcSE_CHD427/CHD427", "stat_CHD427/CHD427", "pvalue_CHD427/CHD427","padj_CHD427/CHD427","threshold_CHD427/CHD427", "threshold_UP_CHD427/CHD427","threshold_DOWN_CHD427/CHD427")  
write.table(commongenesmissense_up, file= "2018_03_07_GenesCommonMissense_UP.txt", row.names = F, quote = F)
write.csv(commongenesmissense_table_up, file= "2018_03_07_GenesNamesCommon_Missense_Table_UP.csv", row.names = T, quote = F)

commongenesmissense_down<- D.cds.sig.down[which(D.cds.sig.down %in% E.cds.sig.down)]
commongenesmissense_table_down<- cbind(res_table_CHD414_combined_df[rownames(res_table_CHD414_combined_df) %in% commongenesmissense_down,],res_table_CHD427_combined_df[rownames(res_table_CHD427_combined_df) %in% commongenesmissense_down,])
colnames(commongenesmissense_table_down)<- c("baseMean_CHD414", "log2FoldChange_CHD414", "lfcSE_CHD414", "stat_CHD414", "pvalue_CHD414","padj_CHD414","threshold_CHD414", "threshold_UP_CHD414","threshold_DOWN_CHD414","baseMean_CHD427/CHD427", "log2FoldChange_CHD427/CHD427", "lfcSE_CHD427/CHD427", "stat_CHD427/CHD427", "pvalue_CHD427/CHD427","padj_CHD427/CHD427","threshold_CHD427/CHD427", "threshold_UP_CHD427/CHD427","threshold_DOWN_CHD427/CHD427")  
write.table(commongenesmissense_down, file= "2018_03_07_GenesCommonMissense_DOWN.txt", row.names = F, quote = F)
write.csv(commongenesmissense_table_down, file= "2018_03_07_GenesNamesCommon_Missense_Table_DOWN.csv", row.names = T, quote = F)

commongenesmissense<- c(commongenesmissense_up,commongenesmissense_down)
commongenes_missense_table<- rbind(commongenesmissense_table_up, commongenesmissense_table_down)
write.table(commongenesmissense, file= "2018_03_07_GenesNamesCommonMissense.txt", row.names = F, quote = F)
write.csv(commongenesmissense, file= "2018_03_07_GenesNamesCommonMissense_TABLE.csv", row.names = T, quote = F)

###STOP HERE
#1G1 vs CHD427/CHD427
commongenes1G1_CHD427_up<- C.cds.sig.up[which(C.cds.sig.up %in% E.cds.sig.up)]
commongenes1G1_CHD427_table_up<- cbind(res_table_G1_combined_df[rownames(res_table_G1_combined_df) %in% commongeneshomomissense_up,],res_table_CHD427_combined_df[rownames(res_table_CHD427_combined_df) %in% commongeneshomomissense_up,])
colnames(commongenes1G1_CHD427_table_up)<- c("baseMean_1G1", "log2FoldChange_1G1", "lfcSE_1G1", "stat_1G1", "pvalue_1G1","padj_1G1","threshold_1G1", "threshold_UP_1G1","threshold_DOWN_1G1","baseMean_CHD427/CHD427", "log2FoldChange_CHD427/CHD427", "lfcSE_CHD427/CHD427", "stat_CHD427/CHD427", "pvalue_CHD427/CHD427","padj_CHD427/CHD427","threshold_CHD427/CHD427", "threshold_UP_CHD427/CHD427","threshold_DOWN_CHD427/CHD427")  
write.table(commongenes1G1_CHD427_up, file= "2018_08_03_GenesCommon1G1_CHD427_UP.txt", row.names = F, quote = F)
write.csv(commongenes1G1_CHD427_table_up, file= "2018_08_03_GenesNames1G1_CHD427_Table_UP.csv", row.names = T, quote = F)

commongenes1G1_CHD427_down<- C.cds.sig.down[which(C.cds.sig.down %in% E.cds.sig.down)]
commongenes1G1_CHD427_table_down<- cbind(res_table_G1_combined_df[rownames(res_table_G1_combined_df) %in% commongeneshomomissense_down,],res_table_CHD427_combined_df[rownames(res_table_CHD427_combined_df) %in% commongeneshomomissense_down,])
colnames(commongenes1G1_CHD427_table_down)<- c("baseMean_1G1", "log2FoldChange_1G1", "lfcSE_1G1", "stat_1G1", "pvalue_1G1","padj_1G1","threshold_1G1", "threshold_UP_1G1","threshold_DOWN_1G1","baseMean_CHD427/CHD427", "log2FoldChange_CHD427/CHD427", "lfcSE_CHD427/CHD427", "stat_CHD427/CHD427", "pvalue_CHD427/CHD427","padj_CHD427/CHD427","threshold_CHD427/CHD427", "threshold_UP_CHD427/CHD427","threshold_DOWN_CHD427/CHD427")  
write.table(commongenes1G1_CHD427_down, file= "2018_08_03_GenesCommon1G1_CHD427_DOWN.txt", row.names = F, quote = F)
write.csv(commongenes1G1_CHD427_table_down, file= "2018_08_03_GenesNames1G1_CHD427_Table_DOWN.csv", row.names = T, quote = F)

commongenes1G1_CHD427<- c(commongenes1G1_CHD427_up,commongenes1G1_CHD427_down)
commongenes1G1_CHD427_table<- rbind(commongenes1G1_CHD427_table_up, commongenes1G1_CHD427_table_down)
write.table(commongenes1G1_CHD427, file= "2018_08_03_GenesNamesCommon1G1_CHD427.txt", row.names = F, quote = F)
write.csv(commongenes1G1_CHD427_table, file= "2018_08_-3_GenesNamesCommon1G1_CHD427_TABLE.csv", row.names = T, quote = F)



#het and D2355 and CHD427
commongeneshetCHD427CHD414_up<- A.cds.sig.up[which(A.cds.sig.up%in% commongenesmissense_up)]
commongeneshetCHD427CHD414_table_up<-cbind(commongeneshet_table[rownames(commongeneshet_table) %in% commongeneshetCHD427CHD414_up,],res_table_CHD427_combined_df[rownames(res_table_CHD427_combined_df) %in% commongeneshetCHD427CHD414_up,],res_table_CHD414_combined_df[rownames(res_table_CHD414_combined_df) %in% commongeneshetCHD427CHD414_up,])
colnames(commongeneshetCHD427CHD414_table_up)<- c("baseMean_1D4", "log2FoldChange_1D4", "lfcSE_1D4", "stat_1D4", "pvalue_1D4","padj_1D4","threshold_1D4", "threshold_UP_1D4","threshold_DOWN_1D4","baseMean_1G6", "log2FoldChange_1G6", "lfcSE_1G6", "stat_1G6", "pvalue_1G6","padj_1G6","threshold_1G6", "threshold_UP_1G6","threshold_DOWN_1G6","baseMean_CHD427", "log2FoldChange_CHD427", "lfcSE_CHD427", "stat_CHD427", "pvalue_CHD427","padj_CHD427","threshold_CHD427", "threshold_UP_CHD427","threshold_DOWN_CHD427", "baseMean_CHD414", "log2FoldChange_CHD414", "lfcSE_CHD414", "stat_CHD414", "pvalue_CHD414","padj_CHD414","threshold_CHD414", "threshold_UP_CHD414","threshold_DOWN_CHD414")  
write.table(commongeneshetCHD427CHD414_table_up, file= "2018_03_14_GenesNamesCommoncommongeneshetCHD427CHD414_UP.txt", row.names = F, quote = F)
write.csv(commongeneshetCHD427CHD414_table_up, file= "2018_03_14_GenesNamesCommoncommongeneshetCHD427CHD414_UP.csv", row.names = T, quote = F)

commongeneshetCHD427CHD414_down<- A.cds.sig.down[which(A.cds.sig.down%in% commongenesmissense_down)]
commongeneshetCHD427CHD414_table_down<-cbind(commongeneshet_table[rownames(commongeneshet_table) %in% commongeneshetCHD427CHD414_down,],res_table_CHD427_combined_df[rownames(res_table_CHD427_combined_df) %in% commongeneshetCHD427CHD414_down,],res_table_CHD414_combined_df[rownames(res_table_CHD414_combined_df) %in% commongeneshetCHD427CHD414_down,])
colnames(commongeneshetCHD427CHD414_table_down)<- c("baseMean_1D4", "log2FoldChange_1D4", "lfcSE_1D4", "stat_1D4", "pvalue_1D4","padj_1D4","threshold_1D4", "threshold_UP_1D4","threshold_DOWN_1D4","baseMean_1G6", "log2FoldChange_1G6", "lfcSE_1G6", "stat_1G6", "pvalue_1G6","padj_1G6","threshold_1G6", "threshold_UP_1G6","threshold_DOWN_1G6","baseMean_CHD427", "log2FoldChange_CHD427", "lfcSE_CHD427", "stat_CHD427", "pvalue_CHD427","padj_CHD427","threshold_CHD427", "threshold_UP_CHD427","threshold_DOWN_CHD427", "baseMean_CHD414", "log2FoldChange_CHD414", "lfcSE_CHD414", "stat_CHD414", "pvalue_CHD414","padj_CHD414","threshold_CHD414", "threshold_UP_CHD414","threshold_DOWN_CHD414")  
write.table(commongeneshetCHD427CHD414_table_down, file= "2018_03_14_GenesNamesCommoncommongeneshetCHD427CHD414_DOWN.txt", row.names = F, quote = F)
write.csv(commongeneshetCHD427CHD414_table_down, file= "2018_03_14_GenesNamesCommoncommongeneshetCHD427CHD414_DOWN.csv", row.names = T, quote = F)

commongeneshetCHD427CHD414<- c(commongeneshetCHD427CHD414_up,commongeneshetCHD427CHD414_down)
commongeneshetCHD427CHD414_table<- rbind(commongeneshetCHD427CHD414_table_up, commongeneshetCHD427CHD414_down)
write.table(commongeneshetCHD427CHD414_table, file= "2018_03_14_commongeneshetCHD427CHD414.txt", row.names = F, quote = F)
write.csv(commongeneshetCHD427CHD414_table, file= "2018_03_14_commongeneshetCHD427CHD414_TABLE.csv", row.names = T, quote = F)

#Is this overlap statistically significant
phyper(471,1804,(57373-1804),834, lower.tail = F, log.p = TRUE)

##STOP HERE
CDE<-c(C.cds.sig, D.cds.sig, E.cds.sig)
het_uniquegenes<- setdiff(commongeneshet,CDE)
hetuniquegenes_table<- commongeneshet_table[rownames(commongeneshet_table) %in% het_uniquegenes,]
write.table(het_uniquegenes, file= "2018_03_07_GenesUnique_HET.txt", row.names = F, quote = F)
write.csv(hetuniquegenes_table, file= "2018_03_07_GenesNamesUnique_HET_Table.csv", row.names = T, quote = F)

DE<- c(D.cds.sig, E.cds.sig)
lof_uniquegenes<- setdiff(commongeneslof,DE)
lof_uniquegenes_table<- commongeneslof_table[rownames(commongeneslof_table) %in% lof_uniquegenes, ]
write.table(lof_uniquegenes, file= "2018_03_07_GenesUnique_LOF.txt", row.names = F, quote = F)
write.csv(lof_uniquegenes, file= "2018_03_07_GenesNamesUnique_LOF_Table.csv", row.names = T, quote = F)

ABDE<- c(A.cds.sig, B.cds.sig, D.cds.sig, E.cds.sig)
homo_uniquegenes<- setdiff(C.cds.sig, ABDE)
homo_uniquegenes_table<- res_table_G1_combined_df[rownames(res_table_G1_combined_df) %in% homo_uniquegenes,]
write.table(homo_uniquegenes, file= "2018_03_07_GenesUnique_HOMO.txt", row.names = F, quote = F)
write.csv(homo_uniquegenes_table, file= "2018_03_07_GenesNamesUnique_HOMO_Table.csv", row.names = T, quote = F)

ABC<- c(A.cds.sig,B.cds.sig,C.cds.sig)
missense_uniquegenes<- setdiff(commongenesmissense,ABC)
missense_uniquegenes_table<- commongenes_CHD427_hets_table[rownames(commongenes_CHD427_hets_table) %in% missense_uniquegenes,]
write.table(missense_uniquegenes, file= "2018_03_07_GenesUnique_MISSENSE.txt", row.names = F, quote = F)
write.csv(missense_uniquegenes_table, file= "2018_03_07_GenesNamesUnique_MISSENSE_Table.csv", row.names = T, quote = F)

# Subsetting dataset to only include significant genes with padj < 0.05

#Run gprofiler
### Functional analysis using gProfileR (some of these are defaults; check help pages) 
install.packages("gProfileR")
library(gProfileR)

# Running gprofiler to identify enriched processes among significant genes
gprofiler_results_hets <- gprofiler(query = A.cds.sig, 
                                    organism = "hsapiens",
                                    ordered_query = F, 
                                    exclude_iea = F, 
                                    max_p_value = 0.05, 
                                    max_set_size = 0,
                                    correction_method = "fdr",
                                    hier_filtering = "none", 
                                    domain_size = "annotated",
                                    custom_bg = "")

write.table(gprofiler_results_hets,"2018_08_02_CHD7_D3U_bcbio_gprofiler_hets_combined.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_hets <- cbind(gprofiler_results_hets$term.id, gprofiler_results_hets$p.value)
GOs_hets <- allterms_hets[grep('GO:', allterms_hets),]
write.table(GOs_hets, "2018_08_02_CHD7_D3U_bcbio_GOs_hets_combined.txt", sep="\t", quote=F, row.names=F, col.names=F)


gprofiler_results_G1 <- gprofiler(query = C.cds.sig, 
                                  organism = "hsapiens",
                                  ordered_query = F, 
                                  exclude_iea = F, 
                                  max_p_value = 0.05, 
                                  max_set_size = 0,
                                  correction_method = "fdr",
                                  hier_filtering = "none", 
                                  domain_size = "annotated",
                                  custom_bg = "")

write.table(gprofiler_results_G1,"2018_08_02_CHD7_D3U_bcbio_gprofiler_G1_combined.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_G1 <- cbind(gprofiler_results_G1$term.id, gprofiler_results_G1$p.value)
GOs_G1 <- allterms_G1[grep('GO:', allterms_G1),]
write.table(GOs_G1, "2018_08_02_CHD7_D3U_bcbio_GOs_G1_combined.txt", sep="\t", quote=F, row.names=F, col.names=F)

gprofiler_results_D2355 <- gprofiler(query = D.cds.sig, 
                                     organism = "hsapiens",
                                     ordered_query = F, 
                                     exclude_iea = F, 
                                     max_p_value = 0.05, 
                                     max_set_size = 0,
                                     correction_method = "fdr",
                                     hier_filtering = "none", 
                                     domain_size = "annotated",
                                     custom_bg = "")

write.table(gprofiler_results_D2355,"2018_03_08_CHD7_D3U_bcbio_gprofiler_D2355_combined.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_2355 <- cbind(gprofiler_results_D2355$term.id, gprofiler_results_D2355$p.value)
GOs_2355 <- allterms_2355[grep('GO:', allterms_2355),]
write.table(GOs_2355, "2018_03_08_CHD7_D3U_bcbio_GOs_2355_combined.txt", sep="\t", quote=F, row.names=F, col.names=F)

gprofiler_results_R2111 <- gprofiler(query = E.cds.sig, 
                                     organism = "hsapiens",
                                     ordered_query = F, 
                                     exclude_iea = F, 
                                     max_p_value = 0.05, 
                                     max_set_size = 0,
                                     correction_method = "fdr",
                                     hier_filtering = "none", 
                                     domain_size = "annotated",
                                     custom_bg = "")

write.table(gprofiler_results_R2111,"2018_01_25_CHD7_D3U_bcbio_gprofiler_CHD427_combined.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_2111 <- cbind(gprofiler_results_R2111$term.id, gprofiler_results_R2111$p.value)
GOs_2111 <- allterms_2111[grep('GO:', allterms_2111),]
write.table(GOs_2111, "2018_01_25_CHD7_D3U_bcbio_GOs_R2111_combined.txt", sep="\t", quote=F, row.names=F, col.names=F)


gprofiler_results_hets_CHD414 <- gprofiler(query = commongeneshetCHD414, 
                                           organism = "hsapiens",
                                           ordered_query = F, 
                                           exclude_iea = F, 
                                           max_p_value = 0.05, 
                                           max_set_size = 0,
                                           correction_method = "fdr",
                                           hier_filtering = "none", 
                                           domain_size = "annotated",
                                           custom_bg = "")

write.table(gprofiler_results_hets_CHD414,"2017_08_03_CHD7_D3U_bcbio_gprofiler_hets_CHD414.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_hets_CHD414 <- cbind(gprofiler_results_hets_CHD414$term.id, gprofiler_results_hets_CHD414$p.value)
GOs_hets_CHD414 <- allterms_hets_CHD414[grep('GO:', allterms_hets_CHD414),]
write.table(GOs_hets_CHD414, "2018_08_03_CHD7_D3U_bcbio_GOs_hets_CHD414.txt", sep="\t", quote=F, row.names=F, col.names=F)

gprofiler_results_hets_CHD427 <- gprofiler(query = commongenes_CHD427_hets, 
                                           organism = "hsapiens",
                                           ordered_query = F, 
                                           exclude_iea = F, 
                                           max_p_value = 0.05, 
                                           max_set_size = 0,
                                           correction_method = "fdr",
                                           hier_filtering = "none", 
                                           domain_size = "annotated",
                                           custom_bg = "")

write.table(gprofiler_results_hets_CHD427,"2017_08_03_CHD7_D3U_bcbio_gprofiler_hets_CHD427.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_hets_CHD427 <- cbind(gprofiler_results_hets_CHD427$term.id, gprofiler_results_hets_CHD427$p.value)
GOs_hets_CHD427 <- allterms_hets_CHD427[grep('GO:', allterms_hets_CHD427),]
write.table(GOs_hets_CHD427, "2018_08_03_CHD7_D3U_bcbio_GOs_hets_CHD427.txt", sep="\t", quote=F, row.names=F, col.names=F)


#
gprofiler_results_hets_missense <- gprofiler(query = commongeneshetCHD427CHD414, 
                                             organism = "hsapiens",
                                             ordered_query = F, 
                                             exclude_iea = F, 
                                             max_p_value = 0.05, 
                                             max_set_size = 0,
                                             correction_method = "fdr",
                                             hier_filtering = "none", 
                                             domain_size = "annotated",
                                             custom_bg = "")

write.table(gprofiler_results_hets_missense,"2017_03_14_CHD7_D3U_bcbio_gprofiler_hets_missense.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_hets_missense <- cbind(gprofiler_results_hets_missense$term.id, gprofiler_results_hets_missense$p.value)
GOs_hets_missense <- allterms_hets_missense[grep('GO:', allterms_hets_missense),]
write.table(GOs_hets_missense, "2018_03_14_CHD7_D3U_bcbio_GOs_hets_missense_combined.txt", sep="\t", quote=F, row.names=F, col.names=F)

gprofiler_results_lof <- gprofiler(query = commongeneslof, 
                                   organism = "hsapiens",
                                   ordered_query = F, 
                                   exclude_iea = F, 
                                   max_p_value = 0.05, 
                                   max_set_size = 0,
                                   correction_method = "fdr",
                                   hier_filtering = "none", 
                                   domain_size = "annotated",
                                   custom_bg = "")

write.table(gprofiler_results_lof,"2018_08_02_CHD7_D3U_bcbio_gprofiler_lof.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_lof <- cbind(gprofiler_results_lof$term.id, gprofiler_results_lof$p.value)
GOs_lof <- allterms_lof[grep('GO:', allterms_lof),]
write.table(GOs_lof, "2018_08_02_CHD7_D3U_bcbio_GOs_lof.txt", sep="\t", quote=F, row.names=F, col.names=F)


gprofiler_results_missense <- gprofiler(query = commongenesmissense, 
                                        organism = "hsapiens",
                                        ordered_query = F, 
                                        exclude_iea = F, 
                                        max_p_value = 0.05, 
                                        max_set_size = 0,
                                        correction_method = "fdr",
                                        hier_filtering = "none", 
                                        domain_size = "annotated",
                                        custom_bg = "")

write.table(gprofiler_results_missense,"2018_03_14_CHD7_D3U_bcbio_gprofiler_missense.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_missense <- cbind(gprofiler_results_missense$term.id, gprofiler_results_missense$p.value)
GOs_missense <- allterms_missense[grep('GO:', allterms_missense),]
write.table(GOs_missense, "2018_03_14_CHD7_D3U_bcbio_GOs_missense.txt", sep="\t", quote=F, row.names=F, col.names=F)

#gprofiler_results_all <- gprofiler(query = rownames(sig_genes_all), 
#                                        organism = "hsapiens",
#                                        ordered_query = F, 
#                                        exclude_iea = F, 
#                                       max_p_value = 0.05, 
#                                      max_set_size = 0,
#                                      correction_method = "fdr",
#                                     hier_filtering = "none", 
#                                        domain_size = "annotated",
#                                        custom_bg = "")

#write.table(gprofiler_results_all,"2017_09_15_CHD7_D3U_bcbio_gprofiler_all.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
#allterms_all <- cbind(gprofiler_results_all$term.id, gprofiler_results_all$p.value)
#GOs_all <- allterms_all[grep('GO:', allterms_all),]
#write.table(GOs_all, "2017_09_15_CHD7_D3U_bcbio_GOs_all.txt", sep="\t", quote=F, row.names=F, col.names=F)

#gprofiler_results_all <- gprofiler(query = rownames(sig_genes_all), 
#                                   organism = "hsapiens",
#                                   ordered_query = F, 
#                                   exclude_iea = F, 
#                                   max_p_value = 0.05, 
#                                   max_set_size = 0,
#                                   correction_method = "fdr",
#                                   hier_filtering = "none", 
#                                   domain_size = "annotated",
#                                   custom_bg = "")

#write.table(gprofiler_results_all,"2017_09_15_CHD7_D3U_bcbio_gprofiler_all.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
#allterms_all <- cbind(gprofiler_results_all$term.id, gprofiler_results_all$p.value)
#GOs_all <- allterms_all[grep('GO:', allterms_all),]
#write.table(GOs_all, "2017_09_15_CHD7_D3U_bcbio_GOs_all.txt", sep="\t", quote=F, row.names=F, col.names=F)

#Gprofiler for unique genes
gprofiler_results_het_unique <- gprofiler(query = het_uniquegenes, 
                                          organism = "hsapiens",
                                          ordered_query = F, 
                                          exclude_iea = F, 
                                          max_p_value = 0.05, 
                                          max_set_size = 0,
                                          correction_method = "fdr",
                                          hier_filtering = "none", 
                                          domain_size = "annotated",
                                          custom_bg = "")
write.table(gprofiler_results_het_unique,"2017_09_15_CHD7_D3U_bcbio_gprofiler_HET_UNIQUE.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_het_unique <- cbind(gprofiler_results_het_unique$term.id, gprofiler_results_het_unique$p.value)
GOs_het_unique <- allterms_het_unique[grep('GO:', allterms_het_unique),]
write.table(GOs_het_unique, "2017_09_15_CHD7_D3U_bcbio_GOs_HET_UNIQUE.txt", sep="\t", quote=F, row.names=F, col.names=F)

#gprofiler_results_hetmissense_unique <- gprofiler(query = rownames(sig_genes_het_missense_unique), 
#                                          organism = "hsapiens",
#                                         ordered_query = F, 
#                                        exclude_iea = F, 
#                                       max_p_value = 0.05, 
#                                      max_set_size = 0,
#                                     correction_method = "fdr",
#                                    hier_filtering = "none", 
#                                   domain_size = "annotated",
#                                  custom_bg = "")

#write.table(gprofiler_results_hetmissense_unique,"2017_09_15_CHD7_D3U_bcbio_gprofiler_HETMISSENSE_UNIQUE.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
#allterms_het_missense_unique <- cbind(gprofiler_results_hetmissense_unique$term.id, gprofiler_results_hetmissense_unique$p.value)
#GOs_hetmissense_unique <- allterms_het_missense_unique[grep('GO:', allterms_het_missense_unique),]
#write.table(GOs_hetmissense_unique, "2017_09_15_CHD7_D3U_bcbio_GOs_HETMISSENSE_UNIQUE.txt", sep="\t", quote=F, row.names=F, col.names=F)

gprofiler_results_homomissense_unique <- gprofiler(query = rownames(sig_genes_homo_missense_unique), 
                                                   organism = "hsapiens",
                                                   ordered_query = F, 
                                                   exclude_iea = F, 
                                                   max_p_value = 0.05, 
                                                   max_set_size = 0,
                                                   correction_method = "fdr",
                                                   hier_filtering = "none", 
                                                   domain_size = "annotated",
                                                   custom_bg = "")

write.table(gprofiler_results_homomissense_unique,"2017_09_15_CHD7_D3U_bcbio_gprofiler_HOMOMISSENSE_UNIQUE.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_homo_missense_unique <- cbind(gprofiler_results_homomissense_unique$term.id, gprofiler_results_homomissense_unique$p.value)
GOs_homomissense_unique <- allterms_homo_missense_unique[grep('GO:', allterms_homo_missense_unique),]
write.table(GOs_homomissense_unique, "2017_09_15_CHD7_D3U_bcbio_GOs_HOMOMISSENSE_UNIQUE.txt", sep="\t", quote=F, row.names=F, col.names=F)

gprofiler_results_missense_unique <- gprofiler(query = rownames(sig_genes_missense_unique), 
                                               organism = "hsapiens",
                                               ordered_query = F, 
                                               exclude_iea = F, 
                                               max_p_value = 0.05, 
                                               max_set_size = 0,
                                               correction_method = "fdr",
                                               hier_filtering = "none", 
                                               domain_size = "annotated",
                                               custom_bg = "")

write.table(gprofiler_results_missense_unique,"2017_09_15_CHD7_D3U_bcbio_gprofiler_MISSENSE_UNIQUE.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_missense_unique <- cbind(gprofiler_results_missense_unique$term.id, gprofiler_results_missense_unique$p.value)
GOs_missense_unique <- allterms_missense_unique[grep('GO:', allterms_missense_unique),]
write.table(GOs_missense_unique, "2017_09_15_CHD7_D3U_bcbio_GOs_MISSENSE_UNIQUE.txt", sep="\t", quote=F, row.names=F, col.names=F)

gprofiler_results_homo_unique <- gprofiler(query = rownames(sig_genes_homo_unique), 
                                           organism = "hsapiens",
                                           ordered_query = F, 
                                           exclude_iea = F, 
                                           max_p_value = 0.05, 
                                           max_set_size = 0,
                                           correction_method = "fdr",
                                           hier_filtering = "none", 
                                           domain_size = "annotated",
                                           custom_bg = "")

write.table(gprofiler_results_homo_unique,"2017_09_15_CHD7_D3U_bcbio_gprofiler_HOMO_UNIQUE.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_homo_unique <- cbind(gprofiler_results_homo_unique$term.id, gprofiler_results_homo_unique$p.value)
GOs_homo_unique <- allterms_homo_unique[grep('GO:', allterms_homo_unique),]
write.table(GOs_homo_unique, "2017_09_15_CHD7_D3U_bcbio_GOs_homo_UNIQUE.txt", sep="\t", quote=F, row.names=F, col.names=F)

gprofiler_results_LOF_unique <- gprofiler(query = lof_uniquegenes, 
                                          organism = "hsapiens",
                                          ordered_query = F, 
                                          exclude_iea = F, 
                                          max_p_value = 0.05, 
                                          max_set_size = 0,
                                          correction_method = "fdr",
                                          hier_filtering = "none", 
                                          domain_size = "annotated",
                                          custom_bg = "")

write.table(gprofiler_results_LOF_unique,"2017_09_15_CHD7_D3U_bcbio_gprofiler_LOF_UNIQUE.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_lof_unique <- cbind(gprofiler_results_LOF_unique$term.id, gprofiler_results_LOF_unique$p.value)
GOs_lof_unique <- allterms_lof_unique[grep('GO:', allterms_lof_unique),]
write.table(GOs_lof_unique, "2017_09_15_CHD7_D3U_bcbio_GOs_LOF_UNIQUE.txt", sep="\t", quote=F, row.names=F, col.names=F)

gprofiler_results_PCA1 <- gprofiler(query = pca.1, 
                                    organism = "hsapiens",
                                    ordered_query = F, 
                                    exclude_iea = F, 
                                    max_p_value = 0.05, 
                                    max_set_size = 0,
                                    correction_method = "fdr",
                                    hier_filtering = "none", 
                                    domain_size = "annotated",
                                    custom_bg = "")

write.table(gprofiler_results_PCA1,"2017_03_14_CHD7_D3U_bcbio_gprofiler_PCA1.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_pca1 <- cbind(gprofiler_results_PCA1$term.id, gprofiler_results_PCA1$p.value)
GOs_pca1 <- allterms_pca1[grep('GO:', allterms_pca1),]
write.table(GOs_pca1, "2018_03_14_CHD7_D3U_bcbio_GOs_PCA1.txt", sep="\t", quote=F, row.names=F, col.names=F)

gprofiler_results_PCA2 <- gprofiler(query = pca.2, 
                                    organism = "hsapiens",
                                    ordered_query = F, 
                                    exclude_iea = F, 
                                    max_p_value = 0.05, 
                                    max_set_size = 0,
                                    correction_method = "fdr",
                                    hier_filtering = "none", 
                                    domain_size = "annotated",
                                    custom_bg = "")

write.table(gprofiler_results_PCA2,"2017_03_14_CHD7_D3U_bcbio_gprofiler_PCA2.txt", sep="\t", quote=F, row.names=F)
#Extract GO IDs and p values
allterms_pca2 <- cbind(gprofiler_results_PCA2$term.id, gprofiler_results_PCA2$p.value)
GOs_pca2 <- allterms_pca2[grep('GO:', allterms_pca2),]
write.table(GOs_pca2, "2018_03_14_CHD7_D3U_bcbio_GOs_PCA2.txt", sep="\t", quote=F, row.names=F, col.names=F)



pdf("2017_09_14_CHD7_D3U_combined_bcbio_1G1RNAseqvsChIPoverlap.pdf")
venn.plot <- draw.pairwise.venn(area1 = (939+69), area2 = (1264+69), cross.area=69, fill=c("yellow", "purple"), alpha=c(0.5,0.5), cex = 3, cat.fontface=4, category=c("CHD7 -/- RNA-seq", "CHD7 ChIP-genes"), cat.pos = c(200, 160))
grid.draw(venn.plot)
dev.off()

#Grep specific rows to make a counts table for neural and cardiac
bargraph_CTCF<- c("CTCF", "CTCFL")
bargraph_table_CTCF<-as.data.frame(combined_normalized_counts_genotype[rownames(combined_normalized_counts_genotype) %in% bargraph_CTCF,])
bargraph_table_CTCF<- as.data.frame(t(bargraph_table_CTCF))
write.csv(bargraph_table_CTCF, "bargraph_table_CTCF.csv")
bargraph_table_CTCF<- read.csv("bargraph_table_CTCF_2.csv", header = TRUE,)
bargraph_table_CTCF$Genotype<- factor(bargraph_table_CTCF$Genotype, levels = c("PGP1", "CHD7+/-", "CHD7-/-"))
bargraph_table_CTCF$Gene<- factor(bargraph_table_CTCF$Gene)
bargraph_table_CTCFL<- read.csv("bargraph_table_CTCFL.csv", header = TRUE,)
bargraph_table_CTCFL$Genotype<- factor(bargraph_table_CTCFL$Genotype, levels = c("PGP1", "CHD7+/-", "CHD7-/-"))
bargraph_table_CTCFL$Gene<- factor(bargraph_table_CTCFL$Gene)

#Make a boxplot make graph comparing PGP1 vs LOF mutants
theme_set(theme_gray(base_size = 18))
#CTCF
p10<- ggplot(bargraph_table_CTCF, aes(x= Gene, y = Counts, fill = Genotype )) + 
  geom_boxplot()
p10<- p10 + scale_fill_manual(values=c("#636363","#5b0f82", "#feb24c")) + ggtitle("CTCF normalized mRNA counts") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                         panel.background = element_blank(), axis.line = element_line(colour = "black" )) 
p10<- p10 +theme(plot.title = element_text(hjust = 0.5))
ggsave("2018_11_01_mRNA_levels_of_CTCF.tiff", plot = p10 ,  device = "tiff")

#CTCFL
p10<- ggplot(bargraph_table_CTCFL, aes(x= Gene, y = Counts, fill = Genotype )) + 
  geom_boxplot()
p10<- p10 + scale_fill_manual(values=c("#636363","#5b0f82", "#feb24c")) + ggtitle("CTCFL normalized mRNA counts") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                          panel.background = element_blank(), axis.line = element_line(colour = "black" )) 
p10<- p10 +theme(plot.title = element_text(hjust = 0.5))
ggsave("2018_11_01_mRNA_levels_of_CTCFL.tiff", plot = p10 ,  device = "tiff")




#Grep specific rows to make a counts table for neural and cardiac
bargraph_cardiac<- c("SCN1B", "KCNJ2", "KCNE2", "MYBPC2", "MYOM1", "TTN", "FGF13", "SOX9")
bargraph_table_cardiac<-as.data.frame(combined_normalized_counts_cellline[rownames(combined_normalized_counts_cellline) %in% bargraph_cardiac,])
bargraph_table_cardiac<- as.data.frame(t(bargraph_table_cardiac))
write.csv(bargraph_table_cardiac, "bargraph_table_cardiac.csv")
bargraph_table_cardiac<- read.csv("bargraph_table_cardiac.csv", header = TRUE,)
bargraph_table_cardiac$Genotype<- factor(bargraph_table_cardiac$Genotype, levels = c("PGP1", "Q1599X/+", "K1597X/+", "Q1599X/Q1599X", "CHD414/CHD414", "CHD427/CHD427"))
bargraph_table_cardiac$Gene<- factor(bargraph_table_cardiac$Gene)


bargraph_neural<- c("RORA", "SOX3", "NDRG4")
bargraph_table_neural<-as.data.frame(combined_normalized_counts_cellline[rownames(combined_normalized_counts_cellline) %in% bargraph_neural,])
bargraph_table_neural<- as.data.frame(t(bargraph_table_neural))
write.csv(bargraph_table_neural, "bargraph_table_neural.csv")
bargraph_table_neural<- read.csv("bargraph_table_neural.csv", header = TRUE,)
bargraph_table_neural$Genotype<- factor(bargraph_table_neural$Genotype, levels = c("PGP1", "Q1599X/+", "K1597X/+", "Q1599X/Q1599X", "CHD414/CHD414", "CHD427/CHD427"))


#Make a boxplot make graph comparing PGP1 vs LOF mutants
theme_set(theme_gray(base_size = 18))
#Cardiac
p10<- ggplot(bargraph_table_cardiac, aes(x= Gene, y = Counts, fill = Genotype )) + 
  geom_boxplot()
p10<- p10 + scale_fill_manual(values=c("#636363", "#3182bd", "#de2d26","#feb24c", "#2ca25f", "#dd1c77")) + ggtitle("mRNA Levels of Cardiac Genes") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                                           panel.background = element_blank(), axis.line = element_line(colour = "black" )) 
p10<- p10 +theme(plot.title = element_text(hjust = 0.5))
ggsave("2018_03_15_mRNA_levels_of_cardiac_genes.tiff", plot = p10 ,  device = "tiff")

#Neural
p11<- ggplot(bargraph_table_neural, aes(x= Gene, y = Counts, fill = Genotype )) + 
  geom_boxplot()
p11<- p11 + scale_fill_manual(values=c("#636363", "#3182bd", "#de2d26","#feb24c", "#2ca25f", "#dd1c77")) + ggtitle("mRNA Levels of Neural Genes") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                                          panel.background = element_blank(), axis.line = element_line(colour = "black" )) 
p11<- p11 +theme(plot.title = element_text(hjust = 0.5))
ggsave("2018_03_15_mRNA_levels_of_neural_genes.tiff", plot = p11 ,  device = "tiff")

Radhika<- c("ALPK3", "ALPK2", "ALPK1", "TNNT2")
Radhika_Table<-as.data.frame(combined_normalized_counts_cellline[rownames(combined_normalized_counts_cellline) %in% Radhika,])
Radhika_Table<- Radhika_Table[,12:14]
write.csv(Radhika_Table, file = "PGP1_D0_Counts.csv")


#2018_04_24
#Try to make a heatmap of RNA expression of CHD7 direct target genes unique to the het
#PGP1 vs K1597X/+ direct target genes n=67
CHD7_direct_targets_hetonly<- read.table("2018_04_23_CHD7_directtargets_uniquetohet.txt")
CHD7_direct_targets_hetonly<- as.character(CHD7_direct_targets_hetonly[,1])

matrix3 <- assay(rld)[CHD7_direct_targets_hetonly, ]
matrix3 <- matrix3 - rowMeans(matrix3)
#reorder columns because you arent clustering
matrix3<- as.data.frame(matrix3)
matrix3 <- matrix3[,c(12,13,14,1,4,5,10,3,8,9,11,2,6,7,18,19,20,15,16,17)]
#rename the columns
colnames(matrix3)<- c("PGP1_1","PGP1_2","PGP1_2","Q1599X/+_1","Q1599X/+_2","Q1599X/+_3","Q1599X/+_4","K1597X/+_1","K1597X/+_2","K1597X/+_3","K1597X/+_4","Q1599X/Q1599X_1","Q1599X/Q1599X_2","Q1599X/Q1599X_3","CHD414/CHD414_1","CHD414/CHD414_2","CHD414/CHD414_3","CHD427/CHD427_1","CHD427/CHD427_2","CHD427/CHD427_3")

#Reorder the rows based on hhe and hbe expression
hhe_only<- c("PCDH7", "NUDCD1", "COL14A1", "SSBP2", "TRDN", "SYN3")
hbe_only<- c("DDX18","ZIC1", "SPIRE1", "ZNF275", "C8orf37")
hhe_hbe<-c("PSMD14","LSM3","BRD8","CSTF2T", "SPON1", "ERH", "ADAMTS9", "PDGFRB", "RPL38")  
which(row.names(matrix3)== "ZNF275")
matrix3 <- matrix3[c(8,19,20,48,53,66,3,7,34,55,56,4,6,15,22,24,31,45,50,62,1,2,5,9,10:14,16:18,21,23,25:30,32,33,35:44,46,47,49,51,52,54,57:61,63,64,65,67),]
#Add hhe and hbe expression categories
category<-c(rep("HHE only", 6), rep("HBE only", 5), rep("HHE + HBE", 9), rep("None",47))
matrix3df<-data.frame(row.names = row.names(matrix3), category=category)

annotation_col<- as.data.frame(metadata$genotype2)
row.names(annotation_col)<- row.names(metadata)
colnames(annotation_col)<- "genotype"
pheatmap(matrix3, cluster_cols = F, cluster_rows = F, annotation_row = matrix3df, gaps_row= c(6,11,20), cellheight = 10, cellwidth=20, file="2018_04_24_CHD7_direct_targets_hetonly_heatmap.jpg")
dev.off()

#Make a heat map of direct target genes that APPEAR in the het (n=17)
CHD7_direct_targets_hetonly_appear<- read.table("sigheatappear_peaksgenes.txt")
CHD7_direct_targets_hetonly_appear<- as.character(CHD7_direct_targets_hetonly_appear[,1])

matrix4 <- assay(rld)[CHD7_direct_targets_hetonly_appear, ]
matrix4 <- matrix4 - rowMeans(matrix4)
#reorder columns because you arent clustering
matrix4<- as.data.frame(matrix4)
matrix4 <- matrix4[,c(12,13,14,1,4,5,10,3,8,9,11,2,6,7,18,19,20,15,16,17)]
#rename the columns
colnames(matrix4)<- c("PGP1_1","PGP1_2","PGP1_2","Q1599X/+_1","Q1599X/+_2","Q1599X/+_3","Q1599X/+_4","K1597X/+_1","K1597X/+_2","K1597X/+_3","K1597X/+_4","Q1599X/Q1599X_1","Q1599X/Q1599X_2","Q1599X/Q1599X_3","CHD414/CHD414_1","CHD414/CHD414_2","CHD414/CHD414_3","CHD427/CHD427_1","CHD427/CHD427_2","CHD427/CHD427_3")

#Reorder the rows based on hhe and hbe expression
hhe_only_appear<- "C20orf194"
#hbe_only_appear<- none
hhe_hbe_appear<-c("TXNRD1", "CCND2")  
which(row.names(matrix4)== "CCND2")
matrix4 <- matrix4[c(17,7,14,1:6,8:13,15,16),]
#Add hhe and hbe expression categories
category_appear<-c(rep("HHE only", 1), rep("HHE + HBE", 2), rep("None",14))
matrix4df<-data.frame(row.names = row.names(matrix4), category=category_appear)

annotation_col<- as.data.frame(metadata$genotype2)
row.names(annotation_col)<- row.names(metadata)
colnames(annotation_col)<- "genotype"
pheatmap(matrix4, cluster_cols = F, cluster_rows = F, annotation_row = matrix4df, gaps_row= c(1,3), cellheight = 10, cellwidth=20, file="2018_04_24_CHD7_direct_targets_hetappear_heatmap.jpg")
dev.off()


#Add in the normal heart and brain expression (2018_05_08_)
expressiondata<- read.csv("mouse_brain_rnaseq3.csv")
heartexpressed<- subset(expressiondata, e14.5_mean >= 10)
brainexpressed<- subset(expressiondata, Brain_E9.5.norm. >= 10)
heartexpressedgenes<- as.character(unique(heartexpressed$human.External.Gene.Name))
brainexpressedgenes<- as.character(unique(brainexpressed$human.External.Gene.Name))

#Try to make a heatmap of RNA expression of CHD7 direct target genes unique to the het
#PGP1 vs K1597X/+ direct target genes n=67
CHD7_direct_targets_hetonly<- read.table("2018_04_23_CHD7_directtargets_uniquetohet.txt")
CHD7_direct_targets_hetonly<- as.character(CHD7_direct_targets_hetonly[,1])

matrix5 <- assay(rld)[CHD7_direct_targets_hetonly, ]
matrix5 <- matrix5 - rowMeans(matrix5)
#reorder columns because you arent clustering
matrix5<- as.data.frame(matrix5)
matrix5 <- matrix5[,c(12,13,14,1,4,5,10,3,8,9,11,2,6,7,18,19,20,15,16,17)]
#rename the columns
colnames(matrix5)<- c("PGP1_1","PGP1_2","PGP1_3","Q1599X/+_1","Q1599X/+_2","Q1599X/+_3","Q1599X/+_4","K1597X/+_1","K1597X/+_2","K1597X/+_3","K1597X/+_4","Q1599X/Q1599X_1","Q1599X/Q1599X_2","Q1599X/Q1599X_3","CHD414/CHD414_1","CHD414/CHD414_2","CHD414/CHD414_3","CHD427/CHD427_1","CHD427/CHD427_2","CHD427/CHD427_3")

#Reorder the rows based on heart and brain expression (reads > 10)
he_directtargets<- CHD7_direct_targets_hetonly[CHD7_direct_targets_hetonly %in% heartexpressedgenes]
hb_directtargets<- CHD7_direct_targets_hetonly[CHD7_direct_targets_hetonly %in% brainexpressedgenes]
he_be_directtargets<- he_directtargets[he_directtargets %in% hb_directtargets]
he_directtargets_unique<- setdiff(he_directtargets,hb_directtargets)
hb_directtargets_unique<- setdiff(hb_directtargets,he_directtargets)
noexpression<- Reduce(setdiff, list(A=CHD7_direct_targets_hetonly, B= heartexpressedgenes, C=brainexpressedgenes))

matrix5<- rbind(matrix5[row.names(matrix5) %in% he_directtargets_unique,], matrix5[row.names(matrix5) %in% hb_directtargets_unique,], matrix5[row.names(matrix5) %in% he_be_directtargets,], matrix5[row.names(matrix5) %in% noexpression,] )

#Add hhe and hbe expression categories
category<-c(rep("Heart only", 7), rep("Brain only", 14), rep("Heart + Brain", 31), rep("None",34))
matrix5df<-data.frame(row.names = row.names(matrix5), category=category)
category2<-c(rep("wt", 3,), rep("Q1599X/+", 4), rep("K1597X/+", 4), rep("Q1599X/Q1599X", 3), rep("CHD414/CHD414", 3), rep("CHD427/CHD427", 3))
matrix5df2<-data.frame(row.names = colnames(matrix5), category=category2)

annotation_col<- as.data.frame(metadata$genotype2)
annotation_col <- annotation_col[c(12,13,14,1,4,5,10,3,8,9,11,2,6,7,18,19,20,15,16,17),]
row.names(annotation_col)<- colnames(matrix5)
annotation_col<- annotation_col[,-1]
colnames(matrix5df2)<- "genotype"
pheatmap(matrix5, cluster_cols = F, cluster_rows = F, annotation_row = matrix5df, annotation_col = matrix5df2, gaps_row= c(7,21,52), cellheight = 10, cellwidth=20, file="2018_05_08_CHD7_direct_targets_hetonly_heatmap_heartbrain_test.jpg")
dev.off()






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

You need to make the row names be the gene or transcript name. If you're using an annotated FeatureCounts file, make the following adjustments:
```
row.names(counts_file)<- make.names(counts_file$`symbol`, unique = TRUE)
#Remove the symbol column (in my case, column 9)
counts_file<- counts_file[,-9]
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
# Venn Diagrams
If you want to overlap differential genes for different samples, I recommend using the Venn Diagram Package in R. 

library(gProfileR)

# Running gprofiler to identify enriched processes among significant genes
```
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
```

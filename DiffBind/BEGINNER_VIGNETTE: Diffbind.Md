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
# Perform the analysis:
```
peaks<- dba.analyze(peaks)
plot(peaks, contrast = 1)
peaks.db<- dba.report(peaks)
```
## Plot a PCA
```
dba.plotPCA(peaks,contrast=1, label=DBA_TISSUE)
```
Here you can export any plot by clicking "export" on the plot quadrant (bottom right), or you can directly save it using this code. This code is for a tiff, but can be changed for pdf or jpeg by replacing the relevant text.

```
tiff("filename.tiff")
dba.plotPCA(peaks,contrast=1, label=DBA_TISSUE)
dev.off()
```
## Volcano plot
```
dba.plotVolcano(peaks)
```

## Customizable dot plot
Here I was able to customize the plot a bit using the R "plot" function, which is not DiffBind. To get help for any function type ?plot() and read the vignette that comes up on the right. Here I am plotting the concentration (# reads) versus fold change, and I've colored peaks with a positive fold change (up in WT) in read, and a negative fold change (up in mutant) in blue.

```
plot(peaks.db$Conc, peaks.db$Fold,pch=20,cex=.15,col=ifelse(peaks.db$Fold>=0,"red","blue"))
```

## Make a peak file 
This generates a table that you can view all of the differential peaks, with concentrations, fold changes, p-values and FDR.
```
peaks.db.df<- as.data.frame(peaks.db)
```

## Subset the peak file for up and down
This is another nifty trick in R, subsetting dataframes. The following code makes a new variable in the environment, named "peads.db.PGP1enriched". This variable is going to be a subset of the data frame, where peaks.db.df$Fold>0.
In a dataframe, you can specify a specific column using "$". Here we are looking to subset where the value in the column "Fold" > 0. We can do the same for the other direction.

```
peaks.db.PGP1enriched<- peaks.db.df[peaks.db.df$Fold>0,]
peaks.db.CHD4enriched<- peaks.db.df[peaks.db.df$Fold<0,]
```
Now we can save these tables as text files, and they can either be opened in excel or viewed in R. 
```
write.table(peaks.db.PGP1enriched, file = "PGP1vsCHD4_D30_enriched_H3K4me3_peaks.txt", quote = F, row.names = F, sep = "\t")
write.table(peaks.db.CHD4enriched, file = "CHD4vsPGP1_D30_enriched_H3K4me3_peaks.txt", quote = F, row.names = F)

```
To view this table (or any variable in R) in R, simply type:
```
View(peaks.db.PGP1enriched)
``` 
#Next thing we want to do is just get a bed file containing the chr start and stop, separate from this results table. To do this, we want to isolate three columns from the text file. 

```
write.table(peaks.db.PGP1enriched[,1:3], file = "PGP1vsCHD4_D30_enriched_H3K4me3_peaks.bed", quote = F, row.names = F, sep = "\t", col.names = F)
```

Here we are writing a table of the table peaks.db.PGP1enriched, but we are only taking the first three columns 
peaks.db.PGP1enriched[,1:3]

If we were to take the first three rows, it would look like this
peaks.db.PGP1enriched[1:3,]

The rows are entered before the comma, and the columns are entered after. If you wanted just the first three rows and columns 2-5, it would look like this

peaks.db.PGP1enriched[1:3,2:5]

Lastly, if you want all rows, but columns 1,3,5,and 8, it would look like this:
peaks.db.PGP1enriched[,c(1,3,5,8)]

That's DiffBind. If we want to annotate these peaks, we would need another package, ChIP-seeker. 






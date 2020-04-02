
# ChIP-seeker vignette for Beginners

This is an introductory vignette for ChIP-seeker, which I use to annotate peaks for ChIP-seq and ATAC-seq.
There are many packages that can be used to annotate genes, but I like this one specifically because of how it prioritizes gene annotation.

More detailed information can be found here: http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html

# Install and load Packages
```
#Install (this only has to be done once)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ChIPseeker")

#Load (this must be done every time)
library(ChIPseeker)
```

Next, we need to load the genome that we are using. I use hg19 normally, but this can be changed to hg38 or mm10. The txdb is a UCSC genome package that we need to install and load. I've commented out human and included code for mouse here (https://bioconductor.org/packages/release/data/annotation/html/TxDb.Mmusculus.UCSC.mm10.knownGene.html)

```
#Install and load the genome

#library(TxDb.Hsapiens.UCSC.hg19.knownGene) 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")

library("TxDb.Mmusculus.UCSC.mm10.knownGene")
txdb<- TxDb.Mmusculus.UCSC.mm10.knownGene
```
# Data Analysis

## Establish the bed file
Here you are not loading in the bed file, you are simply pointing to the location of the file. If the file is not in your source directory, please provide the path to the file

```
PGP1_enriched_me3<-"PGP1vsCHD4_D30_enriched_H3K4me3_peaks.bed"
```

## Annotate the bed file
```
peakAnno_PGP1_enriched_me3<- annotatePeak(PGP1_enriched_me3, tssRegion = c(-10000,3000), TxDb = txdb, annoDb = "org.Mm.eg.db")
```

## Make the result a data frame

```
peaktableannotated_PGP1_enriched_me3 <- as.data.frame(as.GRanges(peakAnno_PGP1_enriched_me3))
```

## Write the resultant data frame to a file

```
write.table(peaktableannotated_PGP1_enriched_me3, file = "PGP1vsCHD4_diffbind_me3_DF.txt", sep="\t", quote  = F, row.names = F)
```
Here I am writing the file to a tab-separated text file (sep="\t"). I am also removing quotes from all of the rows (try taking that out and see what happens), and removing the row names. 

You can also write to a .csv file which would be write.csv(). For more information, run ?write.csv()

## Get genomic localization and save it as a pdf. 

If you want to see it before saving it, just run the middle line.

```
pdf("Piechartfilename.pdf")
plotAnnoPie(peakAnno_PGP1_enriched_me3)
dev.off()
```
Now, this is a lot of pies on the pie chart. Jon hates them, so he makes me combine all the promoter ones from 1kb to 10kb, so I do that like this:

```
#Pie Chart with Percentages
#Promoter<1, Promoter1-10, 5'UTR, 3'UTR, Exons, 1st Intron, Other Intron, Downstream, Distal Intergenic
pct <- c(20.7, 9.1, 0.16, 1, 1.63, 9.49, 23.76, 0.78, 33.39)
lbls<- c("", "", "", "", "", "", "", "", "")
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # add % to labels
palette<- c("#77d9f9", "#489e40", "#db7208", "#cc956e", "#c592ef", "#6b3993", "#dbd0b6", "#dda523","#7a421d")
pdf("Newpiechartfilename.pdf")
pie(pct,labels = lbls, col=palette)
dev.off()
```

This is not a ChIP-seeker command, this is just an R command "pie". This is a fully customizable pie chart. I got the numbers to put in pct from the first pie chart above. Combined so that the slices of the pie are what you see in the comment. I wound up having blank labels but you can put them in the quotations above. You can also customize the colors.

## Promoter heatmap
```
peak_PGP1_enriched_me3<- readPeakFile(peakfile = PGP1_enriched_me3)
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak_PGP1_enriched_me3, windows=promoter)

pdf("PGP1vsCHD4_diffbind_me3.pdf")
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
dev.off()
```

## Get the name of the genes and assign them to a variable that you can use later.
```
PGP1vsCHD4genes<- unique(as.character(peaktableannotated_PGP1_enriched_me3$SYMBOL))
```

That's it! There are a lot of other interesting charts that ChIP-seeker can make but I normally just stick to these, especially for ChIP-seq of histone marks and ATAC-seq. I have used other charts for transcription factor ChIP-seq but I use these the most.



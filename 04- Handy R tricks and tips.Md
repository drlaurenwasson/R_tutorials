# Handy R tricks and tips

These are some R tips that are not specific to any package, but more for R in general.

## 1. Manipulating data frames

The data frame is likely the most useful thing about R, it can handle large data sets very easily. 


### Rows and columns

Say we have a data frame of an annotated bed file called ```peaks.db.PGP1enriched```. 
First, to view the data frame in a separate R window: 

```
View(peaks.db.PGP1enriched)
```

To see row and column names
```
row.names(peaks.db.PGP1enriched)
colnames(peaks.db.PGP1enriched)
```

#### Isolating rows and columns 

To take all the rows, but the first three columns of a data frame and assign it to a new variable:
```
newdf<- peaks.db.PGP1enriched[,1:3]
```
To take all the columns, but the first three rows of a data frame and assign it to a new variable:

```
newdf<- peaks.db.PGP1enriched[1:3,]
```

The rows are entered before the comma, and the columns are entered after. If you wanted just the first three rows and columns 2-5, it would look like this

```
newdf<- peaks.db.PGP1enriched[1:3,2:5]
```

Lastly, if you want all rows, but columns 1,3,5,and 8, it would look like this:
```
newdf<- peaks.db.PGP1enriched[,c(1,3,5,8)]
```

### Subsetting by certain values in a dataframe:
Coming soon!

## 2. Overlapping lists
This is very, very useful for when you want to overlap two gene lists, say a ChIP-seq list of genes and an RNA-seq list of genes. For this to work, both gene lists need to be character vectors. In the environment, it should say chr if it's a character vector. 

### Make a column a character vector
Say we annotated our gene list and now have a data frame called df, and the gene names are in a column called SYMBOL.
to get a character vector use this command:
```
ChIP_sig_genes<- unique(as.character(df$SYMBOL))
```
This takes all of the gene names in the df column "SYMBOL", then makes them a character vector, then removes duplicates, and stores them in a vector called sig_genes

### Loading in a gene list
Now say, you received a list of significant genes from me, in a file called RNA_sig_genes.txt. It's just a list of genes like this:

```TBX5
MEF2C
SMYD1
```

You will need to read in this list of genes. To do that:
```
RNA_sig_genes<- read.table("RNA_sig_genes.txt", header = F)
```
If you look in yiour environment you would see RNA_sig_genes: x obs of 1 variable. The default of read.table() is to import into a data frame. If you were to view that data frame using View(RNA_sig_genes), you would see that it is one column with the header V1. 
To convert this into a character vector:
```
RNA_sig_genes<- RNA_sig_genes$V1
```
Now if you were to view it in your environment, it's no longer a data frame, so it moved from "Data" to "Values" in the environment. You would then see that it is now a character vector

### Now you are ready to overlap your lists

to see which of the CHIP-seq genes overlap with RNA-seq genes:
```
ChIP_sig_genes[ChIP_sig_genes %in% RNA_sig_genes]
```

to see which of the RNA-seq genes overlap with ChIP-seq genes:
```
RNA_sig_genes[RNA_sig_genes %in% ChIP_sig_genes]
```
This will just print them out on the console. To assign them to a new variable called overlap:

```
overlap<- ChIP_sig_genes[ChIP_sig_genes %in% RNA_sig_genes]
``` 

If you just want to know how many overlap:
```
length(ChIP_sig_genes[ChIP_sig_genes %in% RNA_sig_genes])
```

If you want the ones that do not overlap and assign them to a variable called unique:
In R, the ! is used for "Not" 
```
unique<- ChIP_sig_genes[!ChIP_sig_genes %in% RNA_sig_genes]
```


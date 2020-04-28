# Handy R tricks and tips

These are some R tips that are not specific to any package, but more for R in general.
## 1. Find and set your working directory

The working directory is where all of your files will load from or where they will go when you save them. It's very important to know and set your working directory. When you start a new project in R, it will ask what directory you want to generate the project in, and it will set that as your working directory. If you are coming back after the first time, and dont know what your current working directory is, or if you want to change it, you have to know what the working directory is.

To see what your working directory is:
```
getwd()
```

To change your working directory:
In R-studio at the top, (where you see File, Edit, Code), you should see an option called Session. From here: Session->Set Working Directory -> To File Panes Location

or you can browse by Session->Set Working Directory-> Choose Directory. Then you can pick any directory you want. When you save your R-scripts or export plots, it will save everything in the working directory.


## 2. Load in a file (or a data frame)

Many times, you will need to load a file in from Excel. This can be done in multiple ways, but the easiest is to convert your Excel file to a .csv file. I try to avoid using Excel as much as possible, it messes with the formatting of gene names, but sometimes it is unavoidable. 

Open the file in excel and File-> Save As. Save the document as a .csv file instead of .xlsx. This only works if there is one sheet in your excel file. If there are multiple, export the sheet you want to a new file. Right-Click the tab on the bottom with the sheet name -> Move or Copy -> Select new book from the dropdown menu (check "Make a copy" to avoid moving the file). This will open your sheet of interest in its own book, which can be saved as a .csv file.

A file must be in your working directory in order to call it (or, you must link to the path, but it's easiest to have it in the working directory).

to load a csv and assign it to a variable called df (dataframe):
```
df<- read.csv("File_Name.csv")
```
If your csv has headers (Sample names, etc) =, use:

```
df<- read.csv("File_Name.csv", header = TRUE)
```
This tells R that the first row of your file is different than the other rows (that they are headers) and will automatically make them your column names.

## 3. Manipulating data frames

The data frame is likely the most useful thing about R, it can handle large data sets very easily. 

### Rows and columns

Say we have a data frame called ```df```. 
First, to view the data frame in a separate R window: 

```
View(df)
```

To see row and column names
```
row.names(df)
colnames(df)
```

#### Isolating rows and columns 

To take all the rows, but the first three columns of a data frame and assign it to a new variable:
```
newdf<- df[,1:3]
```
To take all the columns, but the first three rows of a data frame and assign it to a new variable:

```
newdf<- df[1:3,]
```
The rows are entered before the comma, and the columns are entered after. If you wanted just the first three rows and columns 2-5, it would look like this

```
newdf<- df[1:3,2:5]
```

Lastly, if you want all rows, but columns 1,3,5,and 8, it would look like this:
```
newdf<- df[,c(1,3,5,8)]
```
Because you are not removing sequential numbers, you would have to contain them in a list using c()

If you want to REMOVE rows you would stick a "minus" symbol in front. Here you are taking all the data MINUS the ones you specify. Some examples:

To remove the first row:
```
newdf<- df[-1,]
```
To remove the first column:
```
newdf<- df[,-1]
```

To remove rows 1:10: (Note, you have to put minus in front of both, otherwise you will get an error)
```
newdf<- df[-1:-10,]
```

To remove rows 1:10 and columns 4:6
```
newdf<- df[-1:-10,-4:-6]
```

Lastly, to remove columns 1,3,5,and 8, but keep all rows, it would look like this:
```
newdf<- df[,-c(1,3,5,8)]
```

### Subsetting by certain values in a dataframe:
Lets say we have a dataframe (df), which has columns for fold change (fc) and p-value (pval). We want to subset for certain fold change or p-value. This is building upon the scenarios above. 

To select a certain row of a dataframe:
```
newdf<- df[1,]
```
which would take only the first row (a "1" before the comma), but take all the columns (no value after the comma). This is building on the same idea. So you want to take all of the columns for specific rows, so you put your subsetting in brackets before the comma

newdf<- df[<INSERTPARAMETERSHERE>,]

Now, you have two parameters you want to break it down by. fold change and p value. First, fold change: You want to make sure that you take the absolute value of the fold change to get both up and down (or you can make two tables for up and down separately).  If you were just subsetting for fold change, taking all the columns for each row where abs(p) > 2, that would look like:
```
newdf<- df[abs(df$fc) > 2,]
```
If you want to subset just by p-value:
```
newdf<- df[df$pval < 0.05,]
```

To subset by BOTH, you combine the two commands with & (AND). Both parameters go before the comma, since we're subsetting rows by specific parameters.
```
newdf<- df[abs(df$fc) > 2 & df$pval < 0.05,]
```

## 4. Overlapping lists
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



# BioMart in R for beginners

This is an introductory vignette for bioMart. BioMart has a web version which can be found here: https://m.ensembl.org/info/data/biomart/index.html

However, I find it clunky and difficult to use, and much prefer the flexibility of the R-package.

More detailed information can be found here: https://bioconductor.org/packages/release/bioc/html/biomaRt.html

# Install and load Packages
```
#Install (this only has to be done once)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")

#Load (this must be done every time)
library(biomaRt)
```
# Data Analysis

## Load your mart
A "mart" is a database of information for a particular species. 

```
hmart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
```
Other datasets
dataset = 'xtropicalis_gene_ensembl'
dataset = 'mmusculus_gene_ensembl'


Biomart has a few parameters that are important to understand
## Filters 
Filters are the type of input that you are putting into biomart. This could be something like a list of genes, a region of the genome, a chromosome, as examples. Filters have to be inputted from the list of filters for your specific biomart. To obtain what those are we use the following code

```
# Get a list of possible filters
hfilters<- listFilters(hmart)
View(hfilters)
```
## Attributes
Attributes are the type of output that you want bioMart to collect for you. This could be gene names, homology, strand information, etc. Attributes also have be input from the list of attributes for your specific biomart. To obtain those:

```
#get a list of possible attributes
hattributes<- listAttributes(hmart)
View(hattributes)
```




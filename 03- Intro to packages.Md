<b> Introduction to packages </b><br>

Packages are the backbone of R. Packages are modules that are installed and loaded in R to do a variety of functions. <br>
Packages can be found on github or bioconductor. Google is your friend.

Packages only need to be installed once, but must be loaded every time you start R. 
Let's start by installing "Diffbind". 

Google tells me Diffbind is a bioconductor package:
https://bioconductor.org/packages/release/bioc/html/DiffBind.html

To install:
`if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DiffBind")`

This will likely install other packages that are needed. Sometimes you will get an error saying that you have to install other packages. The error messages are pretty straightforward and you can google and install all the packages you need.

To load the package:

`library("DiffBind")`


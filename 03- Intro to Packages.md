# Introduction to Packages

Packages are the backbone of R. Packages are modules that are installed and loaded in R to do a variety of functions.
Packages can be found on github or bioconductor. Google is your friend.

Packages only need to be installed once, but must be loaded every time you start R. Let's start by installing "Diffbind".

Google tells me Diffbind is a bioconductor package: https://bioconductor.org/packages/release/bioc/html/DiffBind.html

To install:


```R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DiffBind")
```

    Bioconductor version 3.10 (BiocManager 1.30.10), R 3.6.1 (2019-07-05)
    Installing package(s) 'DiffBind'
    

    package 'DiffBind' successfully unpacked and MD5 sums checked
    
    The downloaded binary packages are in
    	C:\Users\Lauren\AppData\Local\Temp\Rtmp4qnWMq\downloaded_packages
    

    Old packages: 'backports', 'batchtools', 'biomaRt', 'boot', 'dbplyr',
      'DelayedArray', 'digest', 'edgeR', 'ellipsis', 'fs', 'GenomeInfoDb',
      'ggplot2', 'ggraph', 'ggrepel', 'glue', 'gtools', 'haven', 'Hmisc', 'locfit',
      'matrixStats', 'pillar', 'plyr', 'Rcpp', 'RcppArmadillo', 'RCurl',
      'reshape2', 'rlang', 'S4Vectors', 'scales', 'tibble', 'tidyr', 'tidyselect',
      'vctrs', 'withr', 'xfun', 'askpass', 'BH', 'broom', 'callr', 'caret',
      'class', 'cli', 'clipr', 'cluster', 'curl', 'data.table', 'DBI', 'dplyr',
      'evaluate', 'fansi', 'forcats', 'foreach', 'formatR', 'glmnet', 'gower',
      'hexbin', 'hms', 'htmltools', 'htmlwidgets', 'httpuv', 'httr', 'ipred',
      'IRkernel', 'iterators', 'jsonlite', 'KernSmooth', 'knitr', 'later',
      'lattice', 'lava', 'lubridate', 'markdown', 'MASS', 'Matrix', 'mgcv', 'mime',
      'ModelMetrics', 'modelr', 'nlme', 'nnet', 'numDeriv', 'openssl', 'pkgconfig',
      'prettyunits', 'processx', 'prodlim', 'progress', 'promises', 'ps', 'purrr',
      'quantmod', 'R6', 'recipes', 'repr', 'reprex', 'rmarkdown', 'rstudioapi',
      'rvest', 'selectr', 'shiny', 'spatial', 'SQUAREM', 'stringi', 'survival',
      'sys', 'tidyverse', 'tinytex', 'TTR', 'uuid', 'whisker', 'xml2', 'xts',
      'yaml', 'zoo'
    

This will likely install other packages that are needed. Sometimes you will get an error saying that you have to install other packages. The error messages are pretty straightforward and you can google and install all the packages you need.

To load packages:


```R
library("DiffBind")
```

    Loading required package: GenomicRanges
    Loading required package: stats4
    Loading required package: BiocGenerics
    Loading required package: parallel
    
    Attaching package: 'BiocGenerics'
    
    The following objects are masked from 'package:parallel':
    
        clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
        clusterExport, clusterMap, parApply, parCapply, parLapply,
        parLapplyLB, parRapply, parSapply, parSapplyLB
    
    The following objects are masked from 'package:stats':
    
        IQR, mad, sd, var, xtabs
    
    The following objects are masked from 'package:base':
    
        anyDuplicated, append, as.data.frame, basename, cbind, colnames,
        dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
        grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
        order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
        rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
        union, unique, unsplit, which, which.max, which.min
    
    Loading required package: S4Vectors
    Warning message:
    "package 'S4Vectors' was built under R version 3.6.2"
    Attaching package: 'S4Vectors'
    
    The following object is masked from 'package:base':
    
        expand.grid
    
    Loading required package: IRanges
    Warning message:
    "package 'IRanges' was built under R version 3.6.2"
    Attaching package: 'IRanges'
    
    The following object is masked from 'package:grDevices':
    
        windows
    
    Loading required package: GenomeInfoDb
    Loading required package: SummarizedExperiment
    Warning message:
    "package 'SummarizedExperiment' was built under R version 3.6.2"Loading required package: Biobase
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    Loading required package: DelayedArray
    Warning message:
    "package 'DelayedArray' was built under R version 3.6.2"Loading required package: matrixStats
    Warning message:
    "package 'matrixStats' was built under R version 3.6.2"
    Attaching package: 'matrixStats'
    
    The following objects are masked from 'package:Biobase':
    
        anyMissing, rowMedians
    
    Loading required package: BiocParallel
    Warning message:
    "package 'BiocParallel' was built under R version 3.6.2"
    Attaching package: 'DelayedArray'
    
    The following objects are masked from 'package:matrixStats':
    
        colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges
    
    The following objects are masked from 'package:base':
    
        aperm, apply, rowsum
    
    
    

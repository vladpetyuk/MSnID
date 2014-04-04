MSnID
=====

utilities for handling MS/MS proteomic identifications

to install `MSnID` copy and run these commands from R prompt:

```r
require("devtools") || install.packages("devtools")
# installing dependencies from GitHub to ensure latest versions
install_github("mzID", "thomasp85", quick=TRUE)
install_github("MSnbase", "lgatto", quick=TRUE)
# installing the MSnID itself
install_github("MSnID", "vladpetyuk", quick=TRUE)
```

download an example c_elegans.R script
```r
library("RCurl")
script.url <- "https://raw2.github.com/vladpetyuk/MSnID/master/demo/c_elegans.R"
script <- getURL(script.url, ssl.verifypeer=0L, followlocation=1L)
writeLines(script, "c_elegans.R")
```


A brief description:
* Input: mzIdentML files with MS/MS search results
* Utilities to explore MS/MS search results and
  assess data confidence metrics (MS/MS match scoring, 
  parent ion mass measurement accuracy, missed cleavages ...)
* Optimization of filtering criteria for MS/MS matches to achieve the most number
  of identifications within a pre-defined FDR limit.
* Convertion to MSnSet object for quantitative analysis

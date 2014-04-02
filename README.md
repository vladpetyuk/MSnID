MSnID
=====

utilities for handling MS/MS proteomic identifications

to install run the `MSnID`:

```r
# add path to Bioconductor repositories
source("http://bioconductor.org/biocLite.R")
options(repos=biocinstallRepos(character()))

install.packages("devtools")
library("devtools")
install_github("MSnID", "vladpetyuk", quick=TRUE)

# installing dependencies from GitHub to ensure latest versions
install_github("MSnbase", "lgatto", quick=TRUE)
install_github("mzID", "thomasp85", quick=TRUE)
```

as an example download and run c_elegans.R script into your working directory
```r
library("RCurl")
script.url <- "https://raw2.github.com/vladpetyuk/MSnID/master/demo/c_elegans.R"
script <- getURL(script.url, ssl.verifypeer=0L, followlocation=1L)
writeLines(script, "c_elegans.R")
```


A brief description:
* Input: mzIdentML files with MS/MS search results
* The package utilities assess the data confidence metrics 
  (MS/MS match scoring, parent ion mass measurement accuracy, missed cleavages ...)
* Optimization of filtering criteria for MS/MS matches to achieve the most number
  of identifications within a pre-defined FDR limit.
* Convertion to MSnSet object for quantitative analysis

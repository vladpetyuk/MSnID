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
install_github("MSnID", "vladpetyuk")
```

A brief description:
* Input: mzIdentML file of LC-MS/MS dataset

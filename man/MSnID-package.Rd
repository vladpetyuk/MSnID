\name{MSnID-package}
\alias{MSnID-package}
\alias{MSnID}
\docType{package}
\title{
MSnID: Utilities for handling MS/MS identifications
}

\description{
Extracts MS/MS ID data from mzIdentML (leveraging mzID package) or text files. After collating the search results from multiple datasets it assesses their identification quality and optimize filtering criteria to achieve the maximal identifications at a user specified false discovery rate. Additional utilities include:
   \enumerate{
      \item{
         post-experimental recalibration of mass measurement accuracy 
      }
      \item{
         assessment of irregular and missed cleavages given the enzyme cleavage pattern 
      }
      \item{
         assessment of false discovery rates at peptide-to-spectrum match, unique peptide and protein levels
      }
      \item{
         leverages brute-force and sophisticated optimization routines (Nelder-Mead and simulated annealing) for finding the filtering criteria that provide the maximum spectrum, peptide or protein identifications while not exceeding a corresponding preset threshold of false discovery rate}
      \item{
         converts the results into MSnSet class object as spectral counting data
      }
   }
}


\details{
\tabular{ll}{
Package: \tab MSnID\cr
Type: \tab Package\cr
Version: \tab 0.0.0\cr
Date: \tab 2013-12-11\cr
License: \tab GPL (>= 2)\cr
}

%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~
}

\author{
   Vladislav A. Petyuk (\email{vladislav.petyuk@pnnl.gov})
}

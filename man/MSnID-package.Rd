\name{MSnID-package}
\alias{MSnID-package}
\docType{package}
\title{
MSnID: Utilities for Handling MS/MS Identifications
}

\description{
    Extracts MS/MS ID data from mzIdentML (leveraging \link{mzID} package) 
    or text files. After collating the search results from multiple datasets 
    it assesses their identification quality and optimize filtering criteria 
    to achieve the maximal identifications at a user specified false discovery 
    rate. Additional utilities include:
    \enumerate{
        \item{
            post-experimental recalibration of mass measurement accuracy 
        }
        \item{
            assessment of irregular and missed cleavages given the 
            enzyme cleavage pattern 
        }
        \item{
            assessment of false discovery rates at peptide-to-spectrum match, 
            unique peptide and protein levels
        }
        \item{
            leverages brute-force and sophisticated optimization routines 
            (Nelder-Mead and simulated annealing) for finding the filtering 
            criteria that provide the maximum spectrum, peptide or protein 
            identifications while not exceeding a corresponding preset 
            threshold of false discovery rate}
        \item{
            converts the results into MSnSet class object as 
            spectral counting data
        }
    }
}


\details{
\tabular{ll}{
Package: \tab MSnID\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2014-04-02\cr
License: \tab Artistic-2.0\cr
}
}

\author{
    Vladislav A. Petyuk (\email{vladislav.petyuk@pnnl.gov})
}

\keyword{package}

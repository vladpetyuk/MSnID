\name{id_quality}
\alias{id_quality}

\title{Identification quality}

\description{
    Reports quality for a given level of 
    identification (spectra, peptide or protein).
}

\usage{
    id_quality(.Object, filter=NULL, level=c("PSM", "Peptide", "Accession"))
}

\arguments{
    \item{.Object}{An instance of class "MSnID".}
    \item{filter}{
        Optional argument.
        Either an instance of \link{MSnIDFilter} class or a \code{"character"}.
    }
    \item{level}{
        Level at which the filter will be evaluated. Possible values are
        "PSM", "Peptide" and "Accession". Default value is "PSM".
    }
}

\value{
    Returns a list with with elements "n" and "fdr". 
    Element "n" contains the number of features (spectra, 
    peptides or proteins/accessions) passing the filter. 
    Element "fdr" is the false discovery rate (i.e. identification
    confidence) for the corresponding features.
}

\author{
    Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}

\seealso{
    \code{\link{MSnID}}
    \code{\link{evaluate_filter}}
}


\examples{
data(c_elegans)
id_quality(msnidObj, level="peptide")
id_quality(msnidObj, filter="`MS-GF:PepQValue` < 0.01", level="peptide")
}




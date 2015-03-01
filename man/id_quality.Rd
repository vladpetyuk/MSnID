\name{id_quality}
\alias{id_quality}

\title{Identification quality}

\description{
    Reports quality for a given level of
    identification (spectra, peptide or protein).
}

\usage{
    id_quality(object, filter=NULL, level=c("PSM", "peptide", "accession"))
}

\arguments{
    \item{object}{An instance of class "MSnID".}
    \item{filter}{
        Optional argument.
        Either an instance of \link{MSnIDFilter} class or a \code{"character"}.
    }
    \item{level}{
        Level at which the filter will be evaluated. Possible values are
        "PSM", "peptide" and "accession". Multiple are OK. Default is all of them.
    }
}

\value{
    Returns a matrix with with column names "fdr" and "n".
    Column "n" contains the number of features (spectra,
    peptides or proteins/accessions) passing the filter.
    Column "fdr" is the false discovery rate (i.e. identification
    confidence) for the corresponding features.
    Row names correspond to the provided levels.
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




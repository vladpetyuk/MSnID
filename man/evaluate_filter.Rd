\name{evaluate_filter}
\alias{evaluate_filter}

\title{Filters the MS/MS identifications}

\description{
    Filter out peptide-to-spectrum MS/MS identifications.
}

\usage{
    evaluate_filter(.Object, filter, level=c("PSM", "Peptide", "Accession"))
}

%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{.Object}{An instance of class "MSnID".}
    \item{filter}{
        Either an instance of \link{MSnIDFilter} class or a \code{"character"}.
    }
    \item{level}{
        Level at which the filter will be evaluated. Possible values are
        "PSM", "Peptide" and "Accession".
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
    \code{\link{id_quality}}
}


\examples{
data(c_elegans)

## Filtering using string:
msnidObj <- assess_termini(msnidObj, validCleavagePattern="[KR]\\\\.[^P]")
table(msnidObj$numIrregCleavages)
evaluate_filter(msnidObj, "numIrregCleavages == 0")

## Filtering using filter object:
# first adding columns that will be used as filters
msnidObj$msmsScore <- -log10(msnidObj$`MS-GF:SpecEValue`)
msnidObj$mzError <- abs(msnidObj$experimentalMassToCharge - 
                        msnidObj$calculatedMassToCharge)
# setting up filter object
filtObj <- MSnIDFilter(msnidObj)
filtObj$msmsScore <- list(comparison=">", threshold=10.0)
filtObj$mzError <- list(comparison="<", threshold=0.1) # 0.1 Thomson
show(filtObj)
evaluate_filter(msnidObj, filtObj)
}

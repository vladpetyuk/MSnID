\name{evaluate_filter}
\alias{evaluate_filter}

\title{Filters the MS/MS identifications}

\description{
    Filter out peptide-to-spectrum MS/MS identifications.
}

\usage{
    evaluate_filter(object, filter, level=c("PSM", "peptide", "accession"))
}

%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{An instance of class "MSnID".}
    \item{filter}{
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

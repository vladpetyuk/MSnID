\name{apply_filter}
\alias{apply_filter}
\alias{apply_filter,MSnID,character-method}
\alias{apply_filter,MSnID,MSnIDFilter-method}


\title{Filters the MS/MS identifications}

\description{
    Filter out peptide-to-spectrum MS/MS identifications.
}

\usage{
    apply_filter(msnidObj, filterObj)
}

\arguments{
    \item{msnidObj}{An instance of class "MSnID".}
    \item{filterObj}{
        Either an instance of \link[MSnID]{MSnIDFilter} 
        class or a \code{"character"}.
    }
}

\details{
    \code{filterObj} argument
    evaluated to a \code{"logical"} for each entry of the
    MS/MS results table.
}

\value{
    Returns an instance of "MSnID" class with with peptide-to-spectrum
    matches that pass criteria defined in \code{filterObj} argument.
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

## Filtering using string:
msnidObj <- assess_termini(msnidObj, validCleavagePattern="[KR]\\\\.[^P]")
table(msnidObj$numIrregCleavages)
# getting rid of any other peptides except fully tryptic
msnidObj <- apply_filter(msnidObj, "numIrregCleavages == 0")
show(msnidObj)

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
# applying filter and comparing MSnID object before and after
show(msnidObj)
msnidObj <- apply_filter(msnidObj, filtObj)
show(msnidObj)
}

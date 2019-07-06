\name{correct_peak_selection}
\alias{correct_peak_selection}


\title{Corrects wrong selection of monoisotopic peak}

\description{
    In a typical setting instruments select ions for fragmentation primarily
    based on ion intensity. For low molecular weight peptides the most
    intense peak usually corresponds to monoisotopic peak 
    (that is only C12 carbon isotopes). With increase of molecular weight,
    instensity of monoisotopic peak becomes smaller relatively to heavier 
    peptide isotopes (that is containing one or a few C13 isotopes).

    The function subtracts or adds the mass difference between
    C13 and C12 isotopes (1.0033548378 Da) if that reduces the mass error. 
    Such a mass error arises from the fact that instrument may peak 
    non-monoisotopic peak for fragmentation and thus report the mass 
    that is different by ~ 1 Da.
}

\usage{
    correct_peak_selection(object)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{An instance of class "MSnID".}
}

\value{
    Returns an instance of "MSnID" class with updated 
    \code{experimentalMassToCharge} value.
}

\author{
    Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}

\seealso{
    \code{\link{MSnID}}
    \code{\link{recalibrate}}
    \code{\link{mass_measurement_error}}
}


\examples{
data(c_elegans)
# plot original mass error
massErr <- (msnidObj$experimentalMassToCharge - 
            msnidObj$calculatedMassToCharge) * 
            msnidObj$chargeState
hist(massErr,xlim=c(-1,+1), breaks=seq(-1.5,+1.5,0.01))
# fixing the problem of picking wrong monoisotopic peak
msnidObj <- correct_peak_selection(msnidObj)
# plot fixed mass error
massErr <- (msnidObj$experimentalMassToCharge - 
            msnidObj$calculatedMassToCharge) * 
            msnidObj$chargeState
hist(massErr,xlim=c(-1,+1), breaks=seq(-1.5,+1.5,0.01))
}


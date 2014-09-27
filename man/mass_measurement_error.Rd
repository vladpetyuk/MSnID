\name{mass_measurement_error}
\alias{mass_measurement_error}

\title{Computes error of the parent ion mass to charge measurement}

\description{
    Computes error of the parent ion mass to charge measurement from 
    \code{experimentalMassToCharge} and \code{calculatedMassToCharge}. 
    The returned value is in points per million (ppm).
}

\usage{
    mass_measurement_error(object)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{
        An instance of class "MSnID".
    }
}

\details{
    It may be more common to compute "mass measurement error". 
    However, the practical difference in "mass measurement error" and
    "mass to charge measurement error" is negligible. Moreover, the
    instruments measure mass/charge ratio, not mass \emph{per se}.
}

\value{
    Returns mass to charge measurement error in "ppm" units.
}

\author{
    Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}

\seealso{
    \code{\link{MSnID}}
    \code{\link{recalibrate}}
    \code{\link{correct_peak_selection}}
}


\examples{
data(c_elegans)
ppm <- mass_measurement_error(msnidObj)
hist(ppm, 100)
}


\name{psms}
\alias{psms}
\alias{psms<-}

\title{Peptide-to-spectrum matches}

\description{
    Returns results of MS/MS search (peptide-to-spectrum)
    matches in the form of \code{data.frame}.
}

\usage{
    psms(object, ...)
    psms(object) <- value
}

%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{An instance of class "MSnID".}
    \item{value}{Value is a \code{data.frame} with MS/MS search results}
    \item{\dots}{ignored for now}
}

\value{
    \code{data.frame}
}

\author{
    Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}

\seealso{
    \code{\link{MSnID}}
}

\examples{
data(c_elegans)
msnidDF <- psms(msnidObj)
head(msnidDF)
}

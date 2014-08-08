\name{peptides}
\alias{peptides}

\title{Non-redundant list of peptides}

\description{
Returns the non-redundant list of peptides from the
\emph{MSnID} object
}

\usage{
    peptides(object)
}

\arguments{
    \item{object}{An instance of class "MSnID".}
}

\value{
    Non-redundant list of peptides.
}

\author{
    Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}

\seealso{
    \code{\link{accessions}}
    \code{\link{proteins}}
}


\examples{
data(c_elegans)
head(peptides(msnidObj))
}


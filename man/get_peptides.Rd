\name{get_peptides}
\alias{get_peptides}

\title{Non-redundant list of peptides}

\description{
Returns the non-redundant list of peptides from the
\emph{MSnID} object
}

\usage{
    get_peptides(.Object)
}

\arguments{
    \item{.Object}{An instance of class "MSnID".}
}

\value{
    Non-redundant list of peptides.
}

\author{
    Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}

\seealso{
    \code{\link{get_accessions}}
}


\examples{
data(c_elegans)
head(get_peptides(msnidObj))
}


\name{get_accessions}
\alias{get_accessions}

\title{Non-redundant list of accession (protein) identifiers}

\description{
Returns the non-redundant list of accession (protein) identifiers from the
\emph{MSnID} object
}

\usage{
    get_accessions(.Object)
}

\arguments{
    \item{.Object}{
        An instance of class "MSnID".
    }
}

\value{
    Non-redundant list of accession (protein) identifiers.
}

\author{
    Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}

\seealso{
    \code{\link{get_peptides}}
}


\examples{
data(c_elegans)
head(get_accessions(msnidObj))
}

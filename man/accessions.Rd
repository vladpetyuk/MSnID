\name{accessions}
\alias{accessions}
\alias{proteins}

\title{Non-redundant list of accession (protein) identifiers}

\description{
Returns the non-redundant list of accession (protein) identifiers from the
\emph{MSnID} object. Most of the times accessions and proteins 
have the same meaning. However, there are cases, for example use of 6-frame 
stop-to-stop translation as FASTA file, where the entries are called with
general term accessions rather then proteins. Currently, accessions and 
proteins have the same meaning in \code{MSnID}.
}

\usage{
    accessions(object, ...)
    proteins(object, ...)
}

\arguments{
    \item{object}{
        An instance of class "MSnID".
    }
    \item{...}{
        ignored parameters
    }

}

\value{
    Non-redundant list of accession (protein) identifiers.
}

\author{
    Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}

\seealso{
    \code{\link{peptides}}
}


\examples{
data(c_elegans)
head(accessions(msnidObj))
head(proteins(msnidObj))
}

\name{report_mods}
\alias{report_mods}
\alias{report_mods,MSnID-method}
\title{Lists modification masses in the MSnID object}

\description{
Parses out masses from the \code{modification} column and return them as 
\code{table} with counts.
}

\usage{
    report_mods(object, ...)
}

\arguments{
    \item{object}{
        An instance of class "MSnID".
    }
    \item{...}{
        currently not used
    }

}

\value{
    Counts of each modification mass listed in \code{modification column}.
}

\author{
    Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}

\seealso{
    \code{\link{add_mod_symbol}}
    \code{\link{map_mod_sites}}
}

\examples{
msnidObj <- MSnID(".")
mzids <- system.file("extdata","phospho.mzid.gz",package="MSnID")
msnidObj <- read_mzIDs(msnidObj, mzids)
# reports the masses and number of their occurences
report_mods(msnidObj)
# clean-up cache
unlink(".Rcache", recursive=TRUE)
}


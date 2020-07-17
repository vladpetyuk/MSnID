
\name{read_mzIDs}
\alias{read_mzIDs}

\title{Populates MS/MS results table from mzIdentML files}

\description{
    Reads mzIdentML files into \code{psms} \code{data.table} slot 
    of \code{object} MSnID instance. The functionality leverage
    \link[mzID]{mzID} package facility. Note, the calls are memoised
    using \link[R.cache]{R.cache} facility. So if the call with the same
    list of files issues again, the results will be read from cache instead of
    re-parsing the mzIdentML files.
}

\usage{
    read_mzIDs(object, mzids, backend)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{An instance of class "MSnID"}
    \item{mzids}{paths to mzIdentML (mzid) files}
    \item{backend}{Package that is leveraged for parsing. 
                Either 'mzID' or 'mzR' corresponding to
                \code{\link[mzID]{mzID-package}} and
                \code{\link[mzR]{mzR-class}} respectively. 
                The 'mzR' parser is much faster, 
                since it is based on C++ code. 'mzID' will be kept in the 
                package for leagacy reasons. Note, the default is 'mzID'.}
}

\details{
    mzIdentML files can be either as is or in gzip compressed form (*.mzid.gz).
}

\value{
    Returns an instance of "MSnID" class with \code{@psms} data.table slot
    populated with MS/MS identifications.
}

\author{
    Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}

\seealso{
    \code{\link[mzID]{flatten}}
    \code{\link[mzID]{mzID}}
}


\examples{
\dontrun{
msnidObj <- MSnID(".")
mzids <- system.file("extdata","c_elegans.mzid.gz",package="MSnID")
msnidObj <- read_mzIDs(msnidObj, mzids)
# clean up the cache directory
unlink(".Rcache", recursive=TRUE)
}
}

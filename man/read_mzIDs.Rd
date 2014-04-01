
\name{read_mzIDs}
\alias{read_mzIDs}

\title{Populates MS/MS results table from mzIdentML files}

\description{
  Reads mzIdentML files into \code{psms} \code{data.table} slot 
  of \code{.Object} MSnID instance. The functionality leverage
  \link[mzID]{mzID} package facility. Note, the calls are memoised
  using \link[R.cache]{R.cache} facility. So if the call with the same
  list of files issues again, the results will be read from cache instead of
  re-parsing the mzIdentML files.
}

\usage{
   read_mzIDs(.Object, mzids)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.Object}{
      An instance of class "MSnID"
   }
  \item{mzids}{
      paths to mzIdentML (mzid) files
   }
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
# see c_elegans demo for details 
demo(c_elegans)
}


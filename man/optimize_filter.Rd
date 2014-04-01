\name{optimize_filter}
\alias{optimize_filter}
\alias{optimize_filter,MSnIDFilter,MSnID-method}

\title{Filter criteria optimization to maximize the number of identifications given the FDR upper threshold}

\description{
  Adjusts parameters in the "MSnIDFilter" instance to achieve the most number of 
  spectra, peptides or proteins/accessions within pre-set FDR upper limit.
}

\usage{
   optimize_filter(.Filter, .Data, fdr.max, method, level, n.iter)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.Filter}{
      An instance of class "MSnIDFilter".
  }
  \item{.Data}{
      An instance of class "MSnID".
  }
  \item{fdr.max}{
      Upper limit on acceptable false discovery rate (FDR).
  }
  \item{method}{
      Optimization method. Possible values are "Grid" or 
      same values as for the \code{method} argument 
      of the \link[stats]{optim} function. Tested and 
      recommended arguments (besides "Grid") 
      of \code{method} are "Nelder-Mead" or "SANN".
  }
  \item{level}{
      Determines at what level to perform optimization. 
      Possible values are "PSM", "Peptides" or "Accession".
  }
  \item{n.iter}{
      For method "Grid" is approxomate number of evaluation point. 
      For "Nelder-Mean" and "SANN" methods see \link[stats]{optim}.
  }
}

\details{
      The "Grid" method is brute-force optimization through evaluation of
      approximately \var{n.iter} combinations of the parameters
      set in the "MSnIDFilter" object. The enumeration of the parameter
      combinations proceeds as follows. The \var{n.iter} number is getting
      split given the dimensionality of the problem (that is the number
      of filter parameters in the "MSnIDFilter" object. For each parameter 
      the evaluation points are equally spaced according to quantiles of the
      parameter distribution. This way we enumerate the grid that has
      more evaluation points in relatively more dense areas.
      
      Note, optimization is computationally expensive. Thus, the 
      \code{optimize_filter} call is memoised using facilities
      from the \link[R.cache]{R.cache} package. Once the same 
      call of \code{optime_filter} function issued second time
      the results will be retrieved from cache.
}

\value{
   Returns an instance of "MSnIDFilter" that is maximized to
   provide the most number of identifications while maintaining 
   a pre-set confidence (FDR).
}

\author{
  Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}

\seealso{
  \code{\link{MSnID}}
  \code{\link{evaluate_filter}}
  \code{\link{apply_filter}}
}


\examples{
# see c_elegans demo for details 
# system.file("demo", "c_elegans.R", package="MSnID")

\dontrun{

demo(c_elegans)
}}



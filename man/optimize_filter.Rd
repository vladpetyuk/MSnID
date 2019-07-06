\name{optimize_filter}
\alias{optimize_filter}
\alias{optimize_filter,MSnIDFilter,MSnID-method}

\title{Filter criteria optimization to maximize the number 
        of identifications given the FDR upper threshold}

\description{
    Adjusts parameters in the "MSnIDFilter" instance 
    to achieve the most number of spectra, peptides or 
    proteins/accessions within pre-set FDR upper limit.
}

\usage{
    optimize_filter(filterObj, msnidObj, fdr.max, method, 
                    level, n.iter, mc.cores=NULL)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{filterObj}{
        An instance of class "MSnIDFilter".
    }
    \item{msnidObj}{
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
        Possible values are "PSM", "peptides" or "accession".
    }
    \item{n.iter}{
        For method "Grid" is approxomate number of evaluation point. 
        For "Nelder-Mean" and "SANN" methods see \link[stats]{optim}.
    }
    \item{mc.cores}{
        Controls the number of enabled cores for parallel processing. 
        Make sense only for "Grid" optimizaton. Default is
        \code{getOption("mc.cores", 2L)}. See \link[parallel]{mclapply}
        for details.
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
data(c_elegans)

# explicitely adding parameters that will be used for data filtering
msnidObj$msmsScore <- -log10(msnidObj$`MS-GF:SpecEValue`)
msnidObj$absParentMassErrorPPM <- abs(mass_measurement_error(msnidObj))

# Setting up filter object
filtObj <- MSnIDFilter(msnidObj)
filtObj$absParentMassErrorPPM <- list(comparison="<", threshold=10.0)
filtObj$msmsScore <- list(comparison=">", threshold=10.0)

system.time({
    filtObj.grid <- optimize_filter(filtObj, msnidObj, fdr.max=0.01, 
                                    method="Grid", level="peptide", n.iter=50)})
show(filtObj.grid)

# Fine tuning. Nelder-Mead optimization.
system.time({
    filtObj.nm <- optimize_filter(filtObj.grid, msnidObj, fdr.max=0.01, 
                                    method="Nelder-Mead", level="peptide", 
                                    n.iter=50)})
show(filtObj.nm)

# Fine tuning. Simulated Annealing optimization.
system.time({
    filtObj.sann <- optimize_filter(filtObj.grid, msnidObj, fdr.max=0.01, 
                                    method="SANN", level="peptide", n.iter=50)})
show(filtObj.sann)
}



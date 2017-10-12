
\name{MSnIDFilter-class}
\docType{class}
\alias{MSnIDFilter-class}
\alias{class:MSnIDFilter}
\alias{MSnIDFilter}



%% methods
\alias{show,MSnIDFilter-method}
\alias{$,MSnIDFilter-method}
\alias{$<-,MSnIDFilter-method}
\alias{as.numeric,MSnIDFilter-method}
\alias{length,MSnIDFilter-method}
\alias{update,MSnIDFilter-method}
\alias{names,MSnIDFilter-method}



\title{
    The "MSnIDFilter" Class for Handling MS/MS Criteria, Relationships and 
    Thresholds for Data Filtration.
}

\description{
    The \code{MSnIDFilter} is a convenience class for manipulating the MS/MS 
    filter for MS/MS results.
}


\section{Objects from the Class}{
    The way to create objects is to call \code{MSnIDFilter}
    constructor function that takes as input the \code{MSnID}
    class instance and (optionally) \code{filterList}. 
}

\section{Slots}{
    \describe{
        \item{\code{MSnIDObj}:}{An instance of class "MSnID".}
        \item{\code{filterList}:}{
            An optional argument. A list with element names corresponding to 
            column names avalable in \code{MSnID} instance. Each element
            contains sub-elements "comparison" and "threshold". "Comparison" is
            one of the relationship operators (e.g. ">") 
            see \link{Comparison} for details. "Threshold" is the corresponding
            parameter value the identification has to be more or less (depending
            on comparison) to pass the filter.
        }
    }
}

\section{Methods}{
    \describe{       

    \item{show}{\code{signature(object="MSnIDFilter")}: \cr 
        Prints \code{MSnIDFilter} object.}

    \item{\code{object$name}, \code{object$name<-value}}{
        Access and set filterList elements.}

    \item{names}{signature(x="MSnIDFilter")\cr
        Returns the names of the criteria.}

    \item{\code{as.numeric}}{\code{signature(x="MSnIDFilter")}\cr
        Converts \code{filterList} into "numeric" vector.
        Vector names are the list element names. 
        Vector values are threshold values.
        Comparison operators are lost.}

    \item{\code{length}}{\code{signature(x="MSnIDFilter")}\cr
        Returns the number of criteria set in the "MSnIDFilter" object.}

    \item{\code{update}}{\code{signature(object="MSnIDFilter", ...)}\cr
        The additional \dots argument is \code{numeric} vector of the same
        length as the number of criteria in \code{MSnIDFilter} object.
        The method update the corresponding thresholds to new provided values.}
    }
}

\author{
    Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}

\seealso{
    \code{\link[MSnbase]{MSnSet}}
    \code{\link[MSnID]{evaluate_filter}}
    \code{\link[MSnID]{apply_filter}}
    \code{\link[MSnID]{optimize_filter}}
}

\examples{
data(c_elegans)

## Filtering using filter object:
# first adding columns that will be used as filters
msnidObj$msmsScore <- -log10(msnidObj$`MS-GF:SpecEValue`)
msnidObj$mzError <- abs(msnidObj$experimentalMassToCharge - 
                        msnidObj$calculatedMassToCharge)
# setting up filter object
filtObj <- MSnIDFilter(msnidObj)
filtObj$msmsScore <- list(comparison=">", threshold=10.0)
filtObj$mzError <- list(comparison="<", threshold=0.1) # 0.1 Thomson
show(filtObj)
# applying filter and comparing MSnID object before and after
show(msnidObj)
msnidObj <- apply_filter(msnidObj, filtObj)
show(msnidObj)
}

\keyword{classes}

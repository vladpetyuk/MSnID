\name{apply_filter}
\alias{apply_filter}
\alias{apply_filter,MSnID,character-method}
\alias{apply_filter,MSnID,MSnIDFilter-method}


\title{Filters the MS/MS identifications}

\description{
  Filter out peptide-to-spectrum MS/MS identifications.
}

\usage{
   apply_filter(.Object, .Filter)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.Object}{
      An instance of class "MSnID".
   }
  \item{.Filter}{
      Either an instance of \link[MSnID]{MSnIDFilter} class or a \code{"character"}.
   }
}

\details{
      \code{.Filter} argument
      evaluated to a \code{"logical"} for each entry of the
      MS/MS results table.
}

\value{
   Returns an instance of "MSnID" class with with peptide-to-spectrum
   matches that pass criteria defined in \code{.Filter} argument.
}

\author{
  Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}

\seealso{
  \code{\link{MSnID}}
  \code{\link{evaluate_filter}}
}


\examples{
# see c_elegans demo for details 
# system.file("demo", "c_elegans.R", package="MSnID")

\dontrun{

demo(c_elegans)
}}

\name{evaluate_filter}
\alias{evaluate_filter}

\title{Filters the MS/MS identifications}

\description{
  Filter out peptide-to-spectrum MS/MS identifications.
}

\usage{
   evaluate_filter(.Object, filter, level=c("PSM", "Peptide", "Accession"))
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.Object}{
      An instance of class "MSnID".
   }
  \item{filter}{
      Either an instance of \link{MSnIDFilter} class or a \code{"character"}.
   }
  \item{level}{
      Level at which the filter will be evaluated. Possible values are
      "PSM", "Peptide" and "Accession".
  }
}

\value{
   Returns a list with with elements "n" and "fdr". 
   Element "n" contains the number of features (spectra, 
   peptides or proteins/accessions) passing the filter. 
   Element "fdr" is the false discovery rate (i.e. identification
   confidence) for the corresponding features.
}

\author{
  Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}

\seealso{
  \code{\link{MSnID}}
  \code{\link{id_quality}}
}


\examples{
# see c_elegans demo for details 
# system.file("demo", "c_elegans.R", package="MSnID")

\dontrun{

demo(c_elegans)
}}

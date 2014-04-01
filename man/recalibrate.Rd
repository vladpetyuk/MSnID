\name{recalibrate}
\alias{recalibrate}

\title{Post-experimental recalibration of observed mass to charge ratios}

\description{
  Mass spectrometry measurements like any other real-worls measurements 
  are prone to systematic errors. 
  Typically they are minimized by instrument calibration prior the analysis.
  Nonetheless, the calibration may drift over time or be affected by
  some adverse factors (temperature or space charge fluctuations).
  
  This function estimates and removes the systematic error from the datasets.
  The side effect is the recalibrated \code{experimentalMassToCharge} values.
}

\usage{
   recalibrate(.Object)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.Object}{
      An instance of class "MSnID".
   }
}

\details{
  Currently it employs a very simple method of zero-centering the 
  histogram of mass measurement errors. In the future it will contain
  more sophisticated recalibration routines.
}

\value{
   "MSnID" class instance with updated \code{experimentalMassToCharge}.
}

\author{
  Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}

\seealso{
  \code{\link{MSnID}}
  \code{\link{mass_measurement_error}}
  \code{\link{correct_peak_selection}}
}

\examples{
# see c_elegans demo for details 
# system.file("demo", "c_elegans.R", package="MSnID")

\dontrun{

demo(c_elegans)
}}

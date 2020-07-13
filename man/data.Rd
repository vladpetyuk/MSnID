\name{msnidObj}
\alias{msnidObj}
\alias{c_elegans}
\docType{data}
\title{Example mzIdenML File and MSnID Object}

\description{
MSnID object from
c_elegans_A_3_1_21Apr10_Draco_10-03-04_msgfplus.mzid.gz dataset
from PeptideAltas repository id PASS00308.
}


\usage{data(c_elegans)}

\examples{
data(c_elegans)
msnidObj


\dontrun{
## How to download the example mzID file from PeptideAltas:
try(setInternet2(FALSE),silent=TRUE)
ftp.loc <- "ftp://PASS00308:PJ5348t@ftp.peptideatlas.org/MSGFPlus_Results/MZID_Files/c_elegans_A_3_1_21Apr10_Draco_10-03-04_msgfplus.mzid.gz"
download.file(ftp.loc, "c_elegans.mzid.gz")}

\dontrun{
## Script for generation of the C. elegans example data:
msnidObj <- MSnID(".")
mzids <- system.file("extdata","c_elegans.mzid.gz",package="MSnID")
msnidObj <- read_mzIDs(msnidObj, mzids)
save(msnidObj, file='c_elegans.RData', compress='xz', compression_level=9)
# MD5 sum for the file is: a7c511a6502a6419127f1e46db48ed92
digest::digest(msnidObj)
# clean up the cache directory
unlink(".Rcache", recursive=TRUE)
}

}

\keyword{datasets}

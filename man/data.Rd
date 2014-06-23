\name{msnidObj}
\alias{msnidObj}
\alias{c_elegans}
\docType{data}
\title{
Example mzIdenML File and MSnID Object
}

\description{
MSnID object from c_elegans_A_3_1_21Apr10_Draco_10-03-04_msgfplus.mzid.gz dataset
from PeptideAltas repository id PASS00308.
}


\usage{
data(c_elegans)
}

\examples{
data(c_elegans)
msnidObj


\dontrun{
## How to download the example mzID file from PeptideAltas:
try(setInternet2(FALSE),silent=TRUE)
ftp.loc <- "ftp://PASS00308:PJ5348t@ftp.peptideatlas.org/MSGFPlus_Results/MZID_Files/c_elegans_A_3_1_21Apr10_Draco_10-03-04_msgfplus.mzid.gz"
download.file(cel.path, "c_elegans.mzid.gz")}

\dontrun{
## Script for generation of the example data:
msnidObj <- MSnID(".")
mzids <- system.file("extdata","c_elegans.mzid.gz",package="MSnID")
msnidObj <- read_mzIDs(msnidObj, mzids)
msnidObj$Accession <- msnidObj$accession
msnidObj$accession <- NULL
msnidObj$isDecoy <- msnidObj$isdecoy
msnidObj$isdecoy <- NULL
msnidObj$calculatedMassToCharge <- msnidObj$calculatedmasstocharge
msnidObj$calculatedmasstocharge <- NULL
msnidObj$experimentalMassToCharge <- msnidObj$experimentalmasstocharge
msnidObj$experimentalmasstocharge <- NULL
msnidObj$chargeState <- msnidObj$chargestate
msnidObj$chargestate <- NULL
msnidObj$spectrumID <- msnidObj$spectrumid
msnidObj$spectrumid <- NULL
msnidObj$Peptide <- paste(msnidObj$pre, msnidObj$pepseq, msnidObj$post, sep='.')
msnidObj$pre <- NULL
msnidObj$post <- NULL
msnidObj$pepseq <- NULL
# msnidObj$rank <- as.integer(msnidObj$rank)
# msnidObj$end <- as.integer(msnidObj$end)
# msnidObj$start <- as.integer(msnidObj$start)
# msnidObj$length <- as.integer(msnidObj$length)
# msnidObj$chargeState <- as.integer(msnidObj$chargeState)
save(msnidObj, file='c_elegans.RData', compress='xz', compression_level=9)
# MD5 sum for actual file is: 71b7ae714a7f83fc5ed91665fa14e8bd
# clean up the cache directory
unlink(".Rcache", recursive=TRUE)}



}
\keyword{datasets}

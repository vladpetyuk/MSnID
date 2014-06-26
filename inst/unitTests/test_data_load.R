library("MSnID")
msnid2 <- MSnID(".")
mzids <- system.file("extdata","c_elegans.mzid.gz",package="MSnID")
msnid2 <- read_mzIDs(msnid2, mzids)
msnid2$Accession <- msnid2$accession
msnid2$accession <- NULL
msnid2$isDecoy <- msnid2$isdecoy
msnid2$isdecoy <- NULL
msnid2$calculatedMassToCharge <- msnid2$calculatedmasstocharge
msnid2$calculatedmasstocharge <- NULL
msnid2$experimentalMassToCharge <- msnid2$experimentalmasstocharge
msnid2$experimentalmasstocharge <- NULL
msnid2$chargeState <- msnid2$chargestate
msnid2$chargestate <- NULL
msnid2$spectrumID <- msnid2$spectrumid
msnid2$spectrumid <- NULL
msnid2$Peptide <- paste(msnid2$pre, msnid2$pepseq, msnid2$post, sep='.')
msnid2$pre <- NULL
msnid2$post <- NULL
msnid2$pepseq <- NULL
# clean up the cache directory
unlink(".Rcache", recursive=TRUE)


# Note, the psms slot is data.table
# the .internal.selfref is session specific.
# So I'll compare psms as data.frame only.

test_data_load <- function() {
   # now, check if it is what it supposed to be
   data(c_elegans)
   checkIdentical(psms(msnidObj), psms(msnid2))
}

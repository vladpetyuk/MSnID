library("MSnID")
msnid2 <- MSnID(".")
mzids <- system.file("extdata","c_elegans.mzid.gz",package="MSnID")
msnid2 <- read_mzIDs(msnid2, mzids)
unlink(".Rcache", recursive=TRUE)


# Note, the psms slot is data.table
# the .internal.selfref is session specific.
# So I'll compare psms as data.frame only.

test_data_load <- function() {
    # now, check if it is what it supposed to be
    data(c_elegans)
    checkIdentical(psms(msnidObj), psms(msnid2))
}

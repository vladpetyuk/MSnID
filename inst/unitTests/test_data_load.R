library("MSnID")
library("digest")

msnid2 <- MSnID(".")
mzids <- system.file("extdata","c_elegans.mzid.gz",package="MSnID")
msnid2 <- read_mzIDs(msnid2, mzids, backend = 'mzID')
unlink(".Rcache", recursive=TRUE)
# the object from mzR parser is a bit different from mzID (at this point)
msnid3 <- read_mzIDs(msnid2, mzids, backend = 'mzR') 
unlink(".Rcache", recursive=TRUE)


# Note, the psms slot is data.table
# the .internal.selfref is session specific.
# So I'll compare psms as data.frame only.

test_data_load <- function() {
    # now, check if it is what it supposed to be
    data(c_elegans)
    checkIdentical(psms(msnidObj), psms(msnid2))
}


test_column_names <- function() {
    # now, check if the column names are the same
    data(c_elegans)
    checkIdentical(names(msnidObj), names(msnid2))
}


test_data_load_mzR <- function() {
    # now, check if it is what it supposed to be
    # checkIdentical(digest(psms(msnid3)),'e5c572c07878673f1165822969f81869')
    # checkIdentical(digest(psms(msnid3)),'0b4e3b61e3fe007ed11651632fa3f1fb')
    checkIdentical(digest(psms(msnid3)),'d1c961e8b3decd00ae7d376ab87af42f')
}

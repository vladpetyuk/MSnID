library("MSnID")
library("digest")
data(c_elegans)
msnidObj <- assess_missed_cleavages(msnidObj, 
                                    missedCleavagePattern="[KR](?=[^P$])")


test_check_name_presence <- function(){
    checkTrue("numMissCleavages" %in% names(msnidObj))
}

test_check_type <- function(){
    checkTrue(is.numeric(msnidObj$numMissCleavages))
}

test_check_result <- function(){
    checkIdentical(digest(msnidObj$numMissCleavages),
                   "5d6f1fbd250a5ee10b0ad9d9c171ed90")
}


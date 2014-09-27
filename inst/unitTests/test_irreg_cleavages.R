library("MSnID")
library("digest")
data(c_elegans)
msnidObj <- assess_termini(msnidObj, validCleavagePattern="[KR]\\.[^P]")

test_check_name_presence <- function(){
    checkTrue("numIrregCleavages" %in% names(msnidObj))
}

test_check_type <- function(){
    checkTrue(is.numeric(msnidObj$numIrregCleavages))
}

test_check_result <- function(){
    checkIdentical(digest(msnidObj$numIrregCleavages),
                   "e0c4ab882ffb0c626c53ab801dc5a7e3")
}



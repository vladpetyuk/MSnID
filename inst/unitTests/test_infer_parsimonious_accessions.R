library("MSnID")
data(c_elegans)


test_infer_parsimonious_accessions_old <- function(){
    # explicitely adding parameters that will be used for data filtering
    msnidObj$msmsScore <- -log10(msnidObj$`MS-GF:SpecEValue`)
    msnidObj$absParentMassErrorPPM <- abs(mass_measurement_error(msnidObj))
    
    # quick-and-dirty filter. The filter is too strong for the sake of saving time
    # at the minimal set of proteins inference step.
    msnidObj <- apply_filter(msnidObj, 'msmsScore > 12 & absParentMassErrorPPM < 2')
    
    msnidObj2 <- infer_parsimonious_accessions(msnidObj)
    checkEqualsNumeric(length(proteins(msnidObj2)), 551)
}



# Above is the old function for testing protein inference.  I'll leave it for
# now.  Below is the new way, where first all the inference will be done
# outside of the test functions.


# explicitely adding parameters that will be used for data filtering
msnidObj$msmsScore <- -log10(msnidObj$`MS-GF:SpecEValue`)
msnidObj$absParentMassErrorPPM <- abs(mass_measurement_error(msnidObj))
# quick-and-dirty filter. The filter is too strong for the sake of saving time
# at the minimal set of proteins inference step.
msnidObj <- apply_filter(msnidObj, 'msmsScore > 12 & absParentMassErrorPPM < 2')

# checking with razor peptide
msnidObj.r <- infer_parsimonious_accessions(msnidObj, unique_only=FALSE)

# checking unique peptides only
msnidObj.u <- infer_parsimonious_accessions(msnidObj, unique_only=TRUE)

test_infer_parsimonious_accessions_number.r <- function(){
   checkEqualsNumeric(length(proteins(msnidObj.r)), 551)
}
test_infer_parsimonious_accessions_number.u <- function(){
    checkEqualsNumeric(length(proteins(msnidObj.u)), 427)
}

test_infer_parsimonious_accessions_hash.r <- function(){
    checkIdentical(digest(psms(msnidObj.r)$accession),
                   'e8ff9bbf36929ce49850733b3640edc8')
}
test_infer_parsimonious_accessions_hash.u <- function(){
    checkIdentical(digest(psms(msnidObj.u)$accession),
                   '2a375d3ec85d022aeba5d210bca55a53')
}

# Future challenges is to come up with tests that check inference that is
# done outside of MSnID object


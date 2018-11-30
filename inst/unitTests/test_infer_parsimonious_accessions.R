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
msnidObj2 <- infer_parsimonious_accessions(msnidObj)

test_infer_parsimonious_accessions_number <- function(){
   checkEqualsNumeric(length(proteins(msnidObj2)), 551)
}

# test_infer_parsimonious_accessions_hash <- function(){
#    checkIdentical(digest(psms(msnidObj2)),'42c967304603b17ef667fae5b8d5657f')
# }

test_infer_parsimonious_accessions_hash <- function(){
    checkIdentical(digest(psms(msnidObj2)$accession),
                   '6a3566a95b2e49a0f966d22ed897a752')
}

# Future challenges is to come up with tests that check inference that is
# done outside of MSnID object


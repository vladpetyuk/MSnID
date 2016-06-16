library("MSnID")
data(c_elegans)

test_infer_parsimonious_accessions <- function(){

        # explicitely adding parameters that will be used for data filtering
    msnidObj$msmsScore <- -log10(msnidObj$`MS-GF:SpecEValue`)
    msnidObj$absParentMassErrorPPM <- abs(mass_measurement_error(msnidObj))
    
    # quick-and-dirty filter. The filter is too strong for the sake of saving time
    # at the minimal set of proteins inference step.
    msnidObj <- apply_filter(msnidObj, 'msmsScore > 12 & absParentMassErrorPPM < 2')
    
    msnidObj2 <- infer_parsimonious_accessions(msnidObj)
    checkEqualsNumeric(length(proteins(msnidObj2)), 551)
}



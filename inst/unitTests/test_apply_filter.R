library("MSnID")
data(c_elegans)

test_apply_filter_1 <- function(){
    # testing string filter
    out <- apply_filter(msnidObj, "`MS-GF:SpecEValue` < 10^-10")
    # cheching against expected values
    
    # these are old values that ignored peptide-to-protein
    # matching redundancy.
    # Left them here just in case.
    # ref <- matrix(c(0.000285103349964362, 7017,
    #                 0.000778816199376947, 2570,
    #                 0.00193986420950533 , 1033),
    #               ncol=2, nrow=3, byrow=TRUE,
    #               dimnames=list(c("PSM","peptide","accession"),
    #                             c("fdr","n")))
    
    # new PSM FDR values that relies on calculation that 
    # removes protein assignments first
    ref <- matrix(c(0.0005005005005005, 3998,
                    0.000778816199376947, 2570,
                    0.00193986420950533 , 1033),
                  ncol=2, nrow=3, byrow=TRUE,
                  dimnames=list(c("PSM","peptide","accession"),
                                c("fdr","n")))
    checkEqualsNumeric(id_quality(out), ref)
}


test_apply_filter_2 <- function(){
    # checking object-based filter

    msnidObj$msmsScore <- -log10(msnidObj$`MS-GF:SpecEValue`)
    msnidObj$mzError <- abs(msnidObj$experimentalMassToCharge -
                            msnidObj$calculatedMassToCharge)
    # setting up filter object
    filtObj <- MSnIDFilter(msnidObj)
    filtObj$msmsScore <- list(comparison=">", threshold=9.0)
    filtObj$mzError <- list(comparison="<", threshold=0.2)
    # appying filter
    out <- apply_filter(msnidObj, filtObj)
    
    # cheching against expected values
    
    # old values (see above)
    # ref <- matrix(c(0.0011079650375477 , 8132,
    #                 0.0013550135501355 , 2956,
    #                 0.00801424755120214, 1132),
    #               ncol=2, nrow=3, byrow=TRUE,
    #               dimnames=list(c("PSM","peptide","accession"),
    #                             c("fdr","n")))
     
    # new values (see above)
    ref <- matrix(c(0.00087260034904014, 4588,
                    0.0013550135501355, 2956,
                    0.00801424755120214, 1132),
                  ncol=2, nrow=3, byrow=TRUE,
                  dimnames=list(c("PSM","peptide","accession"),
                                c("fdr","n")))
    
    checkEqualsNumeric(id_quality(out), ref)
}



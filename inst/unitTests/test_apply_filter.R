library("MSnID")
data(c_elegans)

test_apply_filter_1 <- function(){
    # checking string filter
    out <- apply_filter(msnidObj, "`MS-GF:SpecEValue` < 10^-10")
    checkEquals(id_quality(out, level="PSM"),
                list(fdr=0.000285103349964362, n=7017))
    checkEquals(id_quality(out, level="peptide"),
                list(fdr=0.000778816199376947, n=2570))
    checkEquals(id_quality(out, level="accession"),
                list(fdr=0.00193986420950533, n=1033))
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

    out <- apply_filter(msnidObj, filtObj)

    checkEquals(id_quality(out, level="PSM"),
                list(fdr=0.0011079650375477, n=8132))
    checkEquals(id_quality(out, level="peptide"),
                list(fdr=0.0013550135501355, n=2956))
    checkEquals(id_quality(out, level="accession"),
                list(fdr=0.00801424755120214, n=1132))
}



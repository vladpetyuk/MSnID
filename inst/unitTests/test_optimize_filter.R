library("MSnID")
data(c_elegans)

msnidObj$msmsScore <- -log10(msnidObj$`MS-GF:SpecEValue`)
msnidObj$absParentMassErrorPPM <- abs(mass_measurement_error(msnidObj))

filtObj <- MSnIDFilter(msnidObj)
filtObj$absParentMassErrorPPM <- list(comparison="<", threshold=10.0)
filtObj$msmsScore <- list(comparison=">", threshold=10.0)

set.seed(0)
filtObj.grid <- optimize_filter(filtObj, msnidObj, fdr.max=0.001, 
                                method="Grid", level="peptide", 
                                n.iter=500)

set.seed(0)
filtObj.nm <- optimize_filter(filtObj.grid, msnidObj, fdr.max=0.001, 
                              method="Nelder-Mead", level="peptide", 
                              n.iter=500)

set.seed(0)
filtObj.sann <- optimize_filter(filtObj.grid, msnidObj, fdr.max=0.001, 
                              method="SANN", level="peptide", 
                              n.iter=500)


test_grid_optimization <- function(){
    checkEquals(filtObj.grid@filterList$absParentMassErrorPPM$threshold,
                3.99291882215223)
    checkEquals(filtObj.grid@filterList$msmsScore$threshold,
                9.04628489941093)
}


test_Nelder_Mead_optimization <- function(){
    checkEquals(filtObj.nm@filterList$absParentMassErrorPPM$threshold,
                3.81032565709881)
    checkEquals(filtObj.nm@filterList$msmsScore$threshold,
                8.87912408683179)
}


test_SANN_optimization <- function(){
    checkEquals(filtObj.sann@filterList$absParentMassErrorPPM$threshold,
                3.99800146330237)
    checkEquals(filtObj.sann@filterList$msmsScore$threshold,
                8.8787177674537)
}

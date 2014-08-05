

#--- READING MZID FILE(S) ----------------
library("MSnID")
# start the project by providing work directory
msnid <- MSnID(".")
mzids <- system.file("extdata","c_elegans.mzid.gz",package="MSnID")
msnid <- read_mzIDs(msnid, mzids)
#---------------------------------------


#--- CHECKING WHAT IS INSIDE -----------
head(psms(msnid))
names(msnid)
show(msnid) 
#---------------------------------------


# --- EXTRA INFO ABOUT PEPTIDE SEQUENCES ---
msnid <- assess_termini(msnid, validCleavagePattern="[RK]\\.[^P]")
msnid <- assess_missed_cleavages(msnid, missedCleavagePattern="[KR](?=[^P$])")
# visualize distributions of missed cleavages as an example
pepCleav <- unique(psms(msnid)[,c("numMissCleavages", "isDecoy", "peptide")])
pepCleav <- as.data.frame(table(pepCleav[,c("numMissCleavages", "isDecoy")]))
library("ggplot2")
ggplot(pepCleav, aes(x=numMissCleavages, y=Freq, fill=isDecoy)) + 
    geom_bar(stat="identity", position="dodge") + 
    ggtitle("Number of Missed Cleavages")
# number of cysteins per peptide sequence as an example of sequence analysis
msnid$NumCys <- sapply(lapply(strsplit(msnid$peptide,''),'==','C'),sum)
# calculating peptide lengths
msnid$pepLength <- nchar(msnid$peptide) - 4
pepLen <- unique(psms(msnid)[,c("pepLength", "isDecoy", "peptide")])
# distribution of peptide lengths in the dataset
ggplot(pepLen, aes(x=pepLength, fill=isDecoy)) + 
    geom_histogram(position='dodge', binwidth=3) +
    ggtitle("Distribution on of Peptide Lengths")
#------------------------------------------


# --- TRIM THE DATA -----------------------
# Let's take a look how filtering on irregular and
# missed cleavages will affect the FDR.
show(msnid)
# 1. Leave only fully tryptic
msnid <- apply_filter(msnid, "numIrregCleavages == 0")
show(msnid)
# 2. Retain peptides with at most 2 missed cleavages
msnid <- apply_filter(msnid, "numMissCleavages <= 2")
show(msnid)
#-----------------------------------------


# --- DEALING WITH PARENT ION MASS MEASUREMENT ACCURACY ---
# original mass measurement error in ppm
ppm <- mass_measurement_error(msnid)
ggplot(as.data.frame(ppm), aes(x=ppm)) + 
    geom_histogram(binwidth=100)
# In this particular case the problem that it goes
# far is that there was no certainity that the reported
# masses are the masses of monoisotopic peaks.
# Let's take a look at mass measurement error in Dalton units.
dM <- with(psms(msnid), 
            (experimentalMassToCharge-calculatedMassToCharge)*chargeState)
x <- data.frame(dM, isDecoy=msnid$isDecoy)
ggplot(x, aes(x=dM, fill=isDecoy)) + 
    geom_histogram(position='stack', binwidth=0.1)
# Fixing the peak picking problem
msnid.fixed <- correct_peak_selection(msnid)
# Now the errors are confined within 20 ppm
ppm <- mass_measurement_error(msnid.fixed)
ggplot(as.data.frame(ppm), aes(x=ppm)) + 
    geom_histogram(binwidth=0.25)
# alternatively we can ignore erroneously picked picks
# simply filter the data +/- 20 ppm
msnid.chopped <- apply_filter(msnid, "abs(mass_measurement_error(msnid)) < 20")
ppm <- mass_measurement_error(msnid.chopped)
ggplot(as.data.frame(ppm), aes(x=ppm)) + 
    geom_histogram(binwidth=0.25)
# The data can be recalibrated to remove any 
# systematic components in the mass measurement error.
msnid <- recalibrate(msnid.chopped)
ppm <- mass_measurement_error(msnid)
ggplot(as.data.frame(ppm), aes(x=ppm)) + 
    geom_histogram(binwidth=0.25)
#------------------------------------------------


# ---- MS/MS FILTER ------------------------
# First filtering criteria - MS-GF Spec E-value
msnid$msmsScore <- -log10(msnid$`MS-GF:SpecEValue`)
# Second filtering criteria - absolute mass measurement error
msnid$absParentMassErrorPPM <- abs(mass_measurement_error(msnid))
# visualization
params <- psms(msnid)[,c("msmsScore","absParentMassErrorPPM","isDecoy")]
ggplot(params) + 
    geom_density(aes(x = msmsScore, color = isDecoy, ..count..))
ggplot(params) + 
    geom_density(aes(x = absParentMassErrorPPM, color = isDecoy, ..count..))
# subsetting params to top 10000 to accelerate the plotting
set.seed(0)
params <- params[sample.int(nrow(msnid),10000),]
ggplot(data=params, aes(x=msmsScore, y=absParentMassErrorPPM, color=isDecoy)) +
    geom_point(size=1.5)
# Setting up filter object
filtObj <- MSnIDFilter(msnid)
filtObj$absParentMassErrorPPM <- list(comparison="<", threshold=10.0)
filtObj$msmsScore <- list(comparison=">", threshold=10.0)
# print and visualize filter
show(filtObj)
ggplot(data=params, aes(x=msmsScore, y=absParentMassErrorPPM, color=isDecoy)) +
    geom_point(size=1.5) + 
    geom_hline(yintercept=filtObj$absParentMassErrorPPM$threshold, 
        linetype='dashed') + 
    geom_vline(xintercept=filtObj$msmsScore$threshold, 
        linetype='dashed')
# Evaluate filter. Check how well it performs at different levels.
evaluate_filter(msnid, filtObj, level="PSM")
evaluate_filter(msnid, filtObj, level="peptide")
evaluate_filter(msnid, filtObj, level="accession")
# Brute-force optimization by enumeration the criteria combinations.
# The results should be good starting parameters for follow-up 
# fine tuning optimizations.
system.time({
    filtObj.grid <- optimize_filter(filtObj, msnid, fdr.max=0.01, 
                                    method="Grid", level="peptide", 
                                    n.iter=500)})
show(filtObj.grid)
# Fine tuning. Nelder-Mead optimization.
system.time({
    filtObj.nm <- optimize_filter(filtObj.grid, msnid, fdr.max=0.01, 
                                    method="Nelder-Mead", level="peptide",
                                    n.iter=500)})
show(filtObj.nm)
# visualize filter
ggplot(data=params, aes(x=msmsScore, y=absParentMassErrorPPM, color=isDecoy)) +
    geom_point(size=1.5) + 
    geom_hline(yintercept=filtObj.nm$absParentMassErrorPPM$threshold, 
                linetype='dashed') + 
    geom_vline(xintercept=filtObj.nm$msmsScore$threshold, 
                linetype='dashed')
# compare original and optimized filters
evaluate_filter(msnid, filtObj, level="peptide")
evaluate_filter(msnid, filtObj.nm, level="peptide")
# filtering the MSnID object
msnid <- apply_filter(msnid, filtObj.nm)
show(msnid)
# removing reverse/decoy and Contaminants
msnid <- apply_filter(msnid, "isDecoy == FALSE")
show(msnid)
msnid <- apply_filter(msnid, "!grepl('Contaminant',accession)")
show(msnid)
#---------------------------------------------




# --- CONVERTING TO MSnSet -------------------
msnset <- as(msnid, "MSnSet")
#---------------------------------------------



# --- clean up the cache directory ---
unlink(".Rcache", recursive=TRUE)
#---------------------------------------------

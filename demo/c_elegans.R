# --- DOWNLOADING MS/MS SEARCH RESULTS ---
try(setInternet2(FALSE),silent=TRUE)
ftp <- "ftp://PASS00308:PJ5348t@ftp.peptideatlas.org/"
meta <- read.delim(sprintf("%sSample_Metadata.txt", ftp), as.is=TRUE)
meta <- subset(meta, Age == 'young' & Diet == 'ff') # explore the effect of daf-16
for( dataset in meta$PNNL.Dataset.Name){
   cel.path <- sprintf("%s/MSGFPlus_Results/MZID_Files/%s_msgfplus.mzid.gz", 
                       ftp, dataset)
   download.file(cel.path, sprintf("%s_msgfplus.mzid.gz", dataset))
}
#----------------------------------------


#--- READING MZID FILES ----------------
library("MSnID")
# start the project by providing work directory
msnid <- MSnID(".")
mzids <- list.files(".", pattern=".mzid.gz")
msnid <- read_mzIDs(msnid, mzids)
#---------------------------------------






#--- CHECKING WHAT IS INSIDE -----------
head(psms(msnid))
names(msnid)
show(msnid) # Key columns are not set
#---------------------------------------




# --- UPDATES TO INCLUDE KEY COLUMNS ---
# Note, the anticipated/suggested columns in the
# peptide-to-spectrum matching results are:
#    Peptide, Accession, isDecoy, 
#    calculatedMassToCharge, experimentalMassToCharge, 
#    spectrumFile, spectrumID
msnid$Peptide <- paste(msnid$pre, msnid$pepseq, msnid$post, sep='.')
msnid$Accession <- msnid$accession
msnid$isDecoy <- msnid$isdecoy
msnid$calculatedMassToCharge <- msnid$calculatedmasstocharge
msnid$experimentalMassToCharge <- msnid$experimentalmasstocharge
msnid$chargeState <- msnid$chargestate
msnid$spectrumID <- msnid$spectrumid
# tidying up. removing left-over columns
msnid$accession <- NULL
msnid$pepSeq <- NULL
msnid$pre <- NULL
msnid$post <- NULL
msnid$idecoy <- NULL
msnid$spectrumid <- NULL
msnid$experimentalmasstocharge <- NULL
msnid$calculatedmasstocharge <- NULL
msnid$chargestate <- NULL
#-------------------------------------
show(msnid)


# --- EXTRA INFO ABOUT PEPTIDE SEQUENCES ---
msnid <- assess_termini(msnid, validCleavagePattern="[RK]\\.[^P]")
msnid <- assess_missed_cleavages(msnid, missedCleavagePattern="[KR](?=[^P$])")
# visualize distributions of missed cleavages as an example
library("ggplot2")
pepCleav <- unique(psms(msnid)[,c("NumMissCleavages", "isDecoy", "Peptide")])
pepCleav <- as.data.frame(table(pepCleav[,c("NumMissCleavages", "isDecoy")]))
ggplot(pepCleav, aes(x=NumMissCleavages, y=Freq, fill=isDecoy)) + 
   geom_bar(stat='identity', position='dodge')
#----------------------------------------




# --- TRIM THE DATA -----------------------
show(msnid)
# 1. Leave only fully tryptic
msnid <- apply_filter(msnid, "NumIrregCleavages == 0")
show(msnid)
# 2. Retain peptides with at most 2 missed cleavages
msnid <- apply_filter(msnid, "NumMissCleavages <= 2")
show(msnid)
#-----------------------------------------


# --- CHECKING PARENT MASS MEASUREMENT ACCURACY ---
hist(mass_measurement_error(msnid), 100)
# retain only those PSMs that have 
# parent mass measurement accuracy less then 10 ppm
msnid <- apply_filter(msnid, "abs(mass_measurement_error(msnid)) < 10")
hist(mass_measurement_error(msnid), 100)
# recalibrate parent mass measurement (if necessary)
msnid <- recalibrate(msnid)
hist(mass_measurement_error(msnid), 100)
#-------------------------------------------




# ---- MS/MS FILTER ------------------------
# 1. defining parameters that will be used for filtering
msnid$absParentMassErrorPPM <- abs(mass_measurement_error(msnid))
msnid$msmsScore <- -log10(msnid$`ms-gf:specevalue`)
# 2. visualizing parameter distributions
library("reshape2")
library("ggplot2")
params <- psms(msnid)[,c("msmsScore","absParentMassErrorPPM","isDecoy")]
ggplot(params) +
   geom_density(aes(x = msmsScore, color = isDecoy, ..count..))
ggplot(params) +
   geom_density(aes(x = absParentMassErrorPPM, color = isDecoy, ..count..))
# 2. setting up filter object
filtObj <- MSnIDFilter(msnid)
filtObj$absParentMassErrorPPM <- list(comparison="<", threshold=10.0)
filtObj$msmsScore <- list(comparison=">", threshold=7.0)
# print filter
show(filtObj)
# evaluate filter
evaluate_filter(msnid, filtObj, level="PSM")
evaluate_filter(msnid, filtObj, level="Peptide")
evaluate_filter(msnid, filtObj, level="Accession")
# 3. optimize filter
# brute-force optimization by enumeration all the parameter combinations
# these should be good starting parameters for follow-up optimizations
system.time({
   filtObj.grid <- optimize_filter(filtObj, msnid, fdr.max=0.01, 
                                  method="Grid", level="Peptide", n.iter=1000)})
show(filtObj.grid)
# (absParentMassErrorPPM < 2) & (msmsScore > 7.8) 

# Nelder-Mead optimization
set.seed(0)
system.time({
   filtObj.nm <- optimize_filter(filtObj.grid, msnid, fdr.max=0.01, 
                          method="Nelder-Mead", level="Peptide", n.iter=1000)})
show(filtObj.nm)
# (absParentMassErrorPPM < 3) & (msmsScore > 7.8) 

# simulated annealing optimization
set.seed(0)
system.time({
   filtObj.sann <- optimize_filter(filtObj.grid, msnid, fdr.max=0.01, 
                                  method="SANN", level="Peptide", n.iter=1000)})
show(filtObj.sann)
# (absParentMassErrorPPM < 2.2) & (msmsScore > 7.6)

# check the results
evaluate_filter(msnid, filtObj, level="Peptide")
evaluate_filter(msnid, filtObj.grid, level="Peptide")
evaluate_filter(msnid, filtObj.nm, level="Peptide")
evaluate_filter(msnid, filtObj.sann, level="Peptide")
# 4. filter main msnid object
msnid <- apply_filter(msnid, filtObj.sann)
show(msnid)
# 5. let's remove reverse/decoy and Contaminants
msnid <- apply_filter(msnid, "!isDecoy")
show(msnid)
msnid <- apply_filter(msnid, "!grepl('Contaminant',Accession)")
show(msnid)
#--------------------------------------


# --- CONVERTING TO MSnSet -------------
msnset <- as(msnid, "MSnSet")
# Note, feature data is peptide-centric. Peptide to protein assigments
# stored in feature data.
head(fData(msnset))
# setting up pheno data
# Update sample names. 
# Retain only dataset name portion from spectrum file names, 
# that were used as sample names
sampleNames(msnset) <- sub("(.*)_dta\\.txt","\\1",sampleNames(msnset))
show(msnset)
# prepare pheno data
rownames(meta) <- meta$PNNL.Dataset.Name
meta <- meta[,c("Letter.Replicate","Daf.16.type")]
pData(msnset) <- meta[sampleNames(msnset),]
validObject(msnset)
#---------------------------------------


#--- ROLLING TO PROTEIN LEVEL ----------
# assessing the extent of peptide/protein mapping redundancy problem
redundancy <- table(sapply(fData(msnset)$Accession, length))
barplot(redundancy, main='number of proteins containing peptide sequence')
dim(msnset)
# [1] 6440   10
# summing by ignoring redundancy issue
msnset.prot <- combineFeatures(msnset, fData(msnset)$Accession, 
                                  redundancy.handler="multiple", 
                                  fun="sum", cv=FALSE)
# Combined 10130 features into 1892 using sum
# summing of uniquely matching peptides only
msnset.prot <- combineFeatures(msnset, fData(msnset)$Accession, 
                                  redundancy.handler="unique", 
                                  fun="sum", cv=FALSE)
# Combined 4157 features into 1093 using sum
# We'll pick data from unique peptides for further analysis
#---------------------------------------


#--- SUBSET TO FREQUENTLY PRESENT ------
# at least 6 samples have to have non-zero counts
msnset.prot <- msnset.prot[rowSums(exprs(msnset.prot) > 0) >= 6,]
# 566 proteins left
#---------------------------------------


# --- STATISTICAL TESTS ----------------
if(!require("msmsTests")){
   library("BiocInstaller")
   biocLite("msmsTests")
   library("msmsTests")
}
alt.f <- "y ~ Daf.16.type + 1"
null.f <- "y ~ 1"
div <- colSums(exprs(msnset.prot))
res <- msms.glm.qlll(msnset.prot, alt.f, null.f, div=div)
sum(p.adjust(res$p.value, "BH") < 0.05)
# Post-test filter
lst <- test.results(res,msnset.prot,pData(msnset.prot)$Daf.16.type,"wt","mut",div,
                    alpha=0.05,minSpC=0,minLFC=1,
                    method="BH")
res.volcanoplot(lst$tres, min.LFC=1, max.pval=0.05, ylbls=NULL, maxy=4)
#--------------------------------------





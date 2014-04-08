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
show(msnid) # Key columns are not set yet.
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
show(msnid)
#-------------------------------------


# --- EXTRA INFO ABOUT PEPTIDE SEQUENCES ---
msnid <- assess_termini(msnid, validCleavagePattern="[RK]\\.[^P]")
msnid <- assess_missed_cleavages(msnid, missedCleavagePattern="[KR](?=[^P$])")
# visualize distributions of missed cleavages as an example
pepCleav <- unique(psms(msnid)[,c("NumMissCleavages", "isDecoy", "Peptide")])
pepCleav <- as.data.frame(table(pepCleav[,c("NumMissCleavages", "isDecoy")]))
library("ggplot2")
ggplot(pepCleav, aes(x=NumMissCleavages, y=Freq, fill=isDecoy)) + 
   geom_bar(stat="identity", position="dodge") + 
   ggtitle("Number of Missed Cleavages")
# number of cysteins per peptide sequence as an example of sequence analysis
msnid$NumCys <- sapply(lapply(strsplit(msnid$Peptide,''),'==','C'),sum)
# calculating peptide lengths
msnid$PepLength <- nchar(msnid$Peptide) - 4
pepLen <- unique(psms(msnid)[,c("PepLength", "isDecoy", "Peptide")])
# distribution of peptide lengths in the dataset
ggplot(pepLen, aes(x=PepLength, fill=isDecoy)) + 
   geom_histogram(position='dodge', binwidth=3) +
   ggtitle("Distribution on of Peptide Lengths")
#------------------------------------------


# --- TRIM THE DATA -----------------------
# Let's take a look how filtering on irregular and
# missed cleavages will affect the FDR.
show(msnid)
# 1. Leave only fully tryptic
msnid <- apply_filter(msnid, "NumIrregCleavages == 0")
show(msnid)
# 2. Retain peptides with at most 2 missed cleavages
msnid <- apply_filter(msnid, "NumMissCleavages <= 2")
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
dM <- with(psms(msnid), (experimentalMassToCharge-calculatedMassToCharge)*chargeState)
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
msnid$msmsScore <- -log10(msnid$`ms-gf:specevalue`)
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
evaluate_filter(msnid, filtObj, level="Peptide")
evaluate_filter(msnid, filtObj, level="Accession")
# Brute-force optimization by enumeration the criteria combinations.
# The results should be good starting parameters for follow-up 
# fine tuning optimizations.
system.time({
   filtObj.grid <- optimize_filter(filtObj, msnid, fdr.max=0.01, 
                                   method="Grid", level="Peptide", n.iter=500)})
show(filtObj.grid)
# Fine tuning. Nelder-Mead optimization.
system.time({
   filtObj.nm <- optimize_filter(filtObj.grid, msnid, fdr.max=0.01, 
                                 method="Nelder-Mead", level="Peptide", n.iter=500)})
show(filtObj.nm)
# visualize filter
ggplot(data=params, aes(x=msmsScore, y=absParentMassErrorPPM, color=isDecoy)) +
   geom_point(size=1.5) + 
   geom_hline(yintercept=filtObj.nm$absParentMassErrorPPM$threshold, 
              linetype='dashed') + 
   geom_vline(xintercept=filtObj.nm$msmsScore$threshold, 
              linetype='dashed')
# compare original and optimized filters
evaluate_filter(msnid, filtObj, level="Peptide")
evaluate_filter(msnid, filtObj.nm, level="Peptide")
# filtering the MSnID object
msnid <- apply_filter(msnid, filtObj.nm)
show(msnid)
# removing reverse/decoy and Contaminants
msnid <- apply_filter(msnid, "!isDecoy")
show(msnid)
msnid <- apply_filter(msnid, "!grepl('Contaminant',Accession)")
show(msnid)
#---------------------------------------------


# msnid$msgfpepqvalue <- msnid$`ms-gf:pepqvalue`
# ff <- MSnIDFilter(msnid)
# ff$msgfpepqvalue <- list(comparison="<", threshold=0.01)
# evaluate_filter(msnid, ff, level="Peptide")


# --- CONVERTING TO MSnSet -------------------
msnset <- as(msnid, "MSnSet")
# Note, feature data is peptide-centric. Peptide to protein assigments
# stored in feature data.
head(fData(msnset))
# Note, sample names in msnset are based on file names that were used as input
# MS/MS search engine. Let's trim the file names to make them compatible with
# dataset names.
head(sampleNames(msnset))
head(meta)
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
#-------------------------------------------







#--- ROLLING TO PROTEIN LEVEL ----------
# assessing the extent of peptide/protein mapping redundancy problem
redundancy <- table(sapply(fData(msnset)$Accession, length), dnn="redundancy")
redundancy <- 100 * prop.table(redundancy)
ggplot(as.data.frame(redundancy), aes(x=factor(1), y=Freq, fill=redundancy)) +
   geom_bar(stat='identity', width=1) + 
   coord_polar(theta='y') + 
   xlab('') + ylab('') +
   labs(fill='Redundancy') +
   scale_x_discrete(breaks = NULL)
# original number of features/peptides in the data
length(featureNames(msnset))
# summing of uniquely matching peptides only
msnset.prot <- combineFeatures(msnset, fData(msnset)$Accession, 
                                  redundancy.handler="unique", 
                                  fun="sum", cv=FALSE)
# subset proteins
# at least 6 samples must have non-zero counts
msnset.prot <- msnset.prot[rowSums(exprs(msnset.prot) > 0) >= 6,]
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


<<<<<<< HEAD
# --- work in progress ---
# selecting regulated only
=======

# --- HEATMAP ---
if(!require("Heatplus")){
   library("BiocInstaller")
   biocLite("Heatplus")
   library("Heatplus")
}
>>>>>>> 96dc4166f307d747cdb7d9d9d4b2867c4cd40b45
regulated <- subset(lst$tres, adjp < 0.05 & abs(LogFC) > 1)
# order MSnSet object the daf-16 status
msnset.prot <- msnset.prot[,order(pData(msnset.prot)$Daf.16.type)]
# matrix with regulated proteins
selected.data <- exprs(msnset.prot[rownames(regulated),])
# more meaningful sample names
colnames(selected.data) <- with(pData(msnset.prot), 
                                paste(Daf.16.type, Letter.Replicate, sep='.'))
# scaling counts from 0 to 1
selected.data <- sweep(selected.data, 1, apply(selected.data, 1, min), '-')
selected.data <- sweep(selected.data, 1, apply(selected.data, 1, max), '/')
<<<<<<< HEAD
group.colors <- c('green','red')[as.factor(pData(msnset.prot)$Daf.16.type)]
library("gplots")
heatmap.2(selected.data,
          Colv=FALSE,
          dendrogram='row',
          col=colorRampPalette(c("snow","steelblue"))(10),
          ColSideColors=group.colors,
          key=FALSE,
          trace='none',
          main='Scaled Spectral Counts for ')







=======
heatmap_plus(selected.data,
             scale='none',
             col=colorRampPalette(c("snow","steelblue"))(10))
>>>>>>> 96dc4166f307d747cdb7d9d9d4b2867c4cd40b45

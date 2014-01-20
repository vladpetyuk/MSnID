

flatten <- function(object, no.redundancy=FALSE) {
      flatPSM <- mzID:::flatten(object@psm)
      flatPSM <- flatPSM[, colnames(flatPSM) != 'id']
      flatEviData <- 
         cbind(object@evidence@evidence,
               object@database@database[
                  match(object@evidence@evidence$dBSequence_ref,
                        object@database@database$id), ])
      flatEviData <- flatEviData[,!names(flatEviData) == 'id']
      flatPep <- mzID:::flatten(object@peptides)
      flatPepEviData <- 
         merge( flatPep, flatEviData, 
                by.x="id", by.y="peptide_ref", all=TRUE)
      if(no.redundancy){
         flatPepEviData <- 
            flatPepEviData[!duplicated(flatPepEviData[,'id']),]
      }
      flatAll <- merge(flatPSM, flatPepEviData, 
                       by.x='peptide_ref', by.y='id', all=TRUE)
      flatAll$spectrumFile <- 
         object@parameters@rawFile$name[
            match(flatAll$spectraData_ref,
                  object@parameters@rawFile$id)]
      flatAll$databaseFile <- 
         object@parameters@databaseFile$name[
            match(flatAll$searchDatabase_ref,
                  object@parameters@databaseFile$id)]
      flatAll <- flatAll[, !grepl('_ref$', 
                                  names(flatAll), 
                                  perl=T) & 
                            !names(flatAll) == 'id']
      return(flatAll)
   }


# non-exported function
# extract_mzID <- function(mzIdentMLFilePath,
#                          scoreString = "-log10(`MS-GF:SpecEValue`)")
# {
#    x <- mzID(mzIdentMLFilePath)
#    x <- flatten(x) # I think this still reports only one protein
#    # a little bit of tweaking of flatten msms output
#    y <- with(subset(x, !modified), data.frame(
#       spectrumID = sub("index=","",spectrumID),
#       parentMassErrorPPM = 1e6*(calculatedMassToCharge - 
#                                    experimentalMassToCharge)/calculatedMassToCharge,
#       msmsScore=eval(parse(text=scoreString)),
#       chargeState,
#       peptide = paste(pre,pepSeq,post,sep='.'),
#       isDecoy, # Note, it is not computed from accession strings
#       accession,
#       spectrumFile))
#    return(y)
# }



# extract_mzID <- function(mzIdentMLFilePath)
# {
#    x <- mzID(mzIdentMLFilePath)
#    x <- flatten(x) # I think this still reports only one protein. will fix later
#    return(x)
# }



.read_mzIDs <- function(mzids)
{
   cl <- makeCluster(detectCores(), outfile='')
   registerDoParallel(cl)
   timing <- system.time(
      psms <- foreach( i=icount(length(mzids)), 
                       .combine=rbind,
                       .inorder=FALSE,
                       .packages=c("mzID")) 
      %dopar% {
         # ans <- MSnID:::extract_mzID(mzids[i]) # extract_mzID is non-exported
         ans <- MSnID:::flatten(mzID(mzids[i]))
         print(paste(i, mzids[i], "done", sep=" ... "))
         gc()
         ans})
   print(timing)
   stopCluster(cl)
   return(psms)
}


.read_mzIDs.memoized <- addMemoization(.read_mzIDs)




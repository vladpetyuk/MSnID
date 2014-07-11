
.read_mzIDs.memoized <- function(mzids) {
   # Try to load cached data, if exists
   key <- list(mzids)
   data <- loadCache(key)
   if (!is.null(data)) {
      cat("Loaded cached data\n")
   }else{   
      cat("Reading from mzIdentMLs ...\n")
      data <- data.table(flatten(mzID(mzids)))
      saveCache(data, key=key)
   }
   
   # the piece below will be removed in the future,
   # one the corresponding functionality implemented in mzID
   # do some cleaning/tweaking of column modes and name casing
   
   # Columns that must be present:
   # Peptide, Accession, isDecoy, calculatedMassToCharge, 
   # experimentalMassToCharge, chargeState, spectrumFile, spectrumID
   if(!is.null(data$accession)){
      data$Accession <- data$accession
      data$accession <- NULL
   }

   if(!is.null(data$isdecoy)){
      data$isDecoy <- data$isdecoy
      data$isdecoy <- NULL
   }
   
   if(!is.null(data$calculatedmasstocharge)){
      data$calculatedMassToCharge <- as.numeric(data$calculatedmasstocharge)
      data$calculatedmasstocharge <- NULL
   }
   
   if(!is.null(data$experimentalmasstocharge)){
      data$experimentalMassToCharge <- as.numeric(data$experimentalmasstocharge)
      data$experimentalmasstocharge <- NULL
   }
   
   if(!is.null(data$chargestate)){
      data$chargeState <- data$chargestate
      data$chargestate <- NULL
   }
   
   if(!is.null(data$spectrumid)){
      data$spectrumID <- data$spectrumid
      data$spectrumid <- NULL
   }
   
   if(!is.null(data$pre)){
      data$Pre <- data$pre
      data$pre <- NULL
   }
   
   if(!is.null(data$post)){
      data$Post <- data$post
      data$post <- NULL
   }
   
   if(!is.null(data$pepseq)){
      data$PepSeq <- data$pepseq
      data$pepseq <- NULL
   }
   #---
   data$Peptide <- paste(data$Pre, data$PepSeq, data$Post, sep='.')
   
   return(data)
}





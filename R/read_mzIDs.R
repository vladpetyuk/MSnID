
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
   return(data)
}



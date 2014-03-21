# .read_mzIDs <- function(mzids)
# {
#    numCores <- min(length(mzids), detectCores())
#    cl <- makeCluster(numCores, outfile='')
#    registerDoParallel(cl)
#    timing <- system.time(
#       psms <- foreach( i=icount(length(mzids)), 
#                        .combine=rbind,
#                        .inorder=FALSE,
#                        .packages=c("mzID")) 
#       %dopar% {
#          ans <- data.table(flatten(mzID(mzids[i])))
#          print(paste(i, mzids[i], "done", sep=" ... "))
#          ans})
#    print(timing)
#    stopCluster(cl)
#    return(psms)
# }

# .read_mzIDs <- function(mzids)
# {
#    return(data.table(flatten(mzID(mzids))))
# }
# 
# .read_mzIDs.memoized <- addMemoization(.read_mzIDs)


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



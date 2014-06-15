# Here will be methods for handling filter optimization.
# Since they apply both to MSnID and MSnIDFilter objects 
# it probably make sense not to relate these methods exclusively
# to either class. May be they migrate to MSnIDFilter-methods, though.


.get_num_pep_for_fdr <- function(filterThresholds, msmsdata, filter, fdr.max, ...) 
{
   filter <- update(filter, filterThresholds)
   x <- evaluate_filter(msmsdata, filter, ...) # level should get here through ...
   if(is.nan(x$fdr) || x$fdr > fdr.max){
      return(rnorm(1,sd=0.001)) # 0 is bad because optimization does not move
   }else{
      return(x$n)
   }
}


.construct_optimization_grid <- function(filterObj, msnidObj, n.iter)
{
   #
   # don't really like this -1 hack
   n.iter.per.param <- round(n.iter^(1/length(filterObj))) - 1 
   #
   probs <- seq(0, 1, 1/n.iter.per.param)
   eg <- expand.grid(
      lapply(names(filterObj), 
             function(arg) 
                quantile(msnidObj[[arg]], probs, na.rm=T)))
   colnames(eg) <- names(filterObj)
   return(eg)
}


.optimize_filter.grid <- function(filterObj, msnidObj, fdr.max, level, n.iter)
{
   par.grid <- .construct_optimization_grid(filterObj, msnidObj, n.iter)
   evaluations <- apply(par.grid, 1, .get_num_pep_for_fdr, 
                        msnidObj, filterObj, fdr.max, level)
   optim.pars <- par.grid[which.max(evaluations),]
   newFilter <- update(filterObj, as.numeric(optim.pars))
   return(newFilter)
}


.optimize_filter.grid.mclapply <- function(filterObj, msnidObj, fdr.max, level, n.iter)
{
   par.grid <- .construct_optimization_grid(filterObj, msnidObj, n.iter)
#    browser()
   evaluations <- mclapply(seq_len(nrow(par.grid)), 
                           function(i){
                              .get_num_pep_for_fdr( par.grid[i,],
                                                    msnidObj, filterObj, 
                                                    fdr.max, level)
                           }, 
                           mc.cores=detectCores())
   evaluations <- unlist(evaluations)
   optim.pars <- par.grid[which.max(evaluations),]
   newFilter <- update(filterObj, as.numeric(optim.pars))
   return(newFilter)
}
           
           
setMethod("optimize_filter",
          signature(.Filter="MSnIDFilter", .Data="MSnID"),
          definition=function(.Filter, .Data, fdr.max, method, level, n.iter)
          {
             return(.optimize_filter.memoized(.Filter, .Data, 
                                              fdr.max, method, level, n.iter))
          }
)



.optimize_filter <- function(.Filter, .Data, fdr.max, method, level, n.iter)
{
   method <- match.arg(method, 
                       choices=c("Grid", "Nelder-Mead", "SANN"))
   level <- match.arg(level,
                      choices=c("PSM", "Peptide", "Accession"))
   #
   if(method == "Grid"){
      if(.Platform$OS.type == "unix"){
         optimFilter <- .optimize_filter.grid.mclapply(.Filter, .Data,
                                                       fdr.max, level, n.iter)
      }else{
         # yet to implement effective parallelization on Windows
         optimFilter <- .optimize_filter.grid(.Filter, .Data,
                                              fdr.max, level, n.iter)
      }
   }
   if(method %in% c("Nelder-Mead", "SANN")){
      optim.out <- optim(par=as.numeric(.Filter),
                         fn = .get_num_pep_for_fdr,
                         msmsdata = .Data,
                         filter = .Filter,
                         fdr.max = fdr.max,
                         level = level,
                         method = method,
                         control=list(fnscale=-1, maxit=n.iter))
      optimFilter <- update(.Filter, optim.out$par)
   }
   return(optimFilter)
}   
.optimize_filter.memoized <- addMemoization(.optimize_filter)
# .optimize_filter.memoized <- .optimize_filter



# -- addMemoization
# memArgs <- list(...)
# function(..., envir = parent.frame()) {
#    args <- list(fcn, ..., envir = envir)
#    args <- c(args, memArgs)
#    do.call("memoizedCall", args = args)
# }

# -- memoizedCall
# key <- list(what = what, ...)
# if (!force) {
#    res <- loadCache(key = key, dirs = dirs, sources = sources)
#    if (!is.null(res)) {
#       if (verbose) 
#          cat("Returning cached results!")
#       return(res)
#    }
# }
# res <- do.call(what, args = list(...), quote = FALSE, envir = envir) # <- perhaps problem here!
# saveCache(res, key = key, dirs = dirs, sources = sources)
# res

#===========================================================


# .read_mzIDs <- function(mzids)
# {
#    cl <- makeCluster(detectCores(), outfile='')
#    registerDoParallel(cl)
#    timing <- system.time(
#       psms <- foreach( i=icount(length(mzids)), 
#                        .combine=rbind,
#                        .inorder=FALSE,
#                        .packages=c("mzID")) 
#       %dopar% {
#          # ans <- MSnID:::extract_mzID(mzids[i]) # extract_mzID is non-exported
#          ans <- flatten(mzID(mzids[i]))
#          print(paste(i, mzids[i], "done", sep=" ... "))
#          gc()
#          ans})
#    print(timing)
#    stopCluster(cl)
#    return(psms)
# }
# 
# 
# .read_mzIDs.memoized <- addMemoization(.read_mzIDs)
# 
# setMethod("read_mzIDs", "MSnID",
#           definition=function(.Object, pathToMzIDs)
#           {
#              # to check if files are indeed available by the provided path
#              stopifnot(all(sapply(pathToMzIDs, file.exists)))
#              .Object@psms <- .read_mzIDs.memoized(pathToMzIDs)
#              return(.Object)
#           }
# )



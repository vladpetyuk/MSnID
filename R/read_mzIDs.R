
utils::globalVariables(c("i", "spectrumID", "name", "mass", "location", 
                         "modification", "DatabaseAccess", 
                         "DatabaseDescription", "DBseqLength"))


.read_mzIDs.memoized <- function(mzids)
{
    # Try to load cached data, if exists
    key <- list(mzids)
    data <- loadCache(key)
    if (!is.null(data)){
        message("Loaded cached data\n")
    }else{
        message("Reading from mzIdentMLs ...\n")
        data <- data.table(flatten(mzID(mzids), safeNames=FALSE))
        saveCache(data, key=key)
    }

    # Adding peptide identification in a conventional format:
    # X.XXXXX.X, where pre and post amino acids separated by dots
    # from the main sequence
    # This may be substituted in the future with more flexible
    # peptide S4 object
    data$peptide <- paste(data$pre, data$pepSeq, data$post, sep='.')
    
    # trimming spectrumFile databaseFile down to basenames
    # mzID 1.5.3 reported them differently on POSIX vs Win
    # so that brings it to common denominator
    data$spectrumFile <- basename(gsub('\\\\','/',data$spectrumFile))
    data$databaseFile <- basename(gsub('\\\\','/',data$databaseFile))

    # Columns that must be present:
    # peptide, accession, isDecoy, calculatedMassToCharge,
    # experimentalMassToCharge, chargeState, spectrumFile, spectrumID
    # There is no checking for these columns at this point.
    # If the mzIdentML file conformed with PSI standard, most likely
    # those columns are present.
    # There are checks for columns in some of the other functions
    # that rely on presence of particular columns.

    return(data)
}




factor_to_str_converter <- function(df){
    data.frame(lapply(df, function(x){
        if(is.factor(x))
            x <- as.character(x)
        return(x)}), 
        stringsAsFactors=FALSE)
}



.read_mzIDs.mzR.engine.single.file <- function(mzid){
    mzRidentObj <- openIDfile(mzid)
    x.psms <- psms(mzRidentObj) %>% factor_to_str_converter
    x.scor <- score(mzRidentObj) %>% factor_to_str_converter
    x.mods <- modifications(mzRidentObj) %>% factor_to_str_converter
    x.mods <- group_by(x.mods, spectrumID, sequence, peptideRef) %>%
        summarise(modification = paste(mass,' (',location,')',sep='',collapse=', ')) %>%
        select(spectrumID,sequence,peptideRef,modification)
    #' merging
    stopifnot(all(as.character(x.psms$spectrumID) == as.character(x.scor$spectrumID)))
    x <- cbind(x.psms, x.scor[,setdiff(colnames(x.scor),'spectrumID')])
    x <- left_join(x, x.mods, by=c("spectrumID", "sequence", "peptideRef"))
    x$modified <- ifelse(is.na(x$modification), FALSE, TRUE)
    x$spectrumFile <- fileName(mzRidentObj) # very redundant. not good
    x <- rename(x,
                accession = DatabaseAccess,
                description = DatabaseDescription,
                length = DBseqLength,
                pepSeq = sequence,
                peptideRef = peptideRef)
    x <- data.table(x, safeNames=FALSE)
    return(x)
}


.read_mzIDs.mzR <- function(mzids){
    if(length(mzids) == 1){
        res <- .read_mzIDs.mzR.engine.single.file(mzids)
    }
    else {
        nCores <- detectCores()
        nThreads <- ifelse(length(mzids) < nCores, length(mzids), nCores)
        cl <- makeCluster(nThreads, outfile = '')
        on.exit(stopCluster(cl))
        registerDoParallel(cl)
        res <- foreach(i = icount(length(mzids)),
                       .packages = c("mzR",'dplyr','data.table'),
                       .export=c(".read_mzIDs.mzR.engine.single.file",
                                 'factor_to_str_converter')) %dopar% 
            {
                cat("reading ", basename(mzids[i]), "...\n", sep = "")
                res.i <- .read_mzIDs.mzR.engine.single.file(mzids[i])
                cat(basename(mzids[i]), "DONE!\n")
                res.i
            }
        res <- rbindlist(res)
    }
    return(res)
}

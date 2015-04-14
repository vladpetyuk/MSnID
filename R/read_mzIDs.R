
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
    data$spectrumFile <- basename(data$spectrumFile)
    data$databaseFile <- basename(data$databaseFile)

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





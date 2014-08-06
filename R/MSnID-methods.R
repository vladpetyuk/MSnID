

#----Peptide Sequence Handling--------------------------------------------------

.get_clean_peptide_sequence <- function(peptide){
    return(
        .strip_flanking_AAs(
            .strip_modifications_from_peptide(peptide)))
}

.strip_modifications_from_peptide <- function(peptide)
{
    # INPUT: Peptide sequence with flanking AAs e.g. "K.APEP*TID{34.34}E%.-"
    # OUTPUT: sequence with just AA symbols e.g. "K.APEPTIDE.-"
    # clean from the mods
    # make sure it removes "." in the middle of the sequence
    aminoAcids <- "ARNDCEQGHILKMFPSTWYV"
    # http://en.wikipedia.org/wiki/Proteinogenic_amino_acid
    extendedAAs <- "BUZXO"
    # aminoAcids <- paste(LETTERS,collapse='') # hardcore solution
    aminoAcids <- paste(aminoAcids, extendedAAs, sep='')
    nonPeptideSequenceCharacter <- sprintf("[^%s.-]", aminoAcids)
    peptide <- gsub( nonPeptideSequenceCharacter, '', peptide)
    # stip "." in the middle of the sequence
    # This wont touch the dots denoting cleavages, assuming (!) there is
    # only one amino acid shown after the cleavage point.
    peptide <- gsub( "(?<=.{2})\\.(?=.{2})", '', peptide, perl=TRUE)
    return(peptide)
}

.strip_flanking_AAs <- function(peptide)
{
    # INPUT: Peptide sequence with flanking AAs e.g. "K.APEP*TID{34.34}E%.-"
    # OUPUT: Peptide without flanking "APEP*TID{34.34}E%"
    # 
    # this is quicker, but may not be safe if there are some mods with dots
    # stopifnot(nchar(peptide) - nchar(gsub("\\.","",peptide)) == 2)
    # this one is safer, but likely to be slower
    stopifnot(all(grepl("^.\\..+\\..$", peptide)))
    #
    # remove flanking AAs
    peptide <- substring( peptide, 3, nchar(peptide)-2)
    #
    return(peptide)
}


.assess_missed_cleavages <- function(peptide,
                                        missedCleavagePattern)
{
    # the fastest approach
    peptide <- .get_clean_peptide_sequence(peptide)
    # now no flanking AAs, no mods
    numberOfMissedCleavages <- 
        nchar(peptide) - 
        nchar(gsub(missedCleavagePattern, "", peptide, perl=TRUE))
    return( numberOfMissedCleavages )
}


setMethod(
    "assess_missed_cleavages", 
    signature("MSnID"),
    definition=function(.Object, missedCleavagePattern)
    {
        .check_column(.Object, "peptide")
        .Object@psms$numMissCleavages <-
            .assess_missed_cleavages(as.character(.Object@psms$peptide),
                                    missedCleavagePattern)
        return(.Object)
    }
)

#---

.assess_termini <- function(peptide, 
                            validCleavagePattern)
{
    # "[RK]\\.[^P]" | "-\\."
    # It has to be sequence with flanking AAs e.g. K.XXXXXXXR.X
    # first peptide needs to be cleaned on anything other then 
    # 20 amino acids and . and -
    #
    # make sure all have flanked AAs (this line of code is redundant)
    stopifnot(all(grepl("^.\\..+\\..$", peptide)))
    #
    # strip mods?
    peptide <- .strip_modifications_from_peptide(peptide)
    #
    is.N.valid <- grepl(sprintf("^%s", validCleavagePattern), peptide)
    is.C.valid <- grepl(sprintf("%s$", validCleavagePattern), peptide)
    is.protein.N.terminus <- grepl("^-\\.",peptide)
    is.protein.C.terminus <- grepl("\\.-$",peptide)
    number.of.irregular.termini <- 
        as.numeric(!(is.N.valid | is.protein.N.terminus)) + # valid N-term
        as.numeric(!(is.C.valid | is.protein.C.terminus))   # valid C-term
    return(number.of.irregular.termini)
}



setMethod(
    "assess_termini", 
    signature("MSnID"),
    definition=function(.Object, validCleavagePattern)
    {
        .Object@psms$numIrregCleavages <- 
            .assess_termini(as.character(.Object@psms$peptide), 
                            validCleavagePattern)
        return(.Object)
    }
)

#-------------------------------------------------------------------------------








#----FDRs and unique peptides and accessions------------------------------------
setMethod(
    "get_psm_fdr", 
    signature("MSnID"),
    definition=function(.Object)
    {
        isDecoy <- table(factor(.Object@psms$isDecoy, 
                                levels=c(FALSE, TRUE)))
        fdr <- isDecoy["TRUE"]/isDecoy["FALSE"]
        return(fdr)
    }
)




setMethod(
    "get_peptide_fdr", 
    signature("MSnID"),
    definition=function(.Object)
    {
        peptideDecoy <- unique(subset(.Object@psms, 
                                        select=c("peptide","isDecoy")))
        isDecoy <- table(factor(peptideDecoy$isDecoy, levels=c(FALSE, TRUE)))
        fdr <- isDecoy["TRUE"]/isDecoy["FALSE"]
        return(fdr)
    }
)



setMethod(
    "get_accession_fdr", 
    signature("MSnID"),
    definition=function(.Object)
    {
        accessionDecoy <- unique(subset(.Object@psms, 
                                        select=c("accession","isDecoy")))
        isDecoy <- table(factor(accessionDecoy$isDecoy, levels=c(FALSE, TRUE)))
        fdr <- isDecoy["TRUE"]/isDecoy["FALSE"]
        return(fdr)
    }
)




# todo: change to peptides()
setMethod(
    "get_peptides", 
    "MSnID",
    definition=function(.Object)
    {
        unique(as.character(.Object@psms$peptide))
    }
)



# todo: change to accessions()
setMethod(
    "get_accessions", 
    "MSnID",
    definition=function(.Object)
    {
        unique(as.character(.Object@psms$accession))
    }
)
#-------------------------------------------------------------------------------









#----Filter---------------------------------------------------------------------
setMethod(
    "apply_filter", 
    signature(.Object="MSnID", .Filter="character"),
    definition=function(.Object, .Filter)
    {
        # the infamous eval(parse(text=...))
        exprssn <- parse(text=.Filter, srcfile=NULL, keep.source=FALSE)
        #.Object@psms <- subset(.Object@psms, eval(exprssn))
        .Object@psms <- .Object@psms[eval(exprssn),]
        return(.Object)
    }
)

setMethod(
    "apply_filter", 
    signature(.Object="MSnID", .Filter="MSnIDFilter"),
    definition=function(.Object, .Filter)
    {
        filterString <- as(.Filter, "character")
        return(apply_filter(.Object, filterString))
    }
)



# todo: change to returning a vector with names instead of a list
.id_quality <- function(.Object, .Level=c("PSM", "peptide", "accession"))
{
    #
    .Level <- match.arg(.Level)
    if(.Level != "PSM"){
        temp.dt <- .Object@psms[,c(.Level,"isDecoy"),with=FALSE]
        .Object@psms <- unique(temp.dt)
    }
    stopifnot(is.logical(.Object@psms$isDecoy))
    n <- length(.Object@psms$isDecoy)
    decoy <- sum(.Object@psms$isDecoy)
    #
    # catch the case in case there are zero normal matches
    if(decoy == n)
        return(list(fdr=NaN, n=0))
    # if there are normal matches, proceed
    fdr <- decoy/(n-decoy)
    return(c(fdr=fdr, n=n))
}

setMethod(
    "id_quality",
    signature(.Object="MSnID"),
    definition=function(.Object, 
                        filter=NULL, 
                        level=c("PSM", "peptide", "accession"))
    {
        if(!is.null(filter))
            .Object <- apply_filter(.Object, filter)
        # if no filter has been provided just return the quality of
        # features in original object
        level <- match.arg(level)
        return(.id_quality(.Object, level))
    }
)

# Similar of id_quality, but filter hast to be present. It can not be NULL.
setMethod(
    "evaluate_filter",
    signature(.Object="MSnID"),
    definition=function(.Object, 
                        filter, 
                        level=c("PSM", "peptide", "accession"))
    {
        level <- match.arg(level)
        .Object <- apply_filter(.Object, filter)
        return(.id_quality(.Object, level))
    }
)
#-------------------------------------------------------------------------------



setMethod(
    "names",
    signature(x="MSnID"),
    definition=function(x)
    {
        return(colnames(x@psms))
    }
)


setMethod(
    "dim",
    signature(x="MSnID"),
    definition=function(x)
    {
        return(dim(x@psms))
    }
)

#-------------------------------------------------------------------------------





#----Initialization-------------------------------------------------------------
# initialization of MSnID object
MSnID <- function(workDir='.', cleanCache=FALSE)
{
    cachePath <- file.path(workDir,".Rcache")
    # is that enough or should I getOption("R.cache::rootPath")?
    # or is it the same thing? 
    setCacheRootPath(cachePath)
    if(cleanCache) clearCache(cachePath)
    cat("Note, the anticipated/suggested columns in the\n",
        "peptide-to-spectrum matching results are:\n",
        paste(.mustBeColumns, collapse=', '),
        sep='')
    msnidObj <- new("MSnID", workDir=workDir, psms=data.table())
}
#-------------------------------------------------------------------------------

.mustBeColumns <- c("peptide", "accession", "isDecoy", 
                    "calculatedMassToCharge",
                    "experimentalMassToCharge",
                    "chargeState",
                    "spectrumFile", "spectrumID")


setMethod(
    f="psms",
    signature("MSnID"), 
    definition=function(.Object)
    {
        return(as.data.frame(.Object@psms))
    }
)


setReplaceMethod(
    f="psms",
    signature(.Object="MSnID", value="data.frame"),
    definition=function(.Object, value)
    {
        misCol <- setdiff(.mustBeColumns, colnames(value))
        if(!is.null(misCol) & interactive()){
            promptStr <- 
                paste("The data.frame is missing the following columns:\n",
                    paste(strwrap(paste(misCol, collapse=', ')), collapse='\n'),
                    '.\n', 
                    collapse='', sep='')
            warning(promptStr, call. = FALSE, immediate. = TRUE)

            # if use is concerned, do not modify the object
            ANSWER <- readline("Proceed? (Y/N): ")
            if(substr(ANSWER, 1, 1) %in% c("N", "n"))
                return(.Object)
        }

        # if all required columns are present,
        # then there is no need for warnings and questions.
        .Object@psms <- data.table(value)
        return(.Object)
    }
)


setAs(
    from="MSnID", 
    to="data.table",
    def=function(from)
    {
        return(from@psms)
    }
)










#----Misc-----------------------------------------------------------------------
.check_column <- function(.Object, columnName)
# generic column checking
{
    if(is.null(.Object@psms[[columnName]]))
    {
        if(columnName == "peptide"){
            stop("\"peptide\" column is not present!\n",
                    "Note, peptide should be in format: X.XXXXX.X\n",
                    "with flanking amino acids.",
                    call. = FALSE)
        }
    }
    invisible(NULL)
}


setMethod(
    "show", 
    signature("MSnID"),
    definition=function(object)
    {
        cat("MSnID object\n")
        cat("Working directory: \"", object@workDir, "\"\n", sep='')

        # show number of datasets
        try(
            cat("#Spectrum Files: ", 
                length(unique(as.character(object@psms$spectrumFile))), '\n'),
            silent=TRUE)
        
        # show data quality at three levels
        for(i in c("PSM", "peptide", "accession")){
            try({
                temp <- .id_quality(object, i)
                cat("#", i, "s: ", temp['n'], " at ", 
                    signif(100*temp['fdr'], 2), 
                    " % FDR", '\n', sep='')
                }, 
                silent=TRUE)
        }
    }
)


setMethod(
    "read_mzIDs", 
    signature("MSnID"),
    definition=function(object, mzids)
    {
        # check if files are indeed available by the provided path
        stopifnot(all(sapply(mzids, file.exists)))

        # proceed if they are present
        object@psms <- .read_mzIDs.memoized(mzids)
        return(object)
    }
)
#-------------------------------------------------------------------------------




#-------------------------------------------------------------------------------
# Getters and setters modelled after
# library(Biobase)
# getMethod("$", "eSet") # tricky
# getMethod("$<-", "eSet") # straight
# getMethod("[[", "eSet")
# getMethod("[[<-", "eSet")

# Note, get method in Biobase is a bit different
# library(Biobase)
# getMethod("$", "eSet")
setMethod("$", "MSnID",
            definition=function(x, name)
            {
                return(x@psms[[name]])
            }
)

setMethod("[[", "MSnID",
            definition=function(x, i, j, ...)
            {
                return(x@psms[[i]])
            }
)



# todo: double check set $<- method
# todo: the current problem is that := assignment by reference
# still does not work properly (copies by value). Therefore
# will use either set or find smarter way to use :=

# # data.frame style
# setMethod("$<-", "MSnID",
#           definition=function(x, name, value)
#           {
#              x@psms[[name]] <- value
#              return(x)
#           }
# )
setMethod("$<-", "MSnID",
    definition=function(x, name, value)
    {
        x@psms <- eval(substitute(x@psms[,NAME_ARG:=value],
                                    list(NAME_ARG=name)))
        return(x)
    }
)


# setMethod("$<-", "MSnID",
#     definition=function(x, name, value)
#     {
#         x@psms <- eval(substitute(x@psms[,NAME_ARG:=value],
#                                   list(NAME_ARG=name)))
#         return(x)
#     }
# )



setMethod("[[<-", "MSnID",
            definition=function(x, i, j, ..., value)
            {
                x@psms[[i, ...]] <- value
                return(x)
            }
)
#-------------------------------------------------------------------------------







#--------------------------------------------
setMethod("correct_peak_selection", "MSnID",
    definition=function(.Object)
    {
        deltaMz <- .Object$experimentalMassToCharge - 
                    .Object$calculatedMassToCharge
        deltaMass <- deltaMz * .Object$chargeState
        deltaC13C12 <- 1.0033548378
        deltaMass.corrected <- deltaMass - deltaC13C12 * round(deltaMass)
        deltaMz.corrected <- deltaMass.corrected/.Object$chargeState
        .Object$experimentalMassToCharge <- .Object$calculatedMassToCharge +
                                                deltaMz.corrected
        return(.Object)
    }
)

setMethod("mass_measurement_error", "MSnID",
    definition=function(.Object)
    {
        ppm <- 1e6*(.Object$experimentalMassToCharge - 
                    .Object$calculatedMassToCharge)/
                .Object$calculatedMassToCharge
        return(ppm)
    }
)

setMethod("recalibrate", "MSnID",
    definition=function(.Object)
    {
        # this is just a simple shift by a signle bias
        error.ppm <- mass_measurement_error(.Object)
        #%%%%%%%%%%%%%%%%%%%%%%%
        # modeling systematic error
        sys.error.ppm <- median(error.ppm)
        # how about expectation maximization
        # ...
        # position of the point with most density (that is mode)
        bw <- 1 # 1 ppm bin width
        # tie n with the range. 3 is 2^3=8 is extra precision factor.
        # that is points will go ~ 0.125 of ppm
        n <- 2^ceiling(3 + log2(diff(range(error.ppm))/bw))
        n <- max(512, n)
        d <- density(error.ppm, bw=bw, n=n)
        sys.error.ppm <- d$x[which.max(d$y)]
        #%%%%%%%%%%%%%%%%%%%%%%%
        res.error.ppm <- error.ppm - sys.error.ppm # new error residuals
        # now back-calculate the experimentalMassToCharge
        .Object$experimentalMassToCharge <- 
            .Object$calculatedMassToCharge * (1 + res.error.ppm/1e6)
        return(.Object)
    }
)





.convert_MSnID_to_MSnSet <- function(msnid)
{
    #--- exprs data. peptide level
    # applying unique to remove redundant peptide records
    # because of peptide-protein mapping redundancy
    quant <- unique(psms(msnid)[, c("peptide","spectrumFile","spectrumID")])
    exprs.data <- acast(quant, peptide ~ spectrumFile,
                        value.var="spectrumID",
                        fun.aggregate=length)

    #--- feature data
    # todo. switch from subset to faster data.table operation.
    # may be "with=FALSE"
    feature.data <- unique(subset(msnid@psms,
                            select=c("peptide","accession")))
    mapping <- with(feature.data, tapply(accession, peptide, list))
    feature.data <- data.frame(peptide=names(mapping))
    feature.data$accession <- mapping
    rownames(feature.data) <- feature.data$peptide
    feature.data <- feature.data[rownames(exprs.data),]
    feature.data <- new("AnnotatedDataFrame", data=feature.data)

    #--- make MSnSet
    msnset <- new("MSnSet", exprs=exprs.data, featureData=feature.data)
    return(msnset)
}



setAs("MSnID", "MSnSet",
        def=function(from) .convert_MSnID_to_MSnSet(from))



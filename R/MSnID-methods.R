

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
    definition=function(object, missedCleavagePattern)
    {
        .check_column(object, "peptide")
        object@psms$numMissCleavages <-
            .assess_missed_cleavages(as.character(object@psms$peptide),
                                    missedCleavagePattern)
        return(object)
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
    definition=function(object, validCleavagePattern)
    {
        object@psms$numIrregCleavages <-
            .assess_termini(as.character(object@psms$peptide),
                            validCleavagePattern)
        return(object)
    }
)

#-------------------------------------------------------------------------------





#--- Accessors -----------------------------------------------------------------
setMethod(
    "peptides",
    "MSnID",
    definition=function(object)
    {
        unique(as.character(object@psms$peptide))
    }
)



setMethod(
    "accessions",
    "MSnID",
    definition=function(object, ...)
    {
        unique(as.character(object@psms$accession))
    }
)

setMethod(
    "proteins",
    "MSnID",
    definition=function(object, ...)
    {
        accessions(object)
    }
)
#-------------------------------------------------------------------------------







#----Filter---------------------------------------------------------------------
# # Old implementation. Just as a backup.
# # See https://github.com/vladpetyuk/MSnID/issues/5 for the issue.
# setMethod(
#     "apply_filter",
#     signature(msnidObj="MSnID", filterObj="character"),
#     definition=function(msnidObj, filterObj)
#     {
#         exprssn <- parse(text=filterObj, srcfile=NULL, keep.source=FALSE)
#         msnidObj@psms <- msnidObj@psms[eval(exprssn),]
#         return(msnidObj)
#     }
# )
setMethod(
    "apply_filter",
    signature(msnidObj="MSnID", filterObj="character"),
    definition=function(msnidObj, filterObj)
    {
        exprssn <- parse(text=filterObj, srcfile=NULL, keep.source=FALSE)
        idx <- eval(exprssn, envir = msnidObj@psms, enclos = parent.frame())
        msnidObj@psms <- msnidObj@psms[idx,]
        return(msnidObj)
    }
)

setMethod(
    "apply_filter",
    signature(msnidObj="MSnID", filterObj="MSnIDFilter"),
    definition=function(msnidObj, filterObj)
    {
        filterString <- as(filterObj, "character")
        return(apply_filter(msnidObj, filterString))
    }
)

.id_quality <- function(object, .Level=c("PSM", "peptide", "accession"))
{
    #
    .Level <- match.arg(.Level)
    if(.Level != "PSM"){
        temp.dt <- object@psms[,c(.Level,"isDecoy"),with=FALSE]
        object@psms <- unique(temp.dt)
    }else{
        # deal with peptide to protein assignment redundancy
        temp.dt <- object@psms[,c('spectrumFile',
                                  'spectrumID',
                                  'peptide',
                                  'isDecoy'),with=FALSE]
        object@psms <- unique(temp.dt)
    }
    
    stopifnot(is.logical(object@psms$isDecoy))
    n <- length(object@psms$isDecoy)
    decoy <- sum(object@psms$isDecoy)
    #
    # catch the case in case there are zero normal matches
    if(decoy == n)
        return(c(fdr=NaN, n=0))
    # if there are normal matches, proceed
    fdr <- decoy/(n-decoy)
    return(c(fdr=fdr, n=n))
}

setMethod(
    "id_quality",
    signature(object="MSnID"),
    definition=function(object,
                        filter=NULL,
                        level=c("PSM", "peptide", "accession"))
    {
        if(!is.null(filter))
            object <- apply_filter(object, filter)
        # if no filter has been provided just return the quality of
        # features in original object
        level <- match.arg(level, several.ok = TRUE)
        out <- t(sapply(level, .id_quality, object=object))
        return(out)
    }
)

# Similar of id_quality, but filter hast to be present. It can not be NULL.
setMethod(
    "evaluate_filter",
    signature(object="MSnID"),
    definition=function(object,
                        filter,
                        level=c("PSM", "peptide", "accession"))
    {
        level <- match.arg(level, several.ok = TRUE)
        object <- apply_filter(object, filter)
        out <- t(sapply(level, .id_quality, object=object))
        return(out)
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
    msg <- paste("Note, the anticipated/suggested columns in the",
                 "peptide-to-spectrum matching results are:",
                 "-----------------------------------------------",
                 paste0(sort(.mustBeColumns), collapse = "\n"),
                 sep='\n')
    message(msg)
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
    definition=function(object, ...)
    {
        return(as.data.frame(object@psms))
    }
)


setReplaceMethod(
    f="psms",
    signature(object="MSnID", value="data.frame"),
    definition=function(object, value)
    {
        misCol <- setdiff(.mustBeColumns, colnames(value))
        if((length(misCol) > 0) & interactive()){
            promptStr <-
                paste("The data.frame is missing the following columns:\n",
                    paste(strwrap(paste(misCol, collapse=', ')), collapse='\n'),
                    '.\n',
                    collapse='', sep='')
            warning(promptStr, call. = FALSE, immediate. = TRUE)

            # if use is concerned, do not modify the object
            ANSWER <- readline("Proceed? (Y/N): ")
            if(substr(ANSWER, 1, 1) %in% c("N", "n"))
                return(object)
        }

        # if all required columns are present,
        # then there is no need for warnings and questions.
        object@psms <- data.table(value)
        return(object)
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
.check_column <- function(object, columnName)
# generic column checking
{
    if(is.null(object@psms[[columnName]]))
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
    definition=function(object, mzids, backend)
    {
        # check if files are indeed available by the provided path
        stopifnot(all(sapply(mzids, file.exists)))
        backend <- match.arg(backend)
        stopifnot(backend %in% c('mzID','mzR')) # not sure if this is necessary
        
        # Try to load cached data, if exists
        key <- list(mzids)
        data <- loadCache(key)
        if (!is.null(data)){
            message("Loaded cached data\n")
        }else{
            message("Reading from mzIdentMLs ...\n")
            if(backend == 'mzID'){
                data <- data.table(flatten(mzID(mzids), safeNames=FALSE))
                data$databaseFile <- basename(gsub('\\\\','/',data$databaseFile))
            }else{
                data <- .read_mzIDs.mzR(mzids)
            }
            
            #' extra stuff
            data$peptide <- paste(data$pre, data$pepSeq, data$post, sep='.')
            data$spectrumFile <- basename(gsub('\\\\','/',data$spectrumFile))
            
            #' save for reuse
            saveCache(data, key=key)
        }
        
        object@psms <- data
        
        return(object)
    }
)
#-------------------------------------------------------------------------------




#-------------------------------------------------------------------------------
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

setMethod("$<-", "MSnID",
    definition=function(x, name, value)
    {
        x@psms[[name]] <- value
        return(x)
    }
)

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
    definition=function(object)
    {
        deltaMz <- object$experimentalMassToCharge -
                    object$calculatedMassToCharge
        deltaMass <- deltaMz * object$chargeState
        deltaC13C12 <- 1.0033548378
        deltaMass.corrected <- deltaMass - deltaC13C12 * round(deltaMass)
        deltaMz.corrected <- deltaMass.corrected/object$chargeState
        object$experimentalMassToCharge <- object$calculatedMassToCharge +
                                                deltaMz.corrected
        return(object)
    }
)

setMethod("mass_measurement_error", "MSnID",
    definition=function(object)
    {
        ppm <- 1e6*(object$experimentalMassToCharge -
                    object$calculatedMassToCharge)/
                object$calculatedMassToCharge
        return(ppm)
    }
)

setMethod("recalibrate", "MSnID",
    definition=function(object)
    {
        # this is just a simple shift by a signle bias
        error.ppm <- mass_measurement_error(object)
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
        object$experimentalMassToCharge <-
            object$calculatedMassToCharge * (1 + res.error.ppm/1e6)
        return(object)
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



utils::globalVariables(c("accession", "N", "pepSeq", "num"))

setMethod("infer_parsimonious_accessions", "MSnID",
          definition=function(object, unique_only=FALSE, prior=character(0))
          {
              infer_acc <- function(x){
                  res <- list()
                  while(nrow(x) > 0){
                      top_prot <- x[, .N, by=accession][which.max(N),,]$accession
                      top_peps <- subset(x, accession == top_prot)
                      # top_peps <- x[accession == top_prot, c(1,2), with=FALSE] # slower
                      res <- c(res, list(top_peps))
                      x <- subset(x, !(pepSeq %in% top_peps[[1]]))
                      # x <- x[!top_peps, on=.(pepSeq)] # slower
                  }
                  return(rbindlist(res, use.names=F, fill=FALSE, idcol=NULL))
              }
              
              x <- unique(object@psms[,.(pepSeq, accession)])
              setorder(x, accession, pepSeq)
              
              if(unique_only){
                  redundancy <- x[, .(num = length(accession)), by = pepSeq]
                  non_redundant <- redundancy[num == 1]
                  res <- non_redundant[,.(pepSeq)]
                  setkey(res, pepSeq)
              }else{ # consider razor peptides as well
                  
                  # peptides from proteins justified by prior
                  xp <- x[accession %in% prior][,.(pepSeq)]
                  xp <- unique(xp)
                  #
                  x_prior <- x[xp, on=.(pepSeq)] # semi_join
                  x_current <- x[!xp, on=.(pepSeq)] # anti_join
                  
                  # this removes other proteins mapped to 
                  # peptides from prior justified proteins
                  res_prior <- x_prior[accession %in% prior]
                  # key step. the slowest
                  res_current <- infer_acc(x_current)
                  #
                  res <- rbindlist(list(res_prior, res_current))
                  setkey(res, pepSeq, accession)
              }
              
              old_psms <- psms(object)
              setDT(old_psms)
              setkey(old_psms, pepSeq, accession)
              new_psms <- old_psms[res] # semi_join
              psms(object) <- new_psms
              
              return(object)
          }
)



setMethod("report_mods", 
          signature = signature("MSnID"),
          definition = function(object)
          {
              mods <- psms(object)$modification
              mod_masses <- lapply(mods, strsplit, "\\s\\(\\d+\\),?\\s?")
              return(table(unlist(mod_masses)))
          }
)


setMethod("add_mod_symbol", "MSnID",
          definition=function(object, mod_mass, symbol)
          {
              .add_mod_symbol(object, mod_mass, symbol)
          }
)


setMethod("map_mod_sites", "MSnID",
          definition=function(object, 
                              fasta, 
                              accession_col = "accession", 
                              peptide_mod_col = "peptide_mod", 
                              mod_char = "*", 
                              site_delimiter = "lower")
          {
              ids <- psms(object)
              
              # test if decoys are in fasta
              # if not, then add reversed sequences
              # for now we'll support only reverse as decoys
              decoy_acc <- apply_filter(object, "isDecoy")[[accession_col]]
              if(length(decoy_acc) > 0 & !any(decoy_acc %in% names(fasta))){
                  fasta_rev <- reverse(fasta)
                  names(fasta_rev) <- paste0("XXX_",names(fasta))
                  fasta <- c(fasta, fasta_rev)
              }
              
              res <- .map_mod_sites(ids, 
                                    fasta, 
                                    accession_col, 
                                    peptide_mod_col, 
                                    mod_char)
              
              # make site ID from accession_col and 
              # SiteCollapsedFirst
              if(site_delimiter == "lower"){
                  res <- mutate(res, 
                                SiteID = paste0(!!sym(accession_col), 
                                                "-", 
                                                gsub("([[:upper:]])(\\d+),?",
                                                     "\\1\\2\\L\\1", 
                                                     SiteCollapsedFirst, 
                                                     perl = T)))
              }else{
                  res <- mutate(res, 
                                SiteID = paste0(!!sym(accession_col), 
                                                "-", 
                                                gsub(",",
                                                     site_delimiter,
                                                     SiteCollapsedFirst)))
              }
              psms(object) <- res
              return(object)
          }
)






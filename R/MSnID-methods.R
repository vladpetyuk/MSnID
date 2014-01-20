

#----Peptide Sequence Handling--------------------------------------------------

.get_clean_peptide_sequence <- function(peptide){
   return(.strip_flanking_AAs(
      .strip_modifications_from_peptide(peptide)))
}

.strip_modifications_from_peptide <- function(peptide)
{
   # INPUT: Peptide sequence with flanking AAs e.g. "K.APEP*TID{34.34}E%.-"
   # OUTPUT: sequence with just AA symbols e.g. "K.APEPTIDE.-"
   # clean from the mods
   # make sure it removes "." in the middle of the sequence
   aminoAcids <- "ARNDCEQGHILKMFPSTWYV"
   extendedAAs <- "BUZXO" # http://en.wikipedia.org/wiki/Proteinogenic_amino_acid
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
                                    missedCleavagePattern="[KR](?=[^P$])")
{
   # the fastest approach
   peptide <- .get_clean_peptide_sequence(peptide)
   # now no flanking AAs, no mods
   numberOfMissedCleavages <- 
      nchar(peptide) - 
      nchar(gsub(missedCleavagePattern, "", peptide, perl=TRUE))
   return( numberOfMissedCleavages )
}


setMethod("assess_missed_cleavages", signature("MSnID"),
          definition=function(.Object, ...)
          {
             .check_column(.Object, "Peptide")
             .Object@psms$NumMissCleavages <-
                .assess_missed_cleavages(as.character(.Object@psms$Peptide), ...)
             return(.Object)
          }
)

#---

.assess_termini <- function(peptide, 
                           validCleavagePattern="[RK]\\.[^P]")
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

             

setMethod("assess_termini", signature("MSnID"),
          definition=function(.Object, ...)
          {
             .Object@psms$NumIrregCleavages <- 
                .assess_termini(as.character(.Object@psms$Peptide), ...)
             return(.Object)
          }
)

#-------------------------------------------------------------------------------








#----FDRs and unique peptides and accessions------------------------------------
setMethod("get_psm_fdr", "MSnID",
          definition=function(.Object)
          {
             isDecoy <- table(factor(.Object@psms$isDecoy, levels=c(FALSE, TRUE)))
             fdr <- isDecoy["TRUE"]/isDecoy["FALSE"]
             return(fdr)
          }
)




setMethod("get_peptide_fdr", "MSnID",
          definition=function(.Object)
          {
             peptideDecoy <- unique(subset(.Object@psms, select=c("Peptide","isDecoy")))
             isDecoy <- table(factor(peptideDecoy$isDecoy, levels=c(FALSE, TRUE)))
             fdr <- isDecoy["TRUE"]/isDecoy["FALSE"]
             return(fdr)
          }
)



setMethod("get_accession_fdr", "MSnID",
          definition=function(.Object)
          {
             accessionDecoy <- unique(subset(.Object@psms, select=c("Accession","isDecoy")))
             isDecoy <- table(factor(accessionDecoy$isDecoy, levels=c(FALSE, TRUE)))
             fdr <- isDecoy["TRUE"]/isDecoy["FALSE"]
             return(fdr)
          }
)





setMethod("get_peptides", "MSnID",
          definition=function(.Object)
          {
             unique(as.character(.Object@psms$Peptide))
          }
)




setMethod("get_accessions", "MSnID",
          definition=function(.Object)
          {
             unique(as.character(.Object@psms$Accession))
          }
)
#-------------------------------------------------------------------------------









#----Filter---------------------------------------------------------------------
setMethod("apply_filter", 
          signature(.Object="MSnID", .Filter="character"),
          definition=function(.Object, .Filter)
          {
             # the infamous eval(parse(text=...))
             .Object@psms <- subset(.Object@psms, eval(parse(text=.Filter)))
             return(.Object)
          }
)

setMethod("apply_filter", 
          signature(.Object="MSnID", .Filter="MSnIDFilter"),
          definition=function(.Object, .Filter)
          {
             filterString <- as(.Filter, "character")
             return(apply_filter(.Object, filterString))
          }
)





# How to do same dispatch of MSnIDFilter and character?
.id_quality <- function(.Object, .Level=c("PSM", "Peptide", "Accession"))
{
   #
   .Level <- match.arg(.Level)
   if(.Level != "PSM")
      .Object@psms <- unique(subset(.Object@psms, select=c(.Level,"isDecoy")))
   isDecoyTbl <- table(factor(.Object@psms$isDecoy, levels=c(FALSE,TRUE)))
   #
   # catch the case in case there are zero normal matches
   if(isDecoyTbl["FALSE"] == 0)
      return(list(fdr=NaN, n=0))
   # if there are normal matches, proceed
   fdr <- isDecoyTbl["TRUE"]/isDecoyTbl["FALSE"]
   return(list(fdr=as.numeric(fdr), 
               n=sum(isDecoyTbl)))
}

setMethod("id_quality",
          signature(.Object="MSnID"),
          definition=function(.Object, filter=NULL, level=c("PSM", "Peptide", "Accession"))
          {
             if(!is.null(filter))
                .Object <- apply_filter(.Object, filter)
             level <- match.arg(level)
             return(.id_quality(.Object, level))
          }
)

# synonym of id_quality, but with filter present. It can not be NULL.
setMethod("evaluate_filter",
          signature(.Object="MSnID"),
          definition=function(.Object, filter, level=c("PSM", "Peptide", "Accession"))
          {
             level <- match.arg(level)
             .Object <- apply_filter(.Object, filter)
             return(.id_quality(.Object, level))
          }
)
#-------------------------------------------------------------------------------





#----Getters and Setters--------------------------------------------------------
setMethod("set_psm_parameter", 
          signature(.Object="MSnID"),
          definition=function(.Object, ...)
          {
             e <- eval(substitute(list(...)), .Object@psms, parent.frame())
             stopifnot(length(e) == 1)
             stopifnot(length(newPar) == nrow(.Object@psms))
             .Object@psms[[names(e)]] <- e[[1]]
             return(.Object)
          }
)



setMethod("get_psm_parameter", 
          signature(.Object="MSnID", parName="character"),
          definition=function(.Object, parName)
          {
             return(.Object@psms[[parName]])
          }
)



setMethod("delete_psm_parameter", 
          signature(.Object="MSnID"),
          definition=function(.Object, ...)
          {
             parNames <- sapply(as.list(substitute(list(...)))[-1L], as.character)
             for( i in parNames) {.Object@psms[[i]] <- NULL}
             return(.Object)
          }
)



setMethod("names",
          signature(x="MSnID"),
          definition=function(x)
          {
             return(colnames(x@psms))
          }
)


setMethod("dim",
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
   msnidObj <- new("MSnID", workDir=workDir)
}
#-------------------------------------------------------------------------------

.mustBeColumns <- c("Peptide", "Accession", "isDecoy", 
                    "calculatedMassToCharge",
                    "experimentalMassToCharge", 
                    "SpectrumFile")

setGeneric("psms",
           function(.Object, ...) standardGeneric("psms"))
setMethod("psms","MSnID", 
          function(.Object) .Object@psms)

setGeneric("psms<-",
           function(.Object, value) standardGeneric("psms<-"))
setReplaceMethod("psms",
                 signature(.Object="MSnID",
                           value="data.frame"),
                 function(.Object, value) {
                    misCol <- setdiff(.mustBeColumns, colnames(value))
                    if(!is.null(misCol) & interactive()){
                       promptStr <- paste("The data.frame missing the following columns:\n",
                                          paste(strwrap(paste(misCol, collapse=', ')), 
                                                collapse='\n'),
                                          '.\n',
                                          collapse='', sep='')
                       warning(promptStr, call. = FALSE, immediate. = TRUE)
                       ANSWER <- readline("Proceed? (Y/N): ")
                       if (substr(ANSWER, 1, 1) == "n")
                          return(.Object)
                    }
                    .Object@psms <- value
                    return(.Object)
                 })


setAs("MSnID", "data.frame",
      def=function(from) return(.Object@psms))










   
#----Misc-----------------------------------------------------------------------
.check_column <- function(.Object, columnName)
   # generic column checking
{
   if(is.null(.Object@psms[[columnName]])){
      if(columnName == "Peptide"){
         stop("Peptide column is not present!\n",
              "Note, Peptide should be in format: X.XXXXX.X\n",
              "with flanking amino acids.",
              call. = FALSE)
      }
   }
   invisible(NULL)
}


setMethod("show", "MSnID",
          definition=function(object)
          {
             cat("MSnID object\n")
             cat("Working directory: \"", object@workDir, "\"\n", sep='')
             try(
                cat("#Spectrum Files: ", 
                    length(unique(as.character(object@psms$SpectrumFile))), '\n'),
                silent=TRUE)
             #
             for(i in c("PSM", "Peptide", "Accession")){
                try({
                   temp <- .id_quality(object, i)
                   cat("#", i, "s: ", temp$n, " at ", signif(100*temp$fdr, 2), 
                       " % FDR", '\n', sep='')
                }, silent=TRUE)
             }
          }
)


setMethod("read_mzIDs", "MSnID",
          definition=function(.Object, pathToMzIDs)
          {
             # to check if files are indeed available by the provided path
             stopifnot(all(sapply(pathToMzIDs, file.exists)))
             .Object@psms <- .read_mzIDs.memoized(pathToMzIDs)
             return(.Object)
          }
)
#-------------------------------------------------------------------------------






setMethod("$", "MSnID",
          definition=function(x, name)
          {
             eval(substitute(x@psms$NAME_ARG, list(NAME_ARG = name)))
             # return(x@psms[[name]])
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

# getGeneric("[[<-")
# getMethod("[[<-", "eSet")
setMethod("[[<-", "MSnID",
          definition=function(x, i, j, ..., value)
          {
             x@psms[[i, ...]] <- value
             return(x)
          }
)



# library(Biobase)
# getMethod("$", "eSet")
# getMethod("$<-", "eSet")
# getMethod("[[", "eSet")
# getMethod("[[<-", "eSet")


#--------------------------------------------
setGeneric("correct_peak_selection", 
           function(.Object, ...) standardGeneric("correct_peak_selection"))
           
setMethod("correct_peak_selection", "MSnID",
          definition=function(.Object, ...)
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

setGeneric("mass_measurement_error",
           function(.Object, ...) standardGeneric("mass_measurement_error"))

setMethod("mass_measurement_error", "MSnID",
          definition=function(.Object, ...)
          {
             ppm <- 1e6*(.Object$experimentalMassToCharge - 
                            .Object$calculatedMassToCharge)/
                .Object$calculatedMassToCharge
             return(ppm)
          }
)

setGeneric("recalibrate",
           function(.Object, ...) standardGeneric("recalibrate"))

setMethod("recalibrate", "MSnID",
          definition=function(.Object, ...)
          {
             # this is just a simple shift by a signle bias
             error.ppm <- mass_measurement_error(.Object)
             #%%%%%%%%%%%%%%%%%%%%%%%
             # modeling systematic error
             sys.error.ppm <- median(error.ppm)
             # how about expectation maximization
             # ...
             # most dense position (akin mode)
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
             .Object$experimentalMassToCharge <- .Object$calculatedMassToCharge *
                (1 + res.error.ppm/1e6)
             return(.Object)
          }
)
             




.convert_MSnID_to_MSnSet <- function(msnid)
{
   library("reshape2") # this may be taken care of by package
   #--- exprs data. peptide level
   # call to unique to remove inflation due to peptide-protein mapping redundancy
   quant <- unique(msnid@psms[, c("Peptide","SpectrumFile","spectrumID")])
   exprs.data <- acast(quant, Peptide ~ SpectrumFile, 
                       value.var="spectrumID", 
                       fun.aggregate=length)
   
   #--- feature data
   feature.data <- unique(subset(msnid@psms,
                                 select=c("Peptide","Accession")))
   mapping <- with(feature.data, tapply(Accession, Peptide, list))
   feature.data <- data.frame(Peptide=names(mapping))
   feature.data$Accession <- mapping
   rownames(feature.data) <- feature.data$Peptide
   feature.data <- feature.data[rownames(exprs.data),]
   feature.data <- new("AnnotatedDataFrame", data=feature.data)
   
   #--- make MSnSet
   msnset <- new("MSnSet", exprs=exprs.data, featureData=feature.data)
   return(msnset)
}

setAs("MSnID", "MSnSet",
      def=function(from) .convert_MSnID_to_MSnSet(from))



combineFeatures <- function(object, groupBy, redundancy.handler=c("ignore","unique.only"), ...)
{
   # wrapper to combineFeatures to handle redundancy in feature to factor mapping
   # e.g. peptide-to-protein redundancy
   if(!is.list(groupBy)){
      result <- MSnbase::combineFeatures(object, groupBy, ...)
   }else{
      # handling of the redundancy
      if(any(names(groupBy) != rownames(fData(object))))
         stop("names of groupBy list do not match fData of the MSnSet object")
      redundancy.handler <- match.arg(redundancy.handler)
      if(redundancy.handler == "ignore"){
         
         expansion.index <- rep(seq_len(nrow(object)), sapply(groupBy, length))
         new.exprs <- exprs(object)[expansion.index,]
         rownames(new.exprs) <- NULL
         groupBy.idx <- sapply(fData(object), identical, groupBy)
         new.feature.data <- fData(object)[expansion.index,]
         new.feature.data[,groupBy.idx] <- unlist(groupBy)
         rownames(new.feature.data) <- NULL
         new.object <- new("MSnSet", exprs = new.exprs, 
                           featureData = new("AnnotatedDataFrame", 
                                             data = new.feature.data),
                           phenoData=phenoData(object))
         result <- MSnbase::combineFeatures(new.object, unlist(groupBy), ...)
      }else if(redundancy.handler == "unique.only"){
         idx.unique <- sapply(groupBy, length) < 2
         object <- object[idx.unique,]
         groupBy <- unlist(groupBy[idx.unique])
         result <- MSnbase::combineFeatures(object, groupBy, ...)
      }else{
         stop("Method \"", redundancy.handler, 
              "\" for handing the redundancy is not implemented!", sep='')
      }
      return(result)
   }
}

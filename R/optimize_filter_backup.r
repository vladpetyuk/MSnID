# # 
# library("plyr")
# # 
# # loaded.names <- load("psms.RData")
# # msmsdata <- psms
# # rm(loaded.names)
# # msmsdata <- within(msmsdata, {Protein <- accession
# #                               accession <- NULL
# #                               Peptide <- peptide
# #                               peptide <- NULL
# #                               negLog10_SpecProb <- msmsScore
# #                               msmsScore <- NULL
# #                               Dataset <- spectrumFile
# #                               spectrumFile <- NULL})
# # msmsdata <- subset(msmsdata, abs(parentMassErrorPPM) < 10)
# # # data.frame with
# # # Dataset
# # # -- Feature characteristics --
# # # Peptide, Protein, Scan_Number, Charge_State, {perhaps modifications}
# # # -- MSMS identification confidence --
# # # parentMassErrorPPM, negLog10_SpecProb
# # 
# # 
# # # Now augment data with extra stuff
# # #---------------------------------------------------------
# # msmsdata <- transform(msmsdata, 
# #                isContaminant=is_contaminant(Protein, 
# #                                          contaminantPattern="^Contaminant_"),
# #                isDecoy=is_decoy(Protein, decoyPattern="^XXX_"),
# #                NumMissCleav=number_of_missed_cleavages(Peptide),
# #                NumIrregCleav=number_of_irregular_cleavage_termini(Peptide))
# # #--------------------------------------------------------
# # 
# # 
# # 
# # 
# # psm_fdr <- function(msmsdata)
# # {
# #    isDecoy <- table(factor(msmsdata$isDecoy, levels=c(FALSE, TRUE)))
# #    fdr <- isDecoy["TRUE"]/isDecoy["FALSE"]
# #    return(fdr)
# # }
# # 
# # 
# # 
# # peptide_fdr <- function(msmsdata)
# # {
# #    msmsdata <- unique(subset(msmsdata, select=c("Peptide","isDecoy")))
# #    isDecoy <- table(factor(msmsdata$isDecoy, levels=c(FALSE, TRUE)))
# #    fdr <- isDecoy["TRUE"]/isDecoy["FALSE"]
# #    return(fdr)
# # }
# # 
# # 
# # 
# # 
# eval_MSMS_ID_filter <- function(msmsdata, 
#                                 filterString, 
#                                 level=c("Peptide","PSM","Protein"))
# {
#    level <- match.arg(level)
#    if(level == "PSM")
#       msmsdata$PSM <- seq_len(nrow(msmsdata))
#    msmsdata <- subset(msmsdata, eval(parse(text=filterString)))
#    msmsdata <- unique(subset(msmsdata, select=c(level,"isDecoy")))
#    isDecoyTbl <- table(factor(msmsdata$isDecoy, levels=c(FALSE,TRUE)))
#    # catch the case in case there are zero normal matches
#    if(isDecoyTbl["FALSE"] == 0)
#       return(list(fdr=NaN, num.pep=0))
#    # if there are normal matches, proceed
#    fdr <- isDecoyTbl["TRUE"]/isDecoyTbl["FALSE"]
#    return(list(fdr=as.numeric(fdr), 
#                num.pep=as.numeric(isDecoyTbl["FALSE"])))
# }
# # 
# # 
# # 
# # 
# # #----------------------------------------
# # # get the number of peptides and FDR
# get_num_pep_for_fdr <- function(pars, msmsdata, FDR_MAX, ...) 
# {
#    progress_bar$step()
#    filterString <- do.call("sprintf",
#                         c("negLog10_SpecProb > %s & abs(parentMassErrorPPM) < %s", 
#                              as.list(pars)))
#    x <- eval_MSMS_ID_filter(msmsdata, filterString, ...)
#    if(is.na(x$fdr) || x$fdr > FDR_MAX){
#       return(rnorm(1,sd=0.001)) # 0 is bad because optimization does not move
#    }else{
#       return(x$num.pep)
#    }
# }
# # 
# # 
# # # Optimizaton Block
# # #---------------------------------------------------------------
# # #(pars)
# # # MSGF and ppm
# # get_good_starting_parameters <- function(msmsdata, FDR_MAX)
# #    # procedure to find approximate cut-offs for given FDR
# #    # Not sure if tihs is PSM or unique peptide level, but
# #    # for now it is fine anyway.
# # {
# #    parStart <- list(negLog10_SpecProb=NA,
# #                     parentMassErrorPPM=NA)
# #    #--
# # #    msmsdata <- msmsdata[order(-msmsdata$negLog10_SpecProb),]
# # #    fdrs <- cumsum(msmsdata$isDecoy)/seq_along(msmsdata$isDecoy)
# # #    parStart$negLog10_SpecProb <- msmsdata[max(which(fdrs < FDR_MAX)),
# # #                                    "negLog10_SpecProb"]
# # #    #--
# # #    msmsdata <- msmsdata[order(abs(msmsdata$parentMassErrorPPM)),]
# # #    fdrs <- cumsum(msmsdata$isDecoy)/seq_along(msmsdata$isDecoy)
# # #    parStart$parentMassErrorPPM <- abs(msmsdata[max(which(fdrs < FDR_MAX)),
# # #                                     "parentMassErrorPPM"])
# #    #--
# #    # or just plain means
# #    parStart$parentMassErrorPPM <- mean(abs(msmsdata$parentMassErrorPPM))
# #    parStart$negLog10_SpecProb <- mean(msmsdata$negLog10_SpecProb)
# #    
# #    return(parStart)
# # }
# # 
# # 
# # 
# msms_id_optimizer <- function( msmsdata, FDR_MAX, level)
# {
#    #
#    # parStart <- get_good_starting_parameters(msmsdata, FDR_MAX)
#    parStart <- list(parentMassErrorPPM=10, negLog10_SpecProb=8) # dummy values
#    #
#    progress_bar <<- create_progress_bar("text")
#    progress_bar$init(100)
#    optim.out.nm <- optim(par=parStart, 
#                       fn = get_num_pep_for_fdr, 
#                       msmsdata = msmsdata,
#                       FDR_MAX = FDR_MAX, 
#                       level=level,
#                       method="Nelder-Mead",
#                       control=list(fnscale=-1, maxit=100))
#    optim.out.nm$method <- "Nelder-Mead"
#    progress_bar$term()
#    #
#    progress_bar <<- create_progress_bar("text")
#    progress_bar$init(100)
#    optim.out.sann <- optim(par=parStart, 
#                       fn = get_num_pep_for_fdr, 
#                       msmsdata = msmsdata,
#                       FDR_MAX = FDR_MAX, 
#                       level=level,
#                       method="SANN",
#                       control=list(fnscale=-1, maxit=100))
#    optim.out.sann$method <- "SANN"
#    progress_bar$term()
#    #
#    # select the best result
#    optim.out.best <- switch(
#                            which.max(c(optim.out.nm$value, 
#                                        optim.out.sann$value)),
#                            optim.out.nm,
#                            optim.out.sann)
#    return(optim.out.best)
# }
# 
# # 
# # 
# # msmsdata <- subset(msmsdata, NumIrregCleav == 0)
# # #
# # 
# # system.time(
# #    optim.out <- msms_id_optimizer(msmsdata, FDR_MAX=0.01, level="Peptide"))
# ##---------------------------------------------------------------------------
# msmsdata <- msnid@psms
# msmsdata <- within(msmsdata, {Protein <- accession
#                               accession <- NULL
#                               Peptide <- peptide
#                               peptide <- NULL
#                               negLog10_SpecProb <- msmsScore
#                               msmsScore <- NULL
#                               Dataset <- spectrumFile
#                               spectrumFile <- NULL})
# # What I need to report?
# # A *character string* along the some basic summary statistics.
# # What statistics?
# # Number of items: psms/peptides/accessions
# # FDR
# # optimization output as list
# system.time(optim.out <- msms_id_optimizer(msmsdata, FDR_MAX=0.01, level="Peptide"))
# 
# 
# 
# 
# 
# 
# 
# # 
# # filterString <- do.call("sprintf",
# #                         c("negLog10_SpecProb > %s & abs(parentMassErrorPPM) < %s", 
# #                           as.list(optim.out$par)))
# # msmsdata.filtered <- subset(msmsdata, eval(parse(text=filterString)))
# # 
# # # filterString = "negLog10_SpecProb > 9 & abs(parentMassErrorPPM) < 10"
# # # eval_MSMS_ID_filter(msmsdata, filterString)
# # 
# # # I need to recover all the PSMs for the confident peptides
# # # the re-evaluate PSM FDR
# # present.peptides <- unique(msmsdata.filtered$Peptide)
# # # recover peptides even if scores are low
# # msmsdata.filtered.2 <- subset(msmsdata, Peptide %in% present.peptides)
# # # view peptide FDR
# # eval_MSMS_ID_filter(msmsdata.filtered.2, filterString)
# # # check PSM FDR
# # psm_fdr(msmsdata.filtered.2)
# # # get rid of decoys and contaminants
# # # filterString <- paste(filterString, "!isDecoy", "!isContaminant", sep=" & ")
# # msmsdata.filtered.3 <- subset(msmsdata.filtered.2, !isDecoy & !isContaminant)

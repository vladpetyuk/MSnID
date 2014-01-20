

setClass(Class="MSnID",
         representation(
            workDir="character", # working directory
            psms="data.frame" # peptide-to-spectra matches
#             filter="MSnIDFilter" # keeps record of the applied filter
            )
         )




# Why does it have to be a seprate slot? Too complex for one?
# Or perhaps, filter manipulations imply that it exists as a separate class. Does it?
# 1) retrieving as string
# 2) setting/updating particular pars...
# ..a) how would you call a particular parameter? By name? *************
# e.g. msnidObj$MSnFilter$ppm.abs$threshold <- 10 ? Not sure about this way of accessing
# h.a? set_filter inteface? It is better, because there are multiple paths depending 
# on what is provided.
# Filter can exist within MSnIDObject, but get filter should return string (as well).
# Because I apply it as string.
# Do I need it as number?

# partitions and parameter transforms to be added to filter class
setClass(Class="MSnIDFilter",
         representation(
            filterList="list",
            validParNames="character"
         )
)


# #######  BELOW ARE THE THOUGHTS ################################################
# 
# 
# # should be able to set it like
# # filtObj <- set_filter(filtObj, name="abs.ppm", comparison="<", threshold=2)
# # filtObj <- set_filter(filtObj, name="abs.ppm") # acceptable too
# # filtObj <- set_filter(filtObj, name="MSGFScore") # acceptable too
# # filtObj <- optimize_filter_threshold( msnidObj, filtObj, FDR=0.05, level="peptide")
# # eval_filter (at given level)? to get FDR and number of peptides
# # msnid <- apply_filter( "MSnID", "character" | "MSnIDFilter")
# # filtStr <- as(filtObj, "character")
# 
# 
# # alternative
# # set_filter(msnidObj, name="abs.ppm")
# # set_filter(msnidObj, name="XCorr", append=TRUE)
# # set_filter(msnidObj, name="msmsScore", append=FALSE)
# # how about if I want to update the threshold is it append/replace
# # action = c("append", "reset", "update")?
# # how about rule:
# # if not exist - append,
# # if exist - update threshold (leave comparison if exists)
# # reset is separate function
# 
# # Why not list of lists? Dual interface for standalone class is fine
# # Because I want it act as a string.
# # The problem is with guessing directionality.
# # How about: msnidFilterObj <- guess_comparison(msnidFilterObj, msnidObj)
# 
# # string filter supercedes list? No, if set as string, then reset list
# # list reconstruct string.
# 
# # I REALLY WANT TO MAKE USE OF APPLY_FILTER( "CHARACTER" ) FUNCTION !!!
# # PERHAPS IT CAN TRACK RECORD OF "CHARACTER" FILTER OR FILTER OBJECT
# 
# 
# 
# # Assumption #1 I can define scoring/partitioning parameters
# # how to handle partitions?
# # Partitioning is based on numMissCleavages,
# # optimization is based on msmsScore and ppmError (after taking absolute)
# # e.g. (msmsScore > 9 & abs(ppmError) < 10 & numMissCleavages == 0) | 
# #      (msmsScore > 10 & abs(ppmError) < 8 & numMissCleavages == 1) |
# #      (msmsScore > 11 & abs(ppmError) < 5 & numMissCleavages == 2)
# #...
# # Let's ignore partitioning for now therefore
# # optimization is based on all 3 parameters
# # (msmsScore > 11 & abs(ppmError) < 5 & numMissCleavages  < 1)
# 
# # It may be nicer to have them as list
# # list(parName="abs.ppmError", relation="<", threshold=8.5)
# 
# # 
# # setClass(Class="MSnIDFilter",
# #          representation(
# #             partition.parameters="character",
# #             parameters="character",
# # #             transformation? e.g. abs ... or do the transformations before
# #             relational.operator="character",
# #             thresholds="numeric"
# #          )
# # )
# # 
# # 
# 
# 
# 
# 
# 
# 
# #fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
# # filtering
# # 1) ultimately I'd like to pass filter as a string
# # 2) however the values should be able to set programatically
# # Note.
# # fortune(106)
# alternative we can store filter as list
# list(parName1=c(operator, threshold), 
#      parName2=c(operator, threshold), 
#      parName3=c(operator, threshold))
# #
# or if we ensure that scoring parameters are higher the better
# operator becomes ">"
# list(parName1=threshold, 
#      parName2=threshold, 
#      parName3=threshold)
# #
# # concatenated with AND operator
# list(
#    list(name="character", 
#         comparison="character", 
#         threshold="numeric"),
#    ...
# )
# 
# ffs <- list(
#    list(name="msmsScore",
#         comparison=">",
#         threshold=9),
#    list(name="abs.parentMassErrorPPM",
#         comparison="<",
#         threshold=5))
# # convert ffs to string with "&" concatenation 
# paste(lapply(ffs, paste, collapse=' '), collapse=" & ")
# # (no template necessary?)
# # why is it ordered then?
# 
# ffs <- list(
#    list(name="msmsScore",
#         comparison=">",
#         threshold=9),
#    list(name="abs.parentMassErrorPPM",
#         threshold=5,
#         comparison="<"))
# # [c("name","comparison","threshold")]
# # this is more robust since it enforces order during string construction
# termsOrder <- c("name", "comparison", "threshold") # just to be sure
# termsOrder <- names(ffs[[1]]) # ultimately can be changed to any order
# paste(lapply(ffs, with, paste(mget(termsOrder), collapse=' ')), collapse=" & ")
# 
# 
# 
# # #--------- alternative ----------------------------------
# # ffs <- list(
# #    msmsScore=list(comparison=">",
# #                   threshold=9),
# #    abs.ppm=list(comparison="<",
# #                 threshold=5))
# # paste(lapply(names(ffs), function(v){
# #    paste(c(v, ffs[[v]]), collapse=' ')}), collapse=" & ")
# 
# 
# 
# 
# 
# 

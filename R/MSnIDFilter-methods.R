



setGeneric("get_starting_parameters",
           function(x, y, ...) standardGeneric("get_starting_parameters"))



setMethod("get_starting_parameters",
          signature(x="MSnIDFilter", y="MSnID"),
          definition=function(x, y, ...)
          {
             # let's assume we are taking median
             # update the filter list
             for( name_i in names(x@filterList)){
                x@filterList[[name_i]]$threshold <- 
                   median(y[[name_i]],na.rm=TRUE)
             }
             return(x)
          }
)




# # how about rule:
# # if not exist - append,
# # if exist - update threshold (leave comparison if exists)
# # reset is separate function

setMethod("set_filter", 
          signature(.Object="MSnIDFilter"), # dispatch just on the first argument
          definition=function(.Object, 
                              parName,          # although this MUST be character
                              comparison=NULL,
                              threshold=NULL)
          {
             parNames <- sapply(.Object@filterList, "[[", "name")
             parNames <- sapply(.Object@filterList, "[[", "name")
          }
)


.is_filterList_valid <- function(filterList)
{
   entries=c("comparison","threshold")
   all(sapply(filterList, function(v) all(entries %in% names(v))))
}



# internal
.get_filterString <- function(.Object, precision=2)
{
   # 1) assuming names in the list match names of list elements
   # 2) assuming the names are in the rigth order
   res <- lapply(.Object@filterList, 
                 within, 
                 {threshold <- signif(threshold, get("precision", parent.frame(n=4)))})
   res <- lapply(names(res), 
                 function(v) 
                    paste(c(list(v), res[[v]]), collapse=' '))
   res <- paste(sprintf("(%s)", res), collapse=" & ")
   return(res)
}


.get_filterValues <- function(.Object, precision=2)
{
   # 1) assuming names in the list match names of list elements
   # 2) assuming the names are in the rigth order
   res <- sapply(.Object@filterList, "[[", "threshold") 
   return(res)
}


setAs("MSnIDFilter", "character",
      def=function(from) .get_filterString(from, precision=Inf))


setAs("MSnIDFilter", "numeric",
      def=function(from) .get_filterValues(from, precision=Inf))

# getGeneric("as.numeric")
setMethod("as.numeric", "MSnIDFilter",
          definition=function(x, ...)
             as(x,"numeric"))

setMethod("as.vector", "MSnIDFilter",
          definition=function(x) as(x,"numeric"))

setMethod("length", "MSnIDFilter",
          definition=function(x) length(x@filterList))

setMethod("names", "MSnIDFilter",
          definition=function(x) names(x@filterList))

setMethod("update", "MSnIDFilter",
          definition=function(object, newThresholds)
          {
             #
             stopifnot(length(object) == length(newThresholds))
             for(i in 1:length(object))
                object@filterList[[i]]$threshold <- newThresholds[i]
             return(object)
          }
)
             




setMethod("show", "MSnIDFilter",
          definition=function(object)
          {
             cat("An object of calss \"",class(object),"\"\n",sep='')
             cat("Filter as string:\n")
             cat(.get_filterString(object),'\n')
             #cat(capture.output(str(filterList))[-1],sep='\n')
             
          }
)


MSnIDFilter <- function(MSnIDObj, filterList=list())
{
   # validations supposed to be *before* the return
   # let's allow only complete filters now
   if(!is.list(filterList)){
      stop(deparse(substitute(filterList)), 
           " is not a list!")
   }
   if(!.is_filterList_valid(filterList)){
      stop("Invalid (non-complete) entries in the filterList\n")
   }
   # update filterList here, so it includes names
   
   return(new("MSnIDFilter", 
              filterList=filterList,
              validParNames=names(MSnIDObj)))
}


# now let's define "$" operator

setMethod("$", "MSnIDFilter",
          definition=function(x, name)
          {
             # inspired by 
             # showMethods("$")
             # getMethod("$", "refObjectGenerator")
             #..
             # library(Biobase)
             # getMethod("$", "eSet")
             # eval(substitute(phenoData(x)$NAME_ARG, list(NAME_ARG = name)))
             # !!! works !!! # eval(substitute(x@filterList$name))
             return(x@filterList[[name]])
          }
)




setMethod("$<-", "MSnIDFilter",
          definition=function(x, name, value)
          {
             # inpired by
             # library(Biobase)
             # getMethod("$<-", "eSet")
             #---------
             # check if "value" is list
             # 
             if(!(name %in% x@validParNames)){
                stop(name, " is not a valid parameter!\n",
                     "See parameter names method for MSnID object.\n",
                     "Valid names are:\n",
                     paste(x@validParNames, collapse='\n'))
                
             }
             # c("comparison","threshold")
             # c("threshold","comparison")
             if(!is.list(value)){
                stop(value, " is not a list!")
             }
             if(!all(sort(names(value)) == c("comparison","threshold"))){
                stop(value, 
                     " names in the list do not match:\n", 
                     "comparison and threshold")
             }
             #
             comparison.operators <- c(">", ">=", "==", "<=", "<")
             if(!(value$comparison %in% comparison.operators)){
                stop(value$comparison, 
                     sprintf(" has to be one of %s!",
                             paste(comparison.operators, collapse=' ')))
             }
             if(!is.numeric(value$threshold)){
                stop(value$threshold, " is not numeric")
             }
             # to make sure they are ordered
             value <- value[c("comparison","threshold")]
             x@filterList[[name]] <- value
             return(x)
          }
)






# 
# ff = MSnIDFilter()
# print(ff)
# 
# filterList <- list(ppm=list(comparison=">", threshold=pi),
#                    score=list(comparison=">", threshold=0.00000000001),
#                    hz=list(comparison="<", threshold=2/3))
# ff=MSnIDFilter(filterList)
# print(ff)
# 
# filterList <- list(ppm=list(threshold=5),
#                    score=list(comparison=">", threshold=0.01),
#                    hz=list())
# # nice catch of error
# ff=MSnIDFilter(filterList)
# 
# 
# 
# filterList <- list(ppm=list(comparison=">", threshold=pi),
#                    score=list(comparison=">", threshold=0.00000000001),
#                    name=list(comparison='=', threshold=0),
#                    hz=list(comparison="<", threshold=2/3))
# ff=MSnIDFilter(filterList)
# print(ff)

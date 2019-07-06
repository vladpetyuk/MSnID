
.is_filterList_valid <- function(filterList)
{
    entries=c("comparison","threshold")
    all(sapply(filterList, function(v) all(entries %in% names(v))))
}



# internal
.get_filterString <- function(object, precision=2)
{
    # 1) assuming names in the list match names of list elements
    # 2) assuming the names are in the rigth order
    res <- lapply(object@filterList, 
                    within, 
                    {threshold <- signif(threshold, 
                                        get("precision", parent.frame(n=4)))})
    res <- lapply(names(res), 
                    function(v) 
                        paste(c(list(v), res[[v]]), collapse=' '))
    res <- paste(sprintf("(%s)", res), collapse=" & ")
    return(res)
}


.get_filterValues <- function(object, precision=2)
{
    # 1) assuming names in the list match names of list elements
    # 2) assuming the names are in the rigth order
    res <- sapply(object@filterList, "[[", "threshold") 
    return(res)
}


setAs("MSnIDFilter", "character",
        def=function(from) .get_filterString(from, precision=Inf))


setAs("MSnIDFilter", "numeric",
        def=function(from) .get_filterValues(from, precision=Inf))

setMethod("as.numeric", "MSnIDFilter",
            definition=function(x, ...) as(x,"numeric"))

setMethod("length", "MSnIDFilter",
            definition=function(x) length(x@filterList))

setMethod("names", "MSnIDFilter",
            definition=function(x) names(x@filterList))

setMethod("update", "MSnIDFilter",
            definition=function(object, newThresholds)
            {
                stopifnot(length(object) == length(newThresholds))
                for(i in seq_along(object))
                    object@filterList[[i]]$threshold <- newThresholds[i]
                return(object)
            }
)

setMethod("show", "MSnIDFilter",
            definition=function(object)
            {
                cat(class(object)," object\n",sep='')
                cat(.get_filterString(object),'\n')
            }
)


MSnIDFilter <- function(MSnIDObj, filterList=list())
{
    # validations supposed to be *before* the return
    # let's allow only complete filters now
    if(!is.list(filterList)){
        stop(deparse(substitute(filterList)), " is not a list!")
    }
    if(!.is_filterList_valid(filterList)){
        stop("Invalid (non-complete) entries in the filterList\n")
    }
    # update filterList here, so it includes names

    return(new( "MSnIDFilter", 
                filterList=filterList, 
                validParNames=names(MSnIDObj)))
}



setMethod("$", "MSnIDFilter",
            definition=function(x, name) x@filterList[[name]])




setMethod("$<-", "MSnIDFilter",
            definition=function(x, name, value)
            {
                # Has the parameter with the given name
                # been actually present in the MSnID object?
                if(!(name %in% x@validParNames)){
                    stop(name, " is not a valid parameter!\n",
                        "See parameter names method for MSnID object.\n",
                        "Valid names are:\n",
                        paste(x@validParNames, collapse='\n'))
                }

                # value must be a list of comparison operator 
                # and threshold value
                if(!is.list(value)){
                    stop(value, " is not a list!")
                }

                # do I actually have comparison and threshold names?
                if(!all(sort(names(value)) == c("comparison","threshold"))){
                    stop(value, 
                        " names in the list do not match:\n", 
                        "comparison and threshold")
                }

                # is the comparison operator valid?
                comparison.operators <- c(">", ">=", "==", "<=", "<")
                if(!(value$comparison %in% comparison.operators)){
                    stop(value$comparison, 
                        sprintf(" has to be one of %s!",
                                paste(comparison.operators, collapse=' ')))
                }

                # threshold must be numeric (at least currently)
                if(!is.numeric(value$threshold)){
                    stop(value$threshold, " is not numeric")
                }

                # to make sure comparison and threshold are in right order
                # and assign it to the filterList slot
                value <- value[c("comparison","threshold")]
                x@filterList[[name]] <- value
                return(x)
            }
)


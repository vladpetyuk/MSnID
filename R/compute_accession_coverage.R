

.compute_accession_coverage <- function(object,
                                      fasta,
                                      accession_col,
                                      pepSeq_col) {

  # check if fasta entry names are unique
  if(any(duplicated(names(fasta)))){
    message("FASTA entry names are not unique!\n")
    stop()
  }
  
  # check if there is at least some agreement in IDs
  if (!all(object[[accession_col]] %in% names(fasta))) {
    stop("Some accession IDs not found in FASTA entries!\n")
  }
  
  
  get_coverage_for_single_accession <- function(accession_i, ids, fasta) {
    
    
    x <- ids %>% 
      filter(accession == accession_i) %>%
      distinct()
    
    accAAstring <- fasta[[accession_i]]
    
    # main loop
    irl <- IRangesList()
    for(j in 1:nrow(x)) {
      pepSeq_j = x[[pepSeq_col]][j]
      mtch <- regexpr(pepSeq_j, accAAstring)
      start <- as.numeric(mtch)
      width <- attr(mtch, "match.length")
      tgt <- IRanges(start=start, width=width, names=pepSeq_j)
      irl[[j]] <- tgt
    }
    
    fullAACoverage <- sum(IRanges::width(reduce(unlist(irl))))
    percentAACoverage <- 100*fullAACoverage/length(accAAstring)
    return(percentAACoverage)
  }
  
  ids <- object %>% 
    select(accession_col, pepSeq_col) %>%
    distinct()
  
  res <- ids %>% 
    select(accession_col) %>%
    distinct()
  
  res$percentAACoverage <-  map(res[[accession_col]],
                                get_coverage_for_single_accession,
                                ids, fasta)
  res$percentAACoverage <- as.numeric(res$percentAACoverage)
  
  res <- inner_join(object, res, by=accession_col)
  
  return(res)
}

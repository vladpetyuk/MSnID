
utils::globalVariables("seq_length")

remap_fasta_entry_names <- function(path_to_FASTA,
                                    conversion_table,
                                    extraction_pttrn=c("\\|([^|-]+)(-\\d+)?\\|",
                                                       "([A-Z]P_\\d+)",
                                                       "(ENS[A-Z0-9]+)")){
    
    is_compressed <- FALSE
    if(grepl("[.]gz$", path_to_FASTA)){
        is_compressed <- TRUE
    }else if(grepl("[.](bz2|xz|zip)$", path_to_FASTA)){
        stop("The only supported compression is gzip!")
    }
    
    mySequences <- readAAStringSet(path_to_FASTA)
    extraction_pttrn <- paste0(".*", extraction_pttrn, ".*")
    names(mySequences) <- sub(extraction_pttrn, "\\1", names(mySequences))
    prot_lengths <- data.frame(seq_name = names(mySequences),
                               seq_length = width(mySequences),
                               stringsAsFactors = FALSE)
    
    by_ <- "seq_name"
    names(by_) <- colnames(conversion_table)[1]
    conversion_table <- conversion_table %>%
        inner_join(prot_lengths, by = by_) %>%
        group_by_at(colnames(conversion_table)[2]) %>%
        dplyr::slice(which.max(seq_length))
    
    # in case there are multiple "to" IDs, the named vector will return
    # just the first one anyway
    conversion_table <- conversion_table %>%
        select(-seq_length) %>%
        distinct() %>%
        group_by_at(colnames(conversion_table)[1]) %>%
        top_n(1,dplyr::desc(!!as.name(colnames(conversion_table)[2])))
    
    conversion_vec <- conversion_table %>% pull(2)
    names(conversion_vec) <- conversion_table %>% pull(1)
    
    mySequences <- mySequences[names(conversion_vec)]
    names(mySequences) <- conversion_vec[names(mySequences)]
    
    file_no_ext <- tools::file_path_sans_ext(path_to_FASTA, compression=TRUE)
    ext <- sub(file_no_ext, "", path_to_FASTA, fixed=TRUE)
    
    path_to_FASTA_remapped <- paste0(file_no_ext, '_remapped', ext)
    
    writeXStringSet(mySequences, path_to_FASTA_remapped, compress = is_compressed)
    return(path_to_FASTA_remapped)
}





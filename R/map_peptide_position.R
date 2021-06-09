
utils::globalVariables(c(".", "trimmedPeptide", "x",
                         "First_AA", "Last_AA",
                         "ProtSeq", "cleanSeq", "PepLoc", "ModShift",
                         "SiteLoc", "ModAAs", "SiteLoc", "ModAAs",
                         "Site", "SiteCollapsed", "SiteCollapsedFirst"))


.get_non_AA_characters <- function(object, 
                                  column_name, 
                                  extra_allowed_chars = ""){
    # returns any non AA characters from the presumably peptide sequences
    present_chars <- object[[column_name]] %>%
        unique() %>%
        paste0(collapse = '') %>% 
        strsplit(split='') %>% 
        `[[`(1) %>% 
        unique()
    other_chars <- setdiff(present_chars, 
                           c(Biostrings::AA_STANDARD, extra_allowed_chars))
    return(other_chars)
}



.check_presense_of_flanking_AAs <- function(object, column_name){
    # TRUE if peptides in X.XXXXXXX.X format
    flanking_AA_present <- object[[column_name]] %>%
        unique() %>%
        grepl(".\\..+\\..", .) %>%
        all()
    return(flanking_AA_present)
}


# .make_clean_seq_old <- function(object, peptide_seq_column = "peptide"){
#     
#     flanking_AAs_present <- .check_presense_of_flanking_AAs(object, 
#                                                            peptide_seq_column)
#     if(flanking_AAs_present)
#         object$cleanSeq <- sub(".\\.(.+)\\..", "\\1", object[[peptide_seq_column]])
#     else
#         object$cleanSeq <- object[[peptide_seq_column]]
#     
#     non_AA_chars <- .get_non_AA_characters(object, "cleanSeq")
#     
#     if(length(non_AA_chars) > 0){
#         non_AA_chars_pttrn <- non_AA_chars %>% 
#             # map_chr(~paste0("\\",.x)) %>% # no need for \\
#             paste0(collapse='') %>% 
#             paste0("[",.,"]")
#         object$cleanSeq <- str_replace_all(object$cleanSeq, non_AA_chars_pttrn, "")
#     }
#     
#     return(object)
# }


# .make_clean_seq <- function(object, peptide_seq_column = "peptide"){
#     
#     flanking_AAs_present <- .check_presense_of_flanking_AAs(object, 
#                                                            peptide_seq_column)
#     if(flanking_AAs_present){
#         cleanSeq <- sub(".\\.(.+)\\..", "\\1", object[[peptide_seq_column]])
#     }else{
#         cleanSeq <- object[[peptide_seq_column]]
#     }
#     
#     non_AA_chars <- .get_non_AA_characters(object, "cleanSeq")
#     if(length(non_AA_chars) > 0){
#         non_AA_chars_pttrn <- non_AA_chars %>% 
#             paste0(collapse='') %>% 
#             paste0("[",.,"]")
#         cleanSeq <- str_replace_all(cleanSeq, non_AA_chars_pttrn, "")
#     }
#     
#     # assign
#     object_psms <- copy(object@psms) # preventing side-effect on original object
#     object_psms <- data.table(object_psms) # safety step to avoid "Invalid .internal.selfref"
#     object_psms[, cleanSeq := cleanSeq]
#     object@psms <- object_psms
#     
#     return(object)
# }



.make_clean_seq <- function(object, peptide_seq_column = "peptide"){
    
    flanking_AAs_present <- .check_presense_of_flanking_AAs(object, 
                                                           peptide_seq_column)
    if(flanking_AAs_present)
        object$cleanSeq <- sub(".\\.(.+)\\..", "\\1", object[[peptide_seq_column]])
    else
        object$cleanSeq <- object[[peptide_seq_column]]
    
    non_AA_chars <- .get_non_AA_characters(object, "cleanSeq")
    
    if(length(non_AA_chars) > 0){
        non_AA_chars_pttrn <- non_AA_chars %>% 
            paste0(collapse='') %>% 
            paste0("[",.,"]")
        object$cleanSeq <- str_replace_all(object$cleanSeq, non_AA_chars_pttrn, "")
    }
    
    return(object)
}




# .map_peptide_position_old <- function(object, fasta, accession_col = "accession"){
#     
#     object <- .make_clean_seq(object, "peptide")
#     x <- psms(object)
#     
#     # protein ID in `accession`
#     # peptide sequence with flanking AAs in `peptide`
#     
#     # check for decoy entries
#     decoy_acc <- apply_filter(object, "isDecoy")[[accession_col]]
#     if (length(decoy_acc) > 0 & !any(decoy_acc %in% names(fasta))) {
#         fasta_rev <- reverse(fasta)
#         names(fasta_rev) <- paste0("XXX_", names(fasta))
#         fasta <- c(fasta, fasta_rev)
#     }
#     
#     # check if fasta entry names are unique
#     if(any(duplicated(names(fasta)))){
#         stop("FASTA entry names are not unique!\n")
#     }
#     
#     # check if there is at least some agreement in IDs
#     if(length(intersect(x[[accession_col]], names(fasta))) == 0){
#         stop("There is zero overlap in protein IDs and FASTA entry names!\n")
#     }
#     
#     # merger of identifications and FASTA
#     res <- fasta %>%
#         as.data.frame() %>%
#         rownames_to_column(accession_col) %>%
#         dplyr::rename(ProtSeq = x) %>%
#         mutate(ProtLen = str_length(ProtSeq)) %>%
#         inner_join(x, ., by = accession_col)
#     
#     # locating peptide within protein
#     res <- res %>%
#         mutate(First_AA = map2(ProtSeq, 
#                                cleanSeq, 
#                                ~ as.numeric(str_locate_all(.x, .y)[[1]][,1])),
#                Last_AA = map2(First_AA, nchar(cleanSeq) - 1, `+`),
#                First_AA_First = map(First_AA, ~ .[1]) %>% as.numeric(),
#                Last_AA_First = map(Last_AA, ~ .[1]) %>% as.numeric())
#     
#     # drop Protein Sequences
#     res <- res %>%
#         select(-ProtSeq)
#     
#     psms(object) <- as.data.frame(res)
#     
#     return(object)
# }




# .map_peptide_position_old2 <- function(object, fasta, accession_col = "accession"){
#     
#     object <- .make_clean_seq(object, "peptide")
#     
#     # check for decoy entries
#     decoy_acc <- apply_filter(object, "isDecoy")[[accession_col]]
#     if (length(decoy_acc) > 0 & !any(decoy_acc %in% names(fasta))) {
#         fasta_rev <- reverse(fasta)
#         names(fasta_rev) <- paste0("XXX_", names(fasta))
#         fasta <- c(fasta, fasta_rev)
#     }
#     
#     # check if fasta entry names are unique
#     if(any(duplicated(names(fasta)))){
#         stop("FASTA entry names are not unique!\n")
#     }
#     
#     # check if there is at least some agreement in IDs
#     if(length(intersect(object[[accession_col]], names(fasta))) == 0){
#         stop("There is zero overlap in protein IDs and FASTA entry names!\n")
#     }
#     
#     # main loop
#     First_AA <- list()
#     Last_AA <- list()
#     ProtLen <- numeric()
#     for(i in seq_len(nrow(object))){
#         #
#         #
#         if(i %% 1000 == 0) # remove
#             print(i) # remove
#         #
#         #
#         prot_seq_i <- fasta[[object$accession[i]]]
#         if(is.null(prot_seq_i)){
#             x <- NA
#             ProtLen <- c(ProtLen, NA)
#         }else{
#             clean_seq_i <- object$cleanSeq[i]
#             # x <- as.numeric(str_locate_all(prot_seq, clean_seq)[[1]][,1])
#             x <- as.numeric(gregexpr(clean_seq_i, prot_seq_i)[[1]])
#             if(x[1] == -1) 
#                 x <- NA
#             ProtLen <- c(ProtLen, length(prot_seq_i))
#         }
#         First_AA <- c(First_AA, list(x))
#         # add peptide length
#         Last_AA <- c(Last_AA, list(x + nchar(clean_seq_i) - 1))
#     }
#     
#     # assign back to the MSnID object
#     object$ProtLen <- ProtLen
#     object$First_AA <- First_AA
#     object$Last_AA <- Last_AA
#     # derive just the first location
#     object$First_AA_First <- sapply(object$First_AA, "[", 1)
#     object$Last_AA_First <- sapply(object$Last_AA, "[", 1)
#     
#     return(object)
# }





.map_peptide_position <- function(object, fasta, accession_col = "accession"){
    
    object <- .make_clean_seq(object, "peptide")
    x <- psms(object)
    
    # protein ID in `accession`
    # peptide sequence with flanking AAs in `peptide`
    
    # check for decoy entries
    decoy_acc <- apply_filter(object, "isDecoy")[[accession_col]]
    if (length(decoy_acc) > 0 & !any(decoy_acc %in% names(fasta))) {
        fasta_rev <- reverse(fasta)
        names(fasta_rev) <- paste0("XXX_", names(fasta))
        fasta <- c(fasta, fasta_rev)
    }
    
    # check if fasta entry names are unique
    if(any(duplicated(names(fasta)))){
        stop("FASTA entry names are not unique!\n")
    }
    
    # check if there is at least some agreement in IDs
    if(length(intersect(x[[accession_col]], names(fasta))) == 0){
        stop("There is zero overlap in protein IDs and FASTA entry names!\n")
    }
    
    # merger of identifications and FASTA
    prot_pep <- x %>%
        select(!!accession_col, cleanSeq) %>%
        distinct()
    
    res <- fasta %>%
        as.data.frame() %>%
        rownames_to_column(accession_col) %>%
        dplyr::rename(ProtSeq = x) %>%
        mutate(ProtLen = str_length(ProtSeq)) %>%
        left_join(prot_pep, ., by = accession_col)
    
    
    # locating peptide within protein
    res <- res %>%
        mutate(First_AA = map2(ProtSeq, 
                               cleanSeq, 
                               ~ as.numeric(str_locate_all(.x, .y)[[1]][,1])),
               Last_AA = map2(First_AA, nchar(cleanSeq) - 1, `+`),
               First_AA_First = map(First_AA, ~ .[1]) %>% as.numeric(),
               Last_AA_First = map(Last_AA, ~ .[1]) %>% as.numeric())
    
    # drop Protein Sequences
    res <- res %>%
        select(-ProtSeq)
    
    # linking back to the main MSnID object by accession and cleanSeq
    # the linking is two-step to avoid problems with RAM
    res_full <- x %>%
        select(!!accession_col, cleanSeq) %>%
        left_join(res, by=c(accession_col, "cleanSeq"))
    
    # clean-up
    rm(x)
    gc()
    
    # VERSION 1
    columns_to_add <- c("First_AA", "Last_AA", "First_AA_First", "Last_AA_First", "ProtLen")
    
    # assign
    object_psms <- copy(object@psms) # preventing side-effect on original object
    object_psms <- data.table(object_psms) # safety step to avoid "Invalid .internal.selfref"
    object_psms[, (columns_to_add) := res_full[,columns_to_add]]
    object@psms <- object_psms
    
    return(object)
}



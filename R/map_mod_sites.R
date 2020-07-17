

utils::globalVariables(c(".", "TrimmedPeptide", "x",
                         "ProtSeq", "CleanSeq", "PepLoc", "ModShift",
                         "SiteLoc", "ModAAs", "SiteLoc", "ModAAs",
                         "Site", "SiteCollapsed", "SiteCollapsedFirst"))

.map_mod_sites <- function(ids, 
                           fasta, 
                           accession_col, 
                           peptide_mod_col, 
                           mod_char){
    
    # check if fasta entry names are unique
    if(any(duplicated(names(fasta)))){
        message("FASTA entry names are not unique!\n")
        stop()
    }
    
    # check if there is at least some agreement in IDs
    if(length(intersect(ids[[accession_col]], names(fasta))) == 0){
        message("There is zero overlap in protein IDs and FASTA entry names!\n")
        stop()
    }
    
    # check the characters
    # there should be nothing except the AAs and modification character
    # ids <- ids %>%
    #    mutate_(TrimmedPeptide := sub(".\\.(.*)\\..", "\\1", peptide_mod_col))
    ids <- ids %>%
        mutate(TrimmedPeptide = sub(".\\.(.*)\\..", "\\1", .[[peptide_mod_col]]))
    present_chars <- paste0(ids$TrimmedPeptide, collapse = '') %>% 
        strsplit(split='') %>% 
        `[[`(1) %>% 
        unique()
    other_chars <- setdiff(present_chars, c(AA_STANDARD, mod_char))
    if(length(other_chars) > 0){
        message("Detected extra chararacters in the peptide sequences!\n")
        # erase other chars in TrimmedPeptide
        other_chars_pttrn <- other_chars %>% 
            map_chr(~paste0("\\",.x)) %>% 
            paste0(collapse='') %>% 
            paste0("[",.,"]")
        ids <- ids %>%
            mutate(TrimmedPeptide = str_replace_all(TrimmedPeptide, other_chars_pttrn, ""))
    }
    
    # check for missing peptide sequences
    if(any(is.na(ids$TrimmedPeptide))){
        message("Entries with missing peptide sequences were detected. 
The correspoding identifications will be removed.\n")
        ids <- ids %>% filter(!is.na(TrimmedPeptide))
    }
    
    # check if there are peptides without mods
    mod_char_pttrn <- paste0("[\\", mod_char, "]")
    if(any(!str_detect(ids$TrimmedPeptide, mod_char_pttrn))){
        message("Peptides with no modifications in question were detected. 
The correspoding identifications will be removed.\n")
        ids <- ids %>%
            filter(str_detect(TrimmedPeptide, mod_char_pttrn))
    }
    
    # extract clean sequence
    mod_char_pttrn <- paste0("[\\", mod_char, "]")
    ids <- ids %>%
        mutate(CleanSeq = str_replace_all(TrimmedPeptide, mod_char_pttrn, ""))
    
    # merger of identifications and FASTA
    res <- fasta %>%
        as.data.frame() %>%
        rownames_to_column(accession_col) %>%
        rename(ProtSeq = x) %>%
        inner_join(ids, ., by = accession_col)
    
    # the core part
    res <- res %>%
        # locating peptide within protein
        mutate(PepLoc = map2(ProtSeq, CleanSeq, ~ as.numeric(str_locate_all(.x, .y)[[1]][,1])),
               PepLocFirst = map(PepLoc, ~ .[1])) %>%
        # adding protein length
        mutate(ProtLength = str_length(ProtSeq)) %>%
        # drop protein sequence
        select(-ProtSeq) %>%
        # locations of PTM within peptide
        mutate(ModShift = map(TrimmedPeptide, ~ as.numeric(str_locate_all(.,mod_char_pttrn)[[1]][,1])),
               ModShift = map(ModShift, ~ . - seq_along(.) - 1)) %>%
        # extract AAs
        mutate(ModAAs = map2(CleanSeq, ModShift, ~ str_sub(.x, .y+1, .y+1))) %>%
        # calculate site positions: peptide location + mod shift
        mutate(SiteLoc = map2(PepLoc, ModShift, ~ lapply(.x, `+`, .y))) %>%
        # create site notation: aa + location
        mutate(Site = map2(SiteLoc, ModAAs, ~ lapply(.x, function(..) paste0(.y, ..)))) %>%
        # collapsing site notation
        mutate(SiteCollapsed = map(Site, ~ lapply(.x, paste0, collapse =','))) %>%
        # picking the first one
        mutate(SiteCollapsedFirst = map(SiteCollapsed, 1))
    
    # clean-up
    res <- res %>%
        select(-c(TrimmedPeptide,CleanSeq))
    
    return(res)
    
}




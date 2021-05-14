utils::globalVariables(c(".", "trimmedPeptide", "x",
                         "ProtSeq", "CleanSeq", "PepLoc", "ModShift",
                         "SiteLoc", "ModAAs", "SiteLoc", "ModAAs",
                         "Site", "SiteCollapsed", "SiteCollapsedFirst",
                         "sequenceWindow"))


.extract_sequence_window <- function (ids, fasta,
                                     accession_col="accession",
                                     site_loc_col="SiteLoc",
                                     radius=7L,
                                     collapse="|") {
    
    if (!("SiteID" %in% names(ids))) {
        stop("No SiteID found. Call map_mod_sites.")
    }
    
    x <- ids %>%
        select(accession_col, site_loc_col) %>%
        distinct()
    
    x <- fasta %>%
        as.data.frame() %>%
        rownames_to_column(accession_col) %>%
        rename(ProtSeq = x) %>%
        mutate(ProtSeqWidth = nchar(ProtSeq)) %>%
        inner_join(x, ., by=accession_col)
    
    
    f <- function(ProtSeq_i, SiteLoc_i) {
        sequenceWindows <- c()
        # take first occurence of site
        for (k in unlist(SiteLoc_i[[1]])) {
            # get 7 amino acids to the left
            site_left <- substr(ProtSeq_i, max(k-radius,1), k-1)
            # append "-" characters if went over start of protein sequence
            if (k-radius < 1) {
                site_left <- paste0(paste(rep("-", 1-(k-radius)),collapse=""), site_left)
            }
            # get 7 amino acids to the left
            site_right <- substr(ProtSeq_i, k+1, min(k+radius,nchar(ProtSeq_i)))
            # append "-" characters if went over end of protein sequence
            if (k+radius > nchar(ProtSeq_i)) {
                site_right <- site_right <- paste0(site_right, paste(rep("-", k+radius-nchar(ProtSeq_i)),collapse=""))
            }
            # turn phosphorylated amino acid to lowercase
            mod_aa <- tolower(substr(ProtSeq_i, k, k))
            # paste together -7 AA, the modified AA, and +7 AA.
            window <- paste0(site_left, tolower(mod_aa), site_right)
            # add to list
            sequenceWindows <- c(sequenceWindows, window)
        }
        sequenceWindows <- paste(sequenceWindows, collapse=collapse)
        return(sequenceWindows)
    }
    
    
    x$sequenceWindow <-  map2(x$ProtSeq, x[[site_loc_col]], f)
    x$sequenceWindow <- as.character(x$sequenceWindow)
    
    x <- x %>% select(accession_col, site_loc_col, sequenceWindow)
    
    x <- inner_join(ids, x, by=c(accession_col, site_loc_col))
    
    return(x)
}





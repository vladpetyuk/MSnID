

utils::globalVariables(c(".", "TrimmedPeptide", "x",
                         "ProtSeq", "CleanSeq", "PepLoc", "ModShift",
                         "SiteLoc", "ModAAs", "SiteLoc", "ModAAs",
                         "Site", "SiteCollapsed", "SiteCollapsedFirst"))


.map_flanking_sequences <- function (msnid, fasta, radius=7L, collapse="|") {
    
    if (!("SiteID" %in% names(msnid))) {
        stop("No SiteID found. Call map_mod_sites.")
    }
    
    x <- psms(msnid) %>%
        select(accession, SiteLoc) %>%
        distinct()
    
    x <- fasta %>%
        as.data.frame() %>%
        rownames_to_column("accession") %>%
        mutate(accession = sub("^(.P_\\d+\\.\\d+)?\\s.*", "\\1", accession)) %>%
        rename(ProtSeq = x) %>%
        mutate(ProtSeqWidth = nchar(ProtSeq)) %>%
        inner_join(x, ., by="accession")
    
    
    f <- function(ProtSeq_i, SiteLoc_i) {
        flankingSequences <- c()
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
            flank <- paste0(site_left, tolower(mod_aa), site_right)
            # add to list
            flankingSequences <- c(flankingSequences, flank)
        }
        flankingSequences <- paste(flankingSequences, collapse=collapse)
        return(flankingSequences)
    }
    
    
    x$flankingSequence <-  map2(x$ProtSeq, x$SiteLoc, f)
    x$flankingSequence <- as.character(x$flankingSequence)
    
    msnid@psms <- psms(msnid) %>%
        mutate(flankingSequence=NULL) %>%
        left_join(x, by=c("accession", "SiteLoc")) %>% 
        data.table()
    
    return(msnid)
}





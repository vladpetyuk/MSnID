

setMethod("remap_accessions", "MSnID",
          definition=function(object, 
                              conversion_table, 
                              extraction_pttrn, 
                              path_to_FASTA)
          {
              .remap_accessions(object, 
                                conversion_table, 
                                extraction_pttrn, 
                                path_to_FASTA)
          }
)


.remap_accessions <- function(object, 
                             conversion_table,
                             extraction_pttrn,
                             path_to_FASTA){
    
    conv_vec <- conversion_table[,2]
    names(conv_vec) <- conversion_table[,1] # should be msnid object accessions
    
    # if path to FASTA present,
    # retain only those [,2] entries that are present in the FASTA
    if(!is.null(path_to_FASTA)){
        accessions_in_fts <- names(readAAStringSet(path_to_FASTA))
        idx <- conv_vec %in% accessions_in_fts
        conv_vec <- conv_vec[idx]
    }
    
    extraction_pttrn <- paste0(".*", extraction_pttrn, ".*")
    object$accession <- sub(extraction_pttrn, "\\1", object$accession)
    
    # checking that the conversion table covers the accessions
    cvrg <- object$accession %in% names(conv_vec)
    cvrg <- factor(cvrg, levels = c(TRUE, FALSE))
    if(prop.table(table(cvrg))['FALSE'] > 0.5){
        stop("Majority of accessions not in the conversion_table!")
    }
    
    object$accession <- conv_vec[object$accession]
    
    # get rid of entries without annotation
    object <- apply_filter(object, "!is.na(accession)")
    
    # make sure decoy accessions start with XXX
    object$accession <- ifelse(object$isDecoy,
                              paste0("XXX_", object$accession),
                              object$accession)
    
    return(object)
}


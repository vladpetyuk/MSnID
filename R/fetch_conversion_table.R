
utils::globalVariables("REFSEQ")


fetch_conversion_table <- function(organism_name, from, to, backend="AnnotationHub", snapshot_date=NULL){
    
    wrong_column_name <- function(col_name){
        message(paste0(col_name," is not a valid column name in selected OrgDb."))
        message("Valid from/to names:")
        message(paste(AnnotationDbi::columns(db), collapse = ", "))
        stop()
    }
    
    ah <- AnnotationHub()
    
    if(!is.null(snapshot_date)){
        if(snapshot_date %in% AnnotationHub::possibleDates(ah)){
            AnnotationHub::snapshotDate(ah) <- snapshot_date
        }else{
            message(paste0(snapshot_date, " is not valid snapshot date"))
            message("valid dates:")
            message(paste(AnnotationHub::possibleDates(ah), collapse = ", "))
            stop()
        }
    }
    
    orgs <- subset(ah, ah$rdataclass == "OrgDb")
    
    if(!(organism_name %in% orgs$species)){
        stop(paste0(organism_name," is not in AnnotationHub database."))
    }
    
    db <- query(orgs, organism_name)
    db <- db[[1]]
    
    if(!(from %in% AnnotationDbi::columns(db))){
        wrong_column_name(from)
    }
    if(!(to %in% AnnotationDbi::columns(db))){
        wrong_column_name(to)
    }
    
    # main call
    conversion_table <- AnnotationDbi::select(db,
                                              keys = AnnotationDbi::keys(db),
                                              columns = c(from,to),
                                              keytype = "ENTREZID")
    conversion_table <- unique(conversion_table[,c(from,to)])
    conversion_table <- conversion_table[complete.cases(conversion_table),]
    
    # rectify the fact that REFSEQ returns both transcripts and proteins
    # which is a mess
    if("REFSEQ" %in% c(from, to)){
        conversion_table <- conversion_table %>%
            filter(grepl("^.P_.*", REFSEQ))
    }
    
    return(conversion_table)
}




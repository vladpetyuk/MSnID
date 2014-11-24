

setGeneric("assess_missed_cleavages",
            function(object, missedCleavagePattern="[KR](?=[^P$])") 
                standardGeneric("assess_missed_cleavages"))

setGeneric("assess_termini", 
            function(object, validCleavagePattern="[KR]\\.[^P]") 
                standardGeneric("assess_termini"))

setGeneric("accessions", 
            function(object) standardGeneric("accessions"))

setGeneric("proteins", 
            function(object) standardGeneric("proteins"))

setGeneric("peptides", 
            function(object) standardGeneric("peptides"))

setGeneric("read_mzIDs", 
            function(object, mzids) standardGeneric("read_mzIDs"))

setGeneric("apply_filter", 
            function(msnidObj, filterObj) standardGeneric("apply_filter"))

setGeneric("evaluate_filter", 
            function(object, filter, level=c("PSM", "peptide", "accession"))
                standardGeneric("evaluate_filter"))

setGeneric("id_quality",
            function(object, 
                        filter=NULL, 
                        level=c("PSM", "peptide", "accession")) 
                standardGeneric("id_quality"))

setGeneric("correct_peak_selection", 
            function(object) 
                standardGeneric("correct_peak_selection"))

setGeneric("mass_measurement_error",
            function(object) 
                standardGeneric("mass_measurement_error"))

setGeneric("recalibrate",
            function(object) 
                standardGeneric("recalibrate"))

setGeneric("optimize_filter",
            function(filterObj, msnidObj, fdr.max, method, level, n.iter)
                stadardGeneric("optimize_filter"))

setGeneric("psms<-",
            function(object, value) standardGeneric("psms<-"))

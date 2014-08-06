

setGeneric("assess_missed_cleavages",
            function(.Object, missedCleavagePattern="[KR](?=[^P$])") 
                standardGeneric("assess_missed_cleavages"))

setGeneric("assess_termini", 
            function(.Object, validCleavagePattern="[KR]\\.[^P]") 
                standardGeneric("assess_termini"))

setGeneric("accessions", 
            function(.Object) standardGeneric("accessions"))

setGeneric("peptides", 
            function(.Object) standardGeneric("peptides"))

setGeneric("read_mzIDs", 
            function(object, mzids) standardGeneric("read_mzIDs"))

setGeneric("apply_filter", 
            function(.Object, .Filter) standardGeneric("apply_filter"))

setGeneric("evaluate_filter", 
            function(.Object, filter, level=c("PSM", "Peptide", "Accession"))
                standardGeneric("evaluate_filter"))

setGeneric("id_quality",
            function(.Object, 
                        filter=NULL, 
                        level=c("PSM", "Peptide", "Accession")) 
                standardGeneric("id_quality"))

setGeneric("correct_peak_selection", 
            function(.Object) 
                standardGeneric("correct_peak_selection"))

setGeneric("mass_measurement_error",
            function(.Object) 
                standardGeneric("mass_measurement_error"))

setGeneric("recalibrate",
            function(.Object) 
                standardGeneric("recalibrate"))

setGeneric("optimize_filter",
            function(.Filter, .Data, fdr.max, method, level, n.iter)
                stadardGeneric("optimize_filter"))

setGeneric("psms",
            function(.Object) standardGeneric("psms"))

setGeneric("psms<-",
            function(.Object, value) standardGeneric("psms<-"))

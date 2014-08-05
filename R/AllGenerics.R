
setGeneric("set_filter",
            function(.Object, ...) standardGeneric("set_filter"))

setGeneric("assess_missed_cleavages",
            function(.Object, missedCleavagePattern="[KR](?=[^P$])") 
                standardGeneric("assess_missed_cleavages"))

# this may be an odd way to handle assess_termini
setGeneric("assess_termini", 
            function(.Object, validCleavagePattern="[KR]\\.[^P]") 
                standardGeneric("assess_termini"))

setGeneric("get_accession_fdr", 
            function(.Object, ...) standardGeneric("get_accession_fdr"))

setGeneric("get_peptide_fdr", 
            function(.Object, ...) standardGeneric("get_peptide_fdr"))

setGeneric("get_psm_fdr", 
            function(.Object, ...) standardGeneric("get_psm_fdr"))

setGeneric("get_accessions", 
            function(.Object) standardGeneric("get_accessions"))

setGeneric("get_peptides", 
            function(.Object) standardGeneric("get_peptides"))


# setGeneric("read_mzIDs", 
#            function(.Object, mzids) standardGeneric("read_mzIDs"))

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

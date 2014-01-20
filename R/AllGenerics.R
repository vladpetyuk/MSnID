
# setGeneric("parameter_names",
#            function(.Object) standardGeneric("parameter_names"))

setGeneric("set_filter",
           function(.Object, ...) standardGeneric("set_filter"))

setGeneric("delete_psm_parameter",
           function(.Object, ...) standardGeneric("delete_psm_parameter"))

setGeneric("get_psm_parameter",
           function(.Object, parName) standardGeneric("get_psm_parameter"))

setGeneric("set_psm_parameter",
           function(.Object, ...) standardGeneric("set_psm_parameter"))

setGeneric("assess_missed_cleavages",
           function(.Object, ...) standardGeneric("assess_missed_cleavages"))

# this may be an odd way to handle assess_termini
setGeneric("assess_termini", 
           function(.Object, ...) standardGeneric("assess_termini"))

setGeneric("get_accession_fdr", 
           function(.Object, ...) standardGeneric("get_accession_fdr"))

setGeneric("get_peptide_fdr", 
           function(.Object, ...) standardGeneric("get_peptide_fdr"))

setGeneric("get_psm_fdr", 
           function(.Object, ...) standardGeneric("get_psm_fdr"))

setGeneric("get_accessions", 
           function(.Object, ...) standardGeneric("get_accessions"))

setGeneric("get_peptides", 
          function(.Object, ...) standardGeneric("get_peptides"))


setGeneric("read_mzIDs", 
           function(.Object, ...) standardGeneric("read_mzIDs"))


setGeneric("apply_filter", 
           function(.Object, .Filter, ...) standardGeneric("apply_filter"))


setGeneric("evaluate_filter", 
           function(.Object, ...) standardGeneric("evaluate_filter"))
setGeneric("id_quality",
           function(.Object, ...) standardGeneric("id_quality"))



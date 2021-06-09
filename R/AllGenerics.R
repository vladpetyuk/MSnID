
setGeneric("remap_accessions",
           function(object, 
                    conversion_table, 
                    extraction_pttrn=c("\\|([^|-]+)(-\\d+)?\\|",
                                       "([A-Z]P_\\d+)",
                                       "(ENS[A-Z0-9]+)"), 
                    path_to_FASTA=NULL)
               standardGeneric("remap_accessions"))

setGeneric("add_mod_symbol",
           function(object, mod_mass, symbol)
               standardGeneric("add_mod_symbol"))


setGeneric("report_mods",
           function(object, ...) 
               standardGeneric("report_mods"))


setGeneric("map_mod_sites",
           function(object,
                    fasta,
                    accession_col = "accession",
                    peptide_mod_col = "peptide_mod",
                    mod_char = "*",
                    site_delimiter = "lower")
               standardGeneric("map_mod_sites"))


setGeneric("map_peptide_position",
           function(object,
                    fasta,
                    accession_col = "accession")
               standardGeneric("map_peptide_position"))


setGeneric("extract_sequence_window",
           function(object, fasta,
                    accession_col="accession",
                    site_loc_col="SiteLoc",
                    radius=7L,
                    collapse="|")
             standardGeneric("extract_sequence_window"))

setGeneric("compute_accession_coverage",
           function(object,
                    fasta,
                    accession_col="accession",
                    pepSeq_col="pepSeq")
             standardGeneric("compute_accession_coverage"))


setGeneric("infer_parsimonious_accessions",
           function(object, unique_only=FALSE, prior=character(0)) 
               standardGeneric("infer_parsimonious_accessions"))


setGeneric("assess_missed_cleavages",
            function(object, missedCleavagePattern="[KR](?=[^P$])") 
                standardGeneric("assess_missed_cleavages"))

setGeneric("assess_termini", 
            function(object, validCleavagePattern="[KR]\\.[^P]") 
                standardGeneric("assess_termini"))


setGeneric("peptides", 
            function(object) standardGeneric("peptides"))

setGeneric("read_mzIDs", 
            function(object, mzids, backend=c('mzID','mzR')) 
                standardGeneric("read_mzIDs"))

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
            function(filterObj, msnidObj, 
                     fdr.max, method, level, n.iter, mc.cores=NULL)
                standardGeneric("optimize_filter"))

setGeneric("psms<-",
            function(object, value) standardGeneric("psms<-"))


setGeneric("plot_protein_coverage",
           function(object, accession, ...) 
               standardGeneric("plot_protein_coverage"))







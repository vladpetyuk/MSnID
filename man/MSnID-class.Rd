%% ?MSnID

% showMethods(class="MSnID", where="package:MSnID")

% read_mzIDs
% assess_missed_cleavages
% assess_termini
%% get_psm_fdr
%% get_peptide_fdr
%% get_accession_fdr
%% get_peptides
%% get_accessions
% apply_filter
% id_quality
% evaluate_filter
% correct_peak_selection
% mass_measurement_error
% recalibrate
% names
% dim
% show
% psms
% psms<-
% as MSnSet
% $
% $<-
% [[
% [[<-
% optimize_filter


\name{MSnID-class}
\docType{class}
\alias{MSnID-class}
\alias{class:MSnID}
\alias{MSnID}



%% methods
\alias{dim,MSnID-method}
\alias{names,MSnID-method}
\alias{show,MSnID-method}
\alias{$,MSnID-method}
\alias{$<-,MSnID-method}
\alias{[[,MSnID-method}
\alias{[[,MSnID,ANY,ANY-method}
\alias{[[<-,MSnID-method}
\alias{assess_missed_cleavages,MSnID-method}
\alias{assess_termini,MSnID-method}
\alias{correct_peak_selection,MSnID-method}
\alias{evaluate_filter,MSnID-method}
\alias{get_peptides,MSnID-method}
\alias{get_accessions,MSnID-method}
\alias{id_quality,MSnID-method}
\alias{mass_measurement_error,MSnID-method}
\alias{read_mzIDs,MSnID-method}
\alias{recalibrate,MSnID-method}
\alias{psms,MSnID-method}
\alias{psms<-,MSnID,data.frame-method}

% coercion method
\alias{coerce,MSnID,MSnSet-method}
\alias{coerce,MSnID,data.table-method}






\title{
    The "MSnID" Class for Mass Spectrometry Based 
    Proteomics Identification Data
}

\description{
    The \code{MSnID} is a convenience class for manipulating 
    the MS/MS search results.
}

\section{Objects from the Class}{
    The way to create objects is to call \code{MSnID}
    constructor function that takes as an input the project working directory
    \code{workDir} and the second argument if the cache from previous analysis 
    should be cleaned \code{cleanCache}.
}

\section{Slots}{
    \describe{

    \item{\code{workDir}:}{Object of class \code{"character"}. 
        containing working directory for the project. 
        The \code{.Rcache} subdirectory stores the cached 
        resuls form the previous analyses.
        The mechanism of caching relies on \link[R.cache]{R.cache} package.} 

    \item{\code{psms}:}{Object of class \code{data.table} that contains all 
        the MS/MS identification results in the 
        form of peptide(or protein)-spectrum-matches.}

    }
}

\section{Methods}{
    \describe{ 

    \item{read_mzIDs}{\code{signature(.Object, mzids)}: \cr 
        Reads mzIdentML files into \code{psms} \code{data.table} slot 
        of \code{.Object} MSnID instance. The functionality leverage
        \link[mzID]{mzID} package facility. Note, the calls are memoised
        using \link[R.cache]{R.cache} facility. So if the call with the same
        list of files issues again, the results will be read from cache 
        instead of re-parsing the mzIdentML files.}

    \item{\code{psms(object)}, \code{psms(object)<-value}:}{
        Gets and sets MS/MS search results as \code{data.frame}}

    \item{dim}{\code{signature(x = "MSnID")}: \cr Returns the dimensions of 
        the table with MS/MS identification data.}

    \item{get_peptides}{\code{signature(.Object = "MSnID")}: \cr 
        Returns unique peptide list.}

    \item{get_accessions}{\code{signature(.Object = "MSnID")}: \cr 
        Returns unique accessions (typically proteins) list.}

    \item{assess_termini}{Checks the agreement of peptide termini 
        with enzymes cleavage specificity. The return value is 
        the\code{MSnID} object with extra variable \code{numIrregCleavages}.}

    \item{assess_missed_cleavages}{Checks if the peptide sequence contains the
        sites that were not cleaved by the enzyme. 
        For details see \code{\link[MSnID]{assess_missed_cleavages}}}

    \item{mass_measurement_error}{Returns parent ion mass measurement error in 
        parts per million (ppm) units. Note, it requires 
        \code{experimentalMassToCharge} and \code{calculatedMassToCharge} 
        variables to be set.}

    \item{recalibrate}{Recalibrates, that is removes systematic 
        error from \code{experimentalMassToCharge} measurements.}

    \item{correct_peak_selection}{Subtracts or adds the mass difference between
        C13 and C12 isotopes (1.0033548378 Da) if that reduces the mass error. 
        Such a mass error arises from the fact that instrument may peak 
        non-monoisotopic peak for fragmentation and thus report the mass 
        that is different by ~ 1 Da.}

    \item{apply_filter}{
        \code{signature(.Object="MSnID", .Filter="character")}\cr
        \code{signature(.Object="MSnID", .Filter="MSnIDFilter")}\cr
        The \code{.Filter} argument is a \code{"character"} or
        converted to a \code{"character"} text string that is
        evaluated to a \code{"logical"} for each entry of the
        MS/MS results table. Return value is a filtered
        \code{MSnID} object with entries that pass the applied filter.}

    \item{evaluate_filter}{
        \code{evaluate_filter(.Object, filter, 
        level = c("PSM","Peptide", "Accession")}\cr
        Returns a list with \code{fdr} and \code{n} elements. 
        Argument \code{filter} is either \code{"character"} 
        or \code{"MSnIDFilter"} object. Argument \code{level} can 
        take one of the values \code{c("PSM","Peptide","Accession")} 
        and controls the level filter is evaluated.}

    \item{\code{id_quality}}{
        \code{signature(.Object="MSnID", ...)}\cr
        Other optional \dots arguments are \code{filter} is an "MSnIDFilter"
        instance and \code{level}. The \code{level} values are one of "PSM",
        "Peptide", "Accession". The method returns FDR for given 
        level depending of type of identifications.}

    \item{as(,"MSnSet")}{
        \code{signature(x = "MSnID")}: 
        Coerce object from \code{MSnID} to \code{\link{MSnSet}}.}

    \item{names}{
        signature(x="MSnID")\cr
        Returns the column names in the MS/MS results table.}

    \item{\code{object$name}, \code{object$name<-value}}{
        Access and set \code{name} column in MS/MS search results table.}

    \item{\code{object[[i]]}, \code{object[[i]]<-value}}{
        Access and set column \code{i} (character or numeric index) in
        MS/MS search results table.}

    \item{as(,"MSnSet")}{
        \code{signature(from = "MSnID")}: 
        Coerce object from \code{MSnID} to \code{MSnSet}.}    

    \item{as(,"data.table")}{
        \code{signature(from = "MSnID")}: 
        Coerce object from \code{MSnID} to \code{data.table}.}    

    }
}

\author{
    Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}

\seealso{
    \code{\link[MSnbase]{MSnSet}},
    \code{\link[mzID]{mzID}}. 
}

\examples{
\dontrun{
msnidObj <- MSnID(".")
mzids <- system.file("extdata","c_elegans.mzid.gz",package="MSnID")
msnidObj <- read_mzIDs(msnidObj, mzids)
# clean up the cache directory
unlink(".Rcache", recursive=TRUE)
}
}

\keyword{classes}

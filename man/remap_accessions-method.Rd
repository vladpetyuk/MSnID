\name{remap_accessions}
\alias{remap_accessions}
\alias{remap_accessions,MSnID-method}

\title{Changes accessions from one protein id to another}

\description{
    Changes accessions from one protein id to another.
}

\usage{
    remap_accessions(object, 
                     conversion_table, 
                     extraction_pttrn=c("\\\\|([^|-]+)(-\\\\d+)?\\\\|",
                                        "([A-Z]P_\\\\d+)",
                                        "(ENS[A-Z0-9]+)"), 
                     path_to_FASTA=NULL)
}

\arguments{
    \item{object}{
        An instance of class "MSnID".
    }
    \item{conversion_table}{
        (data.frame) first column in the data frame corresponds to
        identifiers in the FASTA file. Second column is the new identifier.
    }
    \item{extraction_pttrn}{
        (string) regex pattern that extract protein identifier from
        FASTA entry name as first group (that is \code{"\\\\1"}).
        The most common patterns are the one corresponding to
        UniProt \code{"\\\\|([^|-]+)(-\\\\d+)?\\\\|"}, 
        RefSeq \code{"^([A-Z]P_\\\\d+)"} and
        ENSEMBL \code{"^(ENS[A-Z0-9]+)"}.
        Other regex patterns can be accepted as well.
        Defaul is UniProt pattern.
    }
    \item{path_to_FASTA}{
        (string) path to FASTA file. If provided only accessions present in
        the given FASTA file will be retained.
    }
}

\value{
    Returns an instance of "MSnID" with updated \code{accessions}.
}

\author{
    Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}

\seealso{
    \code{\link{remap_fasta_entry_names}}
}


\examples{
m <- MSnID(".")
mzids <- system.file("extdata","phospho.mzid.gz",package="MSnID")
m <- read_mzIDs(m, mzids)
head(m$accessions)
conv_tab <- fetch_conversion_table("Homo sapiens", "UNIPROT", "SYMBOL")
m2 <- remap_accessions(m, conv_tab, "\\\\|([^|-]+)(-\\\\d+)?\\\\|")
head(m2$accessions)
unlink(".Rcache", recursive=TRUE)
}



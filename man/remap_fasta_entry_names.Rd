\name{remap_fasta_entry_names}
\alias{remap_fasta_entry_names}
\title{Remapping entries in FASTA file}

\description{
Remaps entries in the FASTA file from one protein identifier to another
according to provided conersion table. Input is a path to FASTA file. Output
is also a path to a new FASTA file with updated entry names.
}

\usage{
remap_fasta_entry_names(path_to_FASTA,
                        conversion_table,
                        extraction_pttrn=c("\\\\|([^|-]+)(-\\\\d+)?\\\\|",
                                                 "([A-Z]P_\\\\d+)",
                                                 "(ENS[A-Z0-9]+)"))
}

\arguments{
    \item{path_to_FASTA}{
        (string) path to FASTA file
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
}

\value{
    path to new FASTA file
}

\author{
    Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}

\examples{
library(Biostrings)
fst_path <- system.file("extdata","for_phospho.fasta.gz",package="MSnID")
readAAStringSet(fst_path)
conv_tab <- fetch_conversion_table("Homo sapiens", "UNIPROT", "SYMBOL")
fst_path_2 <- remap_fasta_entry_names(fst_path, conv_tab, "\\\\|([^|-]+)(-\\\\d+)?\\\\|")
readAAStringSet(fst_path_2)
file.remove(fst_path_2)
}




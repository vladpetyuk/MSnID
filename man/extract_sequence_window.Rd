\name{extract_sequence_window}
\alias{extract_sequence_window}
\alias{extract_sequence_window,MSnID-method}
\title{Extracts +/-7 amino acids to the PTM}

\description{
Given the peptide sequence with modification X.XXXX*XXXX.X and provided
protein sequence FASTA, the method maps the 7 amino acids to the left and right of each PTM, with "-" character appended if PTM is near start or end of protein sequence. The modified AA is turned into lowercase. 
}

\usage{
    extract_sequence_window(object,
                           fasta,
                           accession_col="accession",
                           site_loc_col="SiteLoc",
                           radius=7L,
                           collapse="|")
}


\arguments{
    \item{object}{
        An instance of class MSnID.
    }
    \item{fasta}{
        (AAStringSet object) Protein sequences read from a FASTA file.
        Names must match protein/accesison IDs in the accesson column
        of the MSnID object.
    }
    \item{accession_col}{
        (character) Name of accession column.
    }
    \item{site_loc_col}{
        (character) Name of column containing site locations.
    }
    \item{radius}{
        (integer) How many amino acids to map to left and right of the PTM.
    }
    \item{collapse}{
        (character) Symbol to separate multiple PTMs.
    }
}

\value{
    MSnID object with an extra column \code{sequenceWindow} giving the neighborhood of each PTM.
}

\author{
    Michael Nestor \email{michael.nestor@pnnl.gov}
}


\examples{
m <- MSnID(".")
mzids <- system.file("extdata","phospho.mzid.gz",package="MSnID")
m <- read_mzIDs(m, mzids)

# to know the present mod masses
report_mods(m)

# TMT modification
m <- add_mod_symbol(m, mod_mass="229.1629", symbol="#")
# alkylation
m <- add_mod_symbol(m, mod_mass="57.021463735", symbol="^")
# phosphorylation
m <- add_mod_symbol(m, mod_mass="79.966330925", symbol="*")

# show the mapping
head(unique(subset(psms(m), select=c("modification", "peptide_mod"))))

# read fasta for mapping modifications
fst_path <- system.file("extdata","for_phospho.fasta.gz",package="MSnID")
library(Biostrings)
fst <- readAAStringSet(fst_path)
# to ensure names are the same as in accessions(m)
names(fst) <- sub("(^[^ ]*) .*$", "\\1", names(fst))
# # mapping phosphosites
m <- map_mod_sites(m, fst, "accession", "peptide_mod", "*", "lower")

# # mapping +/-7 amino acids to the PTm
m <- extract_sequence_window(m, fst)

head(unique(subset(psms(m), select=c("accession", "peptide", "sequenceWindow"))))

# clean-up cache
unlink(".Rcache", recursive=TRUE)
}


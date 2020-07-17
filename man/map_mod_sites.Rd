\name{map_mod_sites}
\alias{map_mod_sites}
\alias{map_mod_sites,MSnID-method}
\title{Maps the modifications to protein sequence}

\description{
Given the peptide sequence with modification X.XXXX*XXXX.X and provided
protein sequence FASTA, the method maps the location of the 
modification resulting in \{protein ID\}-\{aa\}\{aa position\}".
}

\usage{
    map_mod_sites(object,
                  fasta, 
                  accession_col = "accession", 
                  peptide_mod_col = "peptide_mod", 
                  mod_char = "*",
                  site_delimiter = "lower")
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
        (string) Name of the column with accession/protein IDs in the 
        MSnID object. Default is "accession".
    }
    \item{peptide_mod_col}{
        (string) Name of the column with modified peptide sequences in the 
        MSnID object. Default is "peptide_mod".
    }
    \item{mod_char}{
        (string) character that annotates the position of the modification.
        Default is "*".
    }
    \item{site_delimiter}{
        (string) either a single character or "lower" (default) meaning it will
        be the same amino acid symbol, but in lower case
    }
}

\value{
    MSnID object with extra columns regarting the modification mapping.
    Most likely, what you need is \code{SiteID}.
    \item{PepLoc}{(list of ints) 
        position of the starting amino acid within protein sequence.
        It is a list, because there may be multiple occurences of the
        same sequence matching the peptide's sequence.}
    \item{PepLocFirst}{(int)
        position of the first occurence of the matching sequence}
    \item{ProtLength}{(int)
        protein length}
    \item{ModShift}{(vector of ints)
        positions of modified amino acids within peptide}
    \item{ModAAs}{(vector of characters)
        single-letter amino acid codes of the modified residues}
    \item{SiteLoc}{(list of vectors of ints)
        positions of the modified amino acids within 
        protein for each occurence of the peptide}
    \item{Site}{(list of vectors of characters)
        modified sites encoded as amino acid symbol follwed by position
        for each occurence of the peptide}
    \item{SiteCollapsed}{(list of characters)
        same as \code{Site}, but modified sites for each peptide occurence
        collapsed into one string}
    \item{SiteCollapsedFirst}{(character)
        first element of the \code{SiteCollapsed} list
    }
    \item{SiteID}{(character)
        accession ID concatenated with \code{SiteCollapsedFirst}}
}

\author{
    Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
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

head(unique(subset(psms(m), select=c("accession", "peptide_mod", "SiteID"))))

# clean-up cache
unlink(".Rcache", recursive=TRUE)
}


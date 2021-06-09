\name{map_peptide_position}
\alias{map_peptide_position}
\alias{map_peptide_position,MSnID-method}
\title{Maps the peptide location(s) to protein sequence}

\description{
The method maps the location of the 
modification resulting in \{protein ID\}-\{aa\}\{aa position\}.
}

\usage{
    map_peptide_position(object,
                        fasta, 
                        accession_col = "accession")
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
}

\value{
    MSnID object with extra columns regarting the peptide location.
    \item{cleanSeq}{(character)
        peptide sequence without modification characters and flanking
        amino acids}
    \item{First_AA}{(list of ints)
        position of the starting amino acid within protein sequence.
        It is a list, because there may be multiple occurences of the
        same sequence matching the peptide's sequence.}
    \item{Last_AA}{(list of ints)
        Same as First_AA, but last amino acid}
    \item{First_AA_First}{(int)
        First element of First_AA. 
        Essentially, a first occurence of the peptide in the sequence.}
    \item{Last_AA_First}{(int)
        First element of Last_AA.}
}

\author{
    Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}


\examples{
m <- MSnID(".")
mzids <- system.file("extdata","phospho.mzid.gz",package="MSnID")
m <- read_mzIDs(m, mzids)
fst_path <- system.file("extdata","for_phospho.fasta.gz",package="MSnID")
library(Biostrings)
fst <- readAAStringSet(fst_path)
# ensure names are the same format as accessions(m)
names(fst) <- sub("(^[^ ]*) .*$", "\\1", names(fst))

# Creating sequences with repeats. This is just for the sake of 
# demonstration of the capability.
for(i in 2:4){
    fst[i] <- paste0(as.character(fst[i]),as.character(fst[i]))
}

# Appending reverse hits. Also, just for the demonstrating the capability.
mr <- mf <- apply_filter(m, "!isDecoy")
library(stringi)
mr$peptide <- stringi::stri_reverse(mr$peptide)
mr$accession <- paste0("XXX_", mr$accession)
mr$isDecoy <- TRUE
psms(m) <- rbind(psms(mr), psms(mf))

# the main call
m2 <- map_peptide_position(m, fst)

head(unique(subset(psms(m2), select=c("accession", 
                                      "peptide", 
                                      "First_AA",
                                      "First_AA_First",
                                      "Last_AA",
                                      "Last_AA_First",
                                      "ProtLen"))))

# clean-up cache
unlink(".Rcache", recursive=TRUE)
}


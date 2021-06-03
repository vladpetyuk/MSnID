\name{plot_protein_coverage}
\alias{plot_protein_coverage}
\alias{plot_protein_coverage,MSnID-method}
\title{Shows the peptide location(s) on the protein sequence}

\description{
Shows the peptide location(s) on the protein sequence.
}

\usage{
    plot_protein_coverage(object, 
                          accession, 
                          save_plot = FALSE)
}


\arguments{
    \item{object}{
        An instance of class MSnID.
    }
    \item{accession}{
        (string) accession name.
    }
    \item{save_plot}{
        (boolean) if TRUE, the plot is saved as png file.
    }
}

\value{
    if save_plot = FALSE, then returns ggplot2 object
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
plot_protein_coverage(m2, "sp|Q9Y2W1|TR150_HUMAN", save_plot = FALSE)

# clean-up cache
unlink(".Rcache", recursive=TRUE)
}


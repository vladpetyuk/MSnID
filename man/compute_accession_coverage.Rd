\name{compute_accession_coverage}
\alias{compute_accession_coverage}
\alias{compute_accession_coverage,MSnID-method}
\title{Compute accession coverage from identified peptides}

\description{
Given a list of accession and peptide pairs, computes the proportion of amino acid sequence covered by the listed peptides for each accession.
}

\usage{
    compute_accession_coverage(object,
                               fasta,
                               accession_col="accession",
                               pepSeq_col="pepSeq")
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
        (character) Name of the accession column.
    }
    \item{pepSeq_col}{
        (character) Name of the peptide sequence column.
    }
}

\value{
    MSnID object with an extra column \code{percentAAcoverage} with values in 0-100.
}

\author{
    Michael Nestor \email{michael.nestor@pnnl.gov}
}


\examples{
m <- MSnID(".")
mzids <- system.file("extdata","phospho.mzid.gz",package="MSnID")
m <- read_mzIDs(m, mzids)

# read fasta
fst_path <- system.file("extdata","for_phospho.fasta.gz",package="MSnID")
library(Biostrings)
fst <- readAAStringSet(fst_path)

# compute accession coverage
m <- compute_accession_coverage(m, fst)

head(unique(subset(psms(m), select=c("accession", "percentAACoverage"))))

# clean-up cache
unlink(".Rcache", recursive=TRUE)
}


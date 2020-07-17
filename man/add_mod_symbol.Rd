\name{add_mod_symbol}
\alias{add_mod_symbol}
\alias{add_mod_symbol,MSnID-method}
\title{Annotates peptide sequences with modification symbols}

\description{
Given the provided modification mass, annotates its position within the peptide
sequence with provided symbol.
}

\usage{
    add_mod_symbol(object, mod_mass, symbol)
}


\arguments{
    \item{object}{
        An instance of class "MSnID".
    }
    \item{mod_mass}{
        (string) modification mass. Must match exactly one of the masses in
        \code{report_mods}.
    }
    \item{symbol}{
        (string) character that annotates the position of the modification
    }    
}

\note{
    In current implementation the method can not distinguish modifications
    with the same mass, but different amino acid specificity as 
    different modifications.
}

\value{
    Returns MSnID object with new or modified \code{peptide_mod} column.
    The column contains peptide sequence with amino acid 
    modifications annotated with symbols.
}

\seealso{
    \code{\link{report_mods}}
    \code{\link{map_mod_sites}}
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

# clean-up cache
unlink(".Rcache", recursive=TRUE)
}


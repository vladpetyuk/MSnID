\name{fetch_conversion_table}
\alias{fetch_conversion_table}
\title{Fetches conversion table form one type of identifiers to another}

\description{
A wrapper function over \pkg{AnnotationHub} that helps to convert
from one protein identifiers to another.
}

\usage{
    fetch_conversion_table(organism_name, 
                           from, to, 
                           backend="AnnotationHub", 
                           snapshot_date=NULL)
}

\arguments{
    \item{organism_name}{
        (string) official organism name. 
        E.g. "Homo sapiens", "Mus musculus" or "Rattus norvegicus".
    }
    \item{from, to}{
        (string) identifier names. Recommended names are
        "SYMBOL", "UNIPROT", "REFSEQ" and "ENSEMBLPROT". Other identifiers
        are possible, but use them at your own risk.
    }
    \item{backend}{
        (string) currently only \pkg{AnnotationHub}
    }
    \item{snapshot_date}{
        (string) snapshot date for \pkg{AnnotationHub}. 
        Default is \code{NULL} meaning latest date.
    }
}

\value{
    data.frame with first column name of from identifier, the second name
    of to identifier
}

\author{
    Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}

\examples{
conv_tbl <- fetch_conversion_table("Rattus norvegicus", "REFSEQ", "SYMBOL")
head(conv_tbl)
conv_tab <- fetch_conversion_table("Homo sapiens", "UNIPROT", "SYMBOL")
head(conv_tbl)
}




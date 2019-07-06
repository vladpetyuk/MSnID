\name{assess_missed_cleavages}
\alias{assess_missed_cleavages}

\title{Counts the missing cleavage sites within the peptides sequence}

\description{
    Bottom-up proteomics approaches utilize endoproteases or chemical agents
    to digest proteins into smaller fragments called peptides. 
    The enzymes recognize short amino acid motifs
    and cleave along the peptide bonds. Chemical agents such as 
    \code{CNBr} also possess amino acid cleavage specificity. In real-world
    not every cleavage site get cleaved during the sample processing. 
    Therefore settings of MS/MS search engines quite often explictly 
    allow up to a certain number missed clevage sites per peptide sequence.

    This function counts the number of missed cleavages in peptide
    sequence given the endoprotease cleavage motif in the form of
    regular expression. The default value for \code{missedCleavagePattern} is
    \code{[KR](?=[^P$])}, which corresponds to trypsin.
}

\usage{
    assess_missed_cleavages(object, missedCleavagePattern="[KR](?=[^P$])")
}

%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{An instance of class "MSnID".}

    \item{missedCleavagePattern}{
        Cleavage pattern in the form of regular expression.
    }
}

\value{
    Returns an instance of "MSnID" class with 
    additional column "numMissCleavages"
}

\author{
    Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}

\section{Warning}{
    If the "MSnID" instance does not contain "peptide" column in MS/MS results
    table then there will be an error.
    E.g. you can check this by \cr \code{"peptide" \%in\% names(msnid)} 
    where \code{msnid} is your "MSnID" instance.
}

\seealso{
    \code{\link{assess_termini}}
}


\examples{
data(c_elegans)
# adding column numMissCleavages containing count of missed cleavages
msnidObj <- assess_missed_cleavages(msnidObj, 
                                    missedCleavagePattern="[KR](?=[^P$])")
# check the distribution
table(msnidObj$numMissCleavages)
}



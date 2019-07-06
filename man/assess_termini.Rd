\name{assess_termini}
\alias{assess_termini}

\title{Checks if the peptide termini conforms with cleavage specificity}

\description{
    Bottom-up proteomics approaches utilize endoproteases or chemical agents
    to digest proteins into smaller fragments called peptides. 
    The enzymes recognize short amino acid motifs
    and cleave along the peptide bonds. Chemical agents such as \code{CNBr} also
    possesses amino acid cleavage specificity. 

    This function checks if peptide termini are as expected given the
    enzymatic/chemical cleavage specificity. 
    The default value for \code{validCleavagePattern} is
    \code{[KR]\\.[^P]}, which corresponds to trypsin.
}

\usage{
    assess_termini(object, validCleavagePattern="[KR]\\\\.[^P]")
}

%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{An instance of class "MSnID".}
    \item{validCleavagePattern}{
        Cleavage pattern in the form of regular expression.
    }
}

\details{
    N- or C- protein termini are not considered as irregular clevages sites.
}

\value{
    Returns an instance of "MSnID" class with additional 
    column "numIrregCleavages". If both termini conforms with 
    cleavage specificity, then value is \code{0}, if one or two 
    termini are irregular then the values are \code{1} and \code{2},
    correspondingly.
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
    \code{\link{assess_missed_cleavages}}
}


\examples{
data(c_elegans)
# adding column numIrregCleavages 
# containing count of irregularly cleaved termini
msnidObj <- assess_termini(msnidObj, validCleavagePattern="[KR]\\\\.[^P]")
# check the distribution
table(msnidObj$numIrregCleavages)
}

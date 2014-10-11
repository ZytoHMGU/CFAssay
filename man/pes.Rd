\name{pes}
\alias{pes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of plating efficiencies from a curve data set containing one specified curve
}
\description{
The function calculates plating efficiencies, i.e. fractions of colonies per cell of untreated cells, for every experimental replicate in a data frame with one specified curve. For that, lines with zero dose (\eqn{dose = 0}) are extracted from the data frame.
}
\usage{
pes(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A data frame which contains columns \code{Exp, dose, ncells, ncolonies}.
}
}
\details{
In the data frame \code{X}, \code{Exp} identifies the experimental replicates and may be numeric or non-numeric. 
}
\value{
The function returns a data frame with two columns \code{pe} and \code{S0}. The first column contains the measured plating efficiencies the second column contains plating efficiencies fitted separately for each repeated experiments. Rows of the data frame are named for the different experiments.
}
\author{
Herbert Braselmann
}
\examples{
datatab<- read.table(system.file("doc", "expl1_cellsurvcurves.txt", package="CFAssay"), header=TRUE, sep="\t")
X<- subset(datatab, cline=="okf6TERT1")  #Specification of curve
pes(X)
pes(subset(datatab, cline=="cal33") )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
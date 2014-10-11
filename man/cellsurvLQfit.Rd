\name{cellsurvLQfit}
\alias{cellsurvLQfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit the linear-quadratic (LQ) model to cell survival data
}
\description{
This function calculates the linear coefficient \eqn{alpha} and the coefficient \eqn{beta} of the dose-squared term (see manual for this R-package) for colony counts measured for a set of irradiation doses and repeated experiments. The function is a wrapper for the R-functions \code{glm} or \code{lm}, which simplifies use of these functions for cell survival data.
}
\usage{
cellsurvLQfit(X, method="ml", PEmethod="fit")
}
\arguments{
  \item{X}{
A data frame which contains at least columns \code{Exp, dose, ncells, ncolonies} and if there is no 0-value in the \code{dose}-column, \code{X} has to contain a further column \code{pe} for plating efficiencies.
}
  \item{method}{
Determines the method used for the fit. \code{"ml"} is for maximum-likelihood, \code{"ls"} for least-squares. \code{"franken"} performs weigthed least-squares with weights as described in Franken et al. (2006).
}
  \item{PEmethod}{
Controls the value of the plating efficiencies, i.e. the colony counts for untreated cells. \code{"fit"} calculates fitted plating efficiencies as model parameters, \code{"fix"} uses fixed ones calculated from the observed zero dose data or from a column named \code{pe} in \code{X}.
}
}
\details{
In the data frame \code{X}, \code{Exp} identifies the experimental replicates and may be numeric or non-numeric. \code{method="ml"} uses R function \code{\link{glm}} with quasipoisson family and link function \code{"log"}. \code{method="ls"} uses R function \code{\link{lm}}. \code{PEmethod="fit"} fits plating efficiencies for every experiments. \code{PEmethod="fix"} uses observed plating efficiencies. If there is no 0-value in the \code{dose}-column, \code{PEmethod} is overwritten with \code{"fix"} and \code{X} has to contain a further column \code{pe} containing the plating efficiencies, i.e. \code{ncolonies/ncells} from untreated cells, not per hundred or percent.
}
\value{
The function returns an object of class \code{cellsurvLQfit}, which is similar to classes \code{glm} or \code{lm}, however containing two additional entries, \code{type} and \code{PEmethod}, which are used for printing and plotting. The full result is returned invisibly, i.e. the function has to be used with \code{print} or \code{plot} or assigned to a variable, say for e.g. \code{fit} as in the example below.
}
\references{
Franken NAP, Rodermond HM, Stap J, et al. Clonogenic assay of cells in vitro. Nature Protoc 2006;1:2315-19.
}
\author{
Herbert Braselmann
}
\seealso{
\code{\link{glm}} and \code{\link{family}} with references for generalized linear modelling, \code{lm}
}
\examples{
datatab<- read.table(system.file("doc", "expl1_cellsurvcurves.txt", package="CFAssay"), header=TRUE, sep="\t")
X<- subset(datatab, cline=="okf6TERT1")
fit<- cellsurvLQfit(X) #using default options
print(fit)
print(fit$type)
print(fit$PEmethod)
#Using other options:
print(cellsurvLQfit(X, method="ls"))
print(cellsurvLQfit(X, PEmethod="fix"))
print(cellsurvLQfit(X, method="ls", PEmethod="fix"))
print(cellsurvLQfit(X, method="franken"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

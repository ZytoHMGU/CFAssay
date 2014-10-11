\name{cellsurvLQdiff}
\alias{cellsurvLQdiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Comparison of two linear-quadratic cell survival curves
}
\description{
The function does an ANOVA test for overall comparison of the parameters \eqn{alpha} and \eqn{beta} of two linear-quadratic cell survival curves. The parameters are fitted simultaneously to the data with this function, i.e. no other function is necessary to derive the fits.
}
\usage{
cellsurvLQdiff(X, curvevar, method="ml", PEmethod="fit")
}
\arguments{
  \item{X}{
A data frame which contains columns \code{Exp, dose, ncells, ncolonies} and a further column containing two different values (character strings), which identify the two curves. Moreover, if there is no 0-value in the \code{dose}-column, \code{X} has to contain a column \code{pe} for plating efficiencies.
}
  \item{curvevar}{
Character string, which has to be one of the column names of the data frame \code{X}, that contains the two different values (character strings that distinguishes between the two curves).
}
  \item{method}{
Determines the method used for the fit. \code{"ml"} is for maximum-likelihood, \code{"ls"} for least-squares. \code{"franken"} performs weigthed least-squares with weights as described in Franken et al. (2006).
}
  \item{PEmethod}{
Controls the value of the plating efficiencies, i.e. the colony counts for untreated cells. \code{"fit"} calculates fitted plating efficiencies as model parameters, \code{"fix"} uses fixed ones calculated from the observed zero dose data or from a column named \code{pe} in \code{X}.
}
}
\details{
In the data frame \code{X}, \code{Exp} identifies the experimental replicates and may be numeric or non-numeric. \code{method="ml"} for maximum-likelihood uses R function \code{\link{glm}} with \code{family} "quasipoisson" and link function \code{"log"}. \code{method="ls"} uses R function \code{\link{lm}}.
}
\value{
The function returns an object of class \code{cellsurvLQdiff} containing three elements, \code{fit1}, \code{fit2} and \code{anv}. \code{fit1} and \code{fit2} are objects of class \code{glm} when \code{method="ml"} or of class \code{lm} when \code{method="ls"}. \code{fit1} has parameters alpha and beta fitted in common for both cell survival curves.  \code{fit2} has parameters alpha and beta fitted differently for both curves. \code{anv} is of class \code{anova} and contains the F-test. Test results are printed, however, the full result inlcuding curve parameters is returned invisibly, i.e. the function has to be used with \code{print} or assigned to a variable, say for e.g. \code{fitcomp} as in the example below.
}
\author{Herbert Braselmann}
\seealso{
\code{\link{glm}} and \code{\link{family}} with references for generalized linear modelling. \code{\link{anova}}, \code{\link{cellsurvLQfit}}.
}
\examples{
datatab<- read.table(system.file("doc", "expl1_cellsurvcurves.txt", package="CFAssay"), header=TRUE, sep="\t")
names(datatab) #contains a column "cline"
table(datatab$cline)
fitcomp<- cellsurvLQdiff(datatab, curvevar="cline") #using default options
print(fitcomp)
plot(cellsurvLQfit(subset(datatab, cline=="okf6TERT1")), col=1)
plot(cellsurvLQfit(subset(datatab, cline=="cal33")), col=2, add=TRUE)
legend(0, 0.02, c("okf6TERT1", "cal33"), text.col=1:2)
#using different options:
print(cellsurvLQdiff(datatab, curvevar="cline", method="ls"))
print(cellsurvLQdiff(datatab, curvevar="cline", PEmethod="fix"))
print(cellsurvLQdiff(datatab, curvevar="cline", method="ls", PEmethod="fix"))
print(cellsurvLQdiff(datatab, curvevar="cline", method="franken"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

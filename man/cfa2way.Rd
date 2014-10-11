\name{cfa2way}
\alias{cfa2way}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Two-way analysis for cell survival data
}
\description{
The function does an ANOVA of cell survival data from experimental 2-way designs where a treatment factor is tested on a control and on an altered cell line or where two different simultaneous treatments are tested on cells from a common unaltered clone. The function is a wrapper for the R-function \code{glm}. quasipoisson family is used with link function \code{"log"}, i.e. dependency of treatment factors is considered as logarithmically additive.
}
\usage{
cfa2way(X, A, B, param="A/B", method="ml")
}
\arguments{
  \item{X}{
a data frame which contains columns \code{Exp, ncells, ncolonies} and two further columns for the treatment variables, see details.
}
  \item{A}{
a character string containing the name of a treatment or cell line variable (first factor in the model)
}
  \item{B}{
a character string containing the name of a treatment or cell line variable (second factor in the model)
}
  \item{param}{
Controls the parametrization of the model. Options are "A/B" for B nested in A , "B/A" for A nested in B and "A*B" for interaction term.
}
  \item{method}{
determines the method used for the fit. \code{"ml"} is for maximum-likelihood, \code{"ls"} for least-squares.
}
}
\details{
In the data frame \code{X}, \code{Exp} identifies the experimental replicates and may be numeric or non-numeric. The two treatment or cell line columns should have numeric values 0, 1, ... for 2, 3, ... levels. For e.g. if a column describes clonal alteration (transfection, knock-down etc.) by a gene then 0 means unaltered or control and 1 means altered. Similar if a column describes treatment with one dose then 0 means untreated and 1 treated. 2 would indicate another dose level from the same treatment drug without taking it as a continuous covariate as for cell survival curves for radiation.
}
\value{
The function returns an object of class \code{cfa2way} containing three elements, \code{fit1}, \code{fit2} and \code{anv}. \code{fit1} and \code{fit2} are objects of class \code{glm} when \code{method="ml"} or of class \code{lm} when \code{method="ls"}. \code{fit1} has logarithmic additive parameters without interaction.  \code{fit2} has logarithmic additive parameters and interaction. \code{anv} is of class \code{\link{anova}} and contains the F-test. The full result is returned invisibly, i.e. the function has to be used with \code{print} or assigned to a variable, say for e.g. \code{fitcomp} as in the example below.
}
\author{Herbert Braselmann}
\seealso{
\code{\link{glm}} and \code{\link{family}} with references for generalized linear modelling.
}
\examples{
datatab<- read.table(system.file("doc", "exp2_2waycfa.txt", package="CFAssay"), header=TRUE, sep="\t")
names(datatab)  # has columns "x5fuCis" and "siRNA" 
fitcomp<- cfa2way(datatab, A="siRNA", B="x5fuCis", param="A/B")
print(fitcomp, labels=c(A="siRNA",B="x5fuCis"))
print(cfa2way(datatab, A="siRNA", B="x5fuCis", param="A/B", method="ls"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

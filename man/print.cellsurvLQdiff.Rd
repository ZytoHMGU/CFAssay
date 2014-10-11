\name{print.cellsurvLQdiff}
\alias{print.cellsurvLQdiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print test results of comparison of two linear-quadratic cell survival curves
}
\description{
The function prints linear coefficients alpha and dose-squared coefficients beta of calculated with function \code{cellsurvLQdiff}. In addition quality statistics of the fit and ANOVA F-test results for overall comparison of the coefficients are printed.
}
\usage{
\method{print}{cellsurvLQdiff}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{cellsurvLQdiff} resulting from function \code{cellsurvLQdiff}.
}
  \item{\dots}{
further arguments to pass to R function \code{print}.
}
}
\author{Herbert Braselmann}
\seealso{
\code{\link{cellsurvLQdiff}}
}
\examples{
datatab<- read.table(system.file("doc", "expl1_cellsurvcurves.txt", package="CFAssay"), header=TRUE, sep="\t")
fitcomp<- cellsurvLQdiff(datatab, curvevar="cline") #using default options
print(fitcomp)
#using different options for cellsurvLQdiff:
print(cellsurvLQdiff(datatab, curvevar="cline", method="ls"))
print(cellsurvLQdiff(datatab, curvevar="cline", PEmethod="fix"))
print(cellsurvLQdiff(datatab, curvevar="cline", method="ls", PEmethod="fix"))
print(cellsurvLQdiff(datatab, curvevar="cline", method="franken"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
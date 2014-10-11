\name{print.cellsurvLQfit}
\alias{print.cellsurvLQfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print summary of an LQ-model fit for cell survival data
}
\description{
The function prints the results of an LQ-model fit for radiation dose dependent cell survival.
}
\usage{
\method{print}{cellsurvLQfit}(x, ...)
}
\arguments{
  \item{x}{
an object of class \code{cellsurvLQfit} resulting from \code{\link{cellsurvLQfit}}.
}
  \item{\dots}{
further arguments to pass to R function \code{print}.
}
}
\details{
In this version of CFAssay the class argument \code{x} is checked by its entry \code{fit$type} and results in an error, when \code{x} results from an independent use of \code{\link{glm}} or \code{\link{lm}}. 
}
\author{Herbert Braselmann}
\seealso{
\code{\link{cellsurvLQfit}}
}
\examples{
datatab<- read.table(system.file("doc", "expl1_cellsurvcurves.txt", package="CFAssay"), header=TRUE, sep="\t")
X<- subset(datatab, cline=="okf6TERT1")
fit<- cellsurvLQfit(X) #using default options
print(fit)
print(fit$type)
print(fit$PEmethod)
#using other options
print(cellsurvLQfit(X, method="ls"))
print(cellsurvLQfit(X, PEmethod="fix"))
print(cellsurvLQfit(X, method="ls", PEmethod="fix"))
print(cellsurvLQfit(X, method="franken"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

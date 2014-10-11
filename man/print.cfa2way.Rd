\name{print.cfa2way}
\alias{print.cfa2way}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print summary of two-way analysis for cell survival data.
}
\description{
The function prints a summary of two-way analysis for cell survival data.
}
\usage{
\method{print}{cfa2way}(x, labels=c(A="A",B="B"), ...)
}
\arguments{
  \item{x}{
an R object of class \code{cfa2way} and resulting from function \code{\link{cfa2way}}.
}
  \item{labels}{
labels for output description. These describe the meaning of \code{A} and of \code{B}.
}
  \item{\dots}{
further arguments to pass to R function \code{print}.
}
}
\author{
Herbert Braselmann}
\seealso{
\code{\link{cfa2way}}
}
\examples{
datatab<- read.table(system.file("doc", "exp2_2waycfa.txt", package="CFAssay"), header=TRUE, sep="\t")
names(datatab)  # has columns "x5fuCis" and "siRNA" 
fitcomp<- cfa2way(datatab, A="siRNA", B="x5fuCis", param="A/B")
print(fitcomp, labels=c(A="siRNA",B="x5fuCis"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

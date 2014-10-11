\name{plotExp.cfa2way}
\alias{plotExp.cfa2way}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diagnostic plot of experimental replicates for two-way analysis of cell survival data
}
\description{
This function plots the estimated means of a two-way ANOVA for each experiment in a series of plots. It is recommended to direct it in a pdf-file.
}
\usage{
\method{plotExp}{cfa2way}(x, labels=c(A="A",B="B"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an R object of class \code{cfa2way} and resulting from function \code{\link{cfa2way}}.
}
  \item{labels}{
Labels for output description. These describe the meaning of \code{A} and of \code{B}.
}
  \item{\dots}{
further arguments to pass to generic function \code{plotExp}.
}
}
\author{Herbert Braselmann}
\seealso{
\code{\link{cfa2way}}
}
\examples{
datatab<- read.table(system.file("doc", "exp2_2waycfa.txt", package="CFAssay"), header=TRUE, sep="\t")
names(datatab)  # has columns "x5fuCis" and "siRNA" 
fitcomp<- cfa2way(datatab, A="siRNA", B="x5fuCis", param="A/B")
print(fitcomp, labels=c(A="siRNA",B="x5fuCis"))
  pdf("TwoWay_experimental_plots.pdf")
    plotExp(fitcomp, labels=c(A="siRNA", B="x5fuCis"))
  dev.off()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

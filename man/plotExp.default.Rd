\name{plotExp.default}
\alias{plotExp.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Default function for plotting of experimental repeats
}
\description{
plotExp.default is the formal default method of the generic \code{\link{plotExp}} function.
}
\usage{
\method{plotExp}{default}(x, ...)
}
\arguments{
  \item{x}{
should be an object of class \code{cellsurvLQfit} or \code{cfa2way}, see details.
}
  \item{\dots}{
further arguments to pass to generic function \code{plotExp}, dependent of the class of the argument \code{x}.
}
}
\details{
In this version of CFAssay other objects than of class \code{cellsurvLQfit} or \code{cfa2way} this will give an error message and a hint to use one of the two defined classes.
}
\author{Herbert Braselmann}
\seealso{
\code{\link{plotExp.cellsurvLQfit}} and \code{\link{plotExp.cfa2way}}
}
\examples{
datatab<- read.table(system.file("doc", "expl1_cellsurvcurves.txt", package="CFAssay"), header=TRUE, sep="\t")
X<- subset(datatab, cline=="okf6TERT1")
fit<- cellsurvLQfit(X)
#pdf("okf6TERT1_experimental_plots.pdf")
   plotExp(fit)
#dev.off()
\dontrun{
x <- 1
plotExp(X) #yields an error for this data type
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
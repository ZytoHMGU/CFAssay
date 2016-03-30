\name{sfpmean}
\alias{sfpmean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pointwise mean survival fractions for curves with several experimental replicates
}
\description{
The function calculates mean survival fractions for curves averaged over experimental replicates. The function is employed by function \code{plot.cellsurvLQfit for plotting observed means}
}
\usage{
sfpmean(X, S0=NULL)
}
\arguments{
  \item{X}{
A data frame which contains columns \code{Exp, dose, ncells, ncolonies} and if \code{S0=NULL}, \code{X} has to contain a further column \code{pe} for plating efficiencies.
}
  \item{S0}{
If not \code{NULL}, a named numerical vector of length equal to the number of different experiments, i.e. \code{length(S0)==length(unique(X$Exp))} has to be \code{TRUE}. Default is \code{S0=NULL}, i.e. undefined.
}
}
\details{
In the data frame \code{X}, \code{Exp} identifies the experimental replicates and may be numeric or non-numeric. \code{S0} may contain plating efficiencies for each replicate, resulting from function \code{\link{pes}} or from \code{\link{cellsurvLQfit}} (fitted). When \code{S0=NULL}, X must have a column with name \code{pe}, containing the plating efficiencies.
}
\value{
A numerical matrix with two rows, the first row containing the survival fractions for each radiation dose, second row the standard deviations.
}
\author{Herbert Braselmann}
\seealso{
\code{\link{pes}}, \code{\link{cellsurvLQfit}}, \code{\link{plot.cellsurvLQfit}}
}
\examples{
datatab <- read.table(system.file("doc", "expl1_cellsurvcurves.txt", package="CFAssay"), header=TRUE, sep="\t")
X <- subset(datatab, cline=="okf6TERT1")
S0 <- pes(X)$pe  #observed plating efficiencies
length(S0)==length(unique(X$Exp)) #length ok?
names(S0) <- pes(X)$Exp
sfpmean(X, S0)
fit <- cellsurvLQfit(X)
fit$coef  #contains fitted log-pe
grep("Exp",names(fit$coef))
S01 <- exp(fit$coef[1:8])  #fitted pe
sfpmean(X, S01)
\dontrun{sfpmean(X)}  #yields an error for this data set
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\name{plot.cellsurvLQfit}
\alias{plot.cellsurvLQfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of an LQ model fit
}
\description{
This function plots a cell survival curve derived from fitting an LQ model with function \code{cellsurvLQfit}
}
\usage{
\method{plot}{cellsurvLQfit}(x, xlim = NULL, ylim = c(0.008, 1), xlab = "Dose (Gy)", ylab = "Survival (1 = 100\%)",col=1, pch=1, add=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{cellsurvLQfit} resulting from function \code{\link{cellsurvLQfit}}.
}
  \item{xlim}{
plot range for the x-axis. Default is the dose range of the data.
}
  \item{ylim}{
plot range for the y-axis. Default is from 0.008 to 1.0
}
  \item{xlab}{
label for the x-axis. Default is "Dose (Gy)".
}
  \item{ylab}{
label for the y-axis. Default is "Survival (1 = 100\%)".
}
  \item{col}{
colour for plot. Default is \code{col = 1}.
}
  \item{pch}{
symbol for plotting points. Default is \code{pch = 1}.
}
  \item{add}{
logical; if TRUE add to an already existing plot, see \code{\link{curve}}.
}
  \item{\dots}{
further arguments to pass to R function \code{plot}.
}
}
\author{
Herbert Braselmann
}
\seealso{
\code{\link{cellsurvLQfit}}
}
\examples{
datatab<- read.table(system.file("doc", "expl1_cellsurvcurves.txt", package="CFAssay"), header=TRUE, sep="\t")
X<- subset(datatab, cline=="okf6TERT1")
fit<- cellsurvLQfit(X)
plot(fit)
S0 <- pes(X)$S0
names(S0) <- pes(X)$Exp
sfpmean(X, S0) #values of plotted mean survival fractions and error bars 
# add second plot
plot(cellsurvLQfit(subset(datatab, cline=="cal33")), col=2, add=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

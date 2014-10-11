\name{plotExp.cellsurvLQfit}
\alias{plotExp.cellsurvLQfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diagnostic plots of LQ model versus experiments
}
\description{
This function plots the fit of an LQ model versus the fits of each experiment in a series of plots. It is recommended to direct it in a pdf-file.
}
\usage{
\method{plotExp}{cellsurvLQfit}(x, xlim = NULL, ylim = c(0.001, 1.5), xlab = "Dose (Gy)", ylab = "Survival (1 = 100\%)", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{cellsurvLQfit} resulting from function \code{\link{cellsurvLQfit}}.
}
  \item{xlim}{
range of x axis (dose).
}
  \item{ylim}{
range of y-axis.
}
  \item{xlab}{
label of x-axis.
}
  \item{ylab}{
label of y-axis.
}
  \item{\dots}{
further arguments to pass to generic function \code{plotExp}.
}
}
\author{Herbert Braselmann}
\seealso{
\code{cellsurvLQfit}
}
\examples{
datatab<- read.table(system.file("doc", "expl1_cellsurvcurves.txt", package="CFAssay"), header=TRUE, sep="\t")
X<- subset(datatab, cline=="okf6TERT1")
fit<- cellsurvLQfit(X)
print(fit)
#pdf("okf6TERT1_experimental_plots.pdf")
   plotExp(fit)
#dev.off()

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

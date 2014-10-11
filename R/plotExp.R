## S3 Generic function for plotting of data from repeated experiments

plotExp <- function(x, ...) {UseMethod("plotExp")}

plotExp.default <- function(x, ...) {
	if(!class(x) %in% c("cellsurvLQfit","cfa2way"))
	stop("argument not of class 'cellsurvLQfit' or 'cfa2way' resulting from CFAssay functions with the same name.")
	plotExp(x)
}

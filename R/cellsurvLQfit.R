cellsurvLQfit <- function(X, method="ml", PEmethod="fit") {
#Options for method: ml, ls, franken; PEmethod: fit, fixed

	catln <- function(...) {cat(...,"\n")}

# Maximum likelihood, ple fitted
fitLQ.ML <- function(X) {     #LQ model, ple fitted
	X$dose2 <- X$dose^2
	X$lcells <- log(X$ncells)
	dose <- NULL   # to suppress R CMD check notes when looking for visible bindings
	lcells <- NULL   # dito
	fit <- glm(ncolonies ~ factor(Exp) -1 + dose + dose2, offset=lcells,
	  family=quasipoisson(link="log"), data=X)
	#print(summary(fit))
	invisible(fit)
}

# Least squares, ple fitted
fitLQ.LS <- function(X) {
	X$dose2 <- X$dose^2
    X$lnS <- log(X$ncolonies/X$ncells)
	fit <- lm(lnS ~ factor(Exp) -1 + dose + dose2, data=X)
	fit$data <- X
	#print(summary(fit))
	invisible(fit)	
}


# Maximum likelihood, ple fixed
fitLQ1.ML <- function(X) {     #LQ model, ple fixed
	uexp <- unique(X$Exp)  #unique experiment numbers
	X1 <- subset(X, dose>0)
	X1$dose2 <- X1$dose^2
	X1$lcells <- log(X1$ncells)
	dose <- NULL   # to suppress R CMD check notes when looking for visible bindings
	lcells <- NULL   # dito
	logPle <- NULL   # dito
	if (0 %in% X$dose) {
	  X0 <- subset(X, dose==0)
	  ple <- sapply(1:length(uexp), function(i) {ind <- which(X0$Exp==uexp[i]); sum(X0$ncolonies[ind])/(sum(X0$ncells[ind]))})  # pooled 0-dose ratio
	  pleExp <- sapply(1:nrow(X1), function(i) ple[which(uexp==X1$Exp[i])]) #assign to experiment
	  X1$logPle <- log(pleExp)
	}
	if (!(0 %in% X$dose)) X1$logPle <- log(X1$pe)
	fit <- glm(ncolonies ~ -1 + dose + dose2, offset=lcells+logPle, family=quasipoisson(link="log"), data=X1)
	#print(summary(fit))
	invisible(fit)
}

# Least squares, ple fixed
fitLQ1.LS <- function(X) {     #LQ model, ple fixed
	uexp <- unique(X$Exp)  #unique experiment numbers
	X1 <- subset(X, dose>0)
	X1$dose2 <- X1$dose^2
    X1$lnS <- log(X1$ncolonies/X1$ncells)
	dose <- NULL   # to suppress R CMD check notes when looking for visible bindings
	logPle <- NULL   # dito
	if (0 %in% X$dose) {
	  X0 <- subset(X, dose==0)
	  ple <- sapply(1:length(uexp), function(i) {ind <- which(X0$Exp==uexp[i]); sum(X0$ncolonies[ind])/(sum(X0$ncells[ind]))})  # pooled 0-dose ratio
	  pleExp <- sapply(1:nrow(X1), function(i) ple[which(uexp==X1$Exp[i])]) #assign to experiment
	  X1$logPle <- log(pleExp)
	}
	if (!(0 %in% X$dose)) X1$logPle <- log(X1$pe)
	fit <- lm(lnS ~  -1 + dose + dose2, offset=logPle, data=X1)
	#print(summary(fit))
	fit$data <- X1
	invisible(fit)
}

# Least squares, ple fixed, Franken weighted
fitLQ1.LSfr <- function(X) {     #LQ model, ple fixed
	uexp <- unique(X$Exp)  #unique experiment numbers
	X1 <- subset(X, dose>0)
	X1$dose2 <- X1$dose^2
    X1$lnS <- log(X1$ncolonies/X1$ncells)
    X1$wt <- X1$ncolonies*X1$ncells/(X1$ncells - X1$ncolonies)
    #X1$wt <- X1$ncolonies
	dose <- NULL   # to suppress R CMD check notes when looking for visible bindings
	logPle <- NULL   # dito
	wt <- NULL   # dito
	if (0 %in% X$dose) {
	  X0 <- subset(X, dose==0)
	  ple <- sapply(1:length(uexp), function(i) {ind <- which(X0$Exp==uexp[i]); sum(X0$ncolonies[ind])/(sum(X0$ncells[ind]))})
	  pleExp <- sapply(1:nrow(X1), function(i) ple[which(uexp==X1$Exp[i])])
	  X1$logPle <- log(pleExp)
	}
	if (!(0 %in% X$dose)) X1$logPle <- log(X1$pe)
	fit <- lm(lnS ~  -1 + dose + dose2, offset=logPle, weights=wt, data=X1)
	#print(summary(fit))
	fit$data <- X1
	invisible(fit)
}

	if(!method %in% c("ml","ls","franken")) stop("Invalid method, use ml, ls or franken!")
	if(!PEmethod %in% c("fit","fix")) stop("Invalid PE handling string, use fit or fix!")
	if(!(0 %in% X$dose) & !("pe" %in% names(X))) stop("No 0-doses and no pe-column, see help!")
	if(!(0 %in% X$dose)) PEmethod <- "fix"
	if(method=="ml" & PEmethod=="fit") fit <- fitLQ.ML(X)
	if(method=="ml" & PEmethod=="fix") fit <- fitLQ1.ML(X)
	if(method=="ls" & PEmethod=="fit") fit <- fitLQ.LS(X)
	if(method=="ls" & PEmethod=="fix") fit <- fitLQ1.LS(X)
	if(method=="franken") {fit <- fitLQ1.LSfr(X); PEmethod <- "fix"}

	fit$type <- c("cfa", method)
	fit$PEmethod <- PEmethod
	catln("method =", method)
	catln("PEmethod =", PEmethod)
	nr <- length(fit$coef)
	print(fit$coef[(nr-1):nr])
	catln("Use 'print' to see detailed results")
	catln()
	class(fit) <- "cellsurvLQfit"
	invisible(fit)
}

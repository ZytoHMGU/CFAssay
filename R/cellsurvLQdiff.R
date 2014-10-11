## F-test for comparison of LQ-models for two cell survival curves
cellsurvLQdiff <- function(X, curvevar, method="ml", PEmethod="fit") {
	catln <- function(...) {cat(...,"\n")}
	#PEmethod <- "fit"

# Maximum likelihood, ple fitted
fitLQ.MLdiff <- function(X) {     #LQ model, ML method, ple fitted
	X$dose2 <- X$dose^2
	X$lcells <- log(X$ncells)
	lcells <- NULL   # to suppress R CMD check notes when looking for visible bindings
	fit1 <- glm(ncolonies ~ (factor(Exp) -1) %in% curves + dose + dose2, offset=lcells,
	  family=quasipoisson(link="log"), data=X)
	#print(summary(fit1))
	fit2 <- glm(ncolonies ~ (factor(Exp) -1 + dose + dose2) %in% curves, offset=lcells,
	  family=quasipoisson(link="log"), data=X)
	#b <- fit2$coef
	#npar <- length(b)
	#print(b[(npar-3):npar])
	#print(summary(fit2)$coef)
	anv <- anova(fit1,fit2,test="F")
	invisible(list(fit1=fit1, fit2=fit2, anv=anv))
}

# Least squares, ple fitted
fitLQ.LSdiff <- function(X) {     #LQ model, LS method, ple fitted
	X$dose2 <- X$dose^2
    X$lnS <- log(X$ncolonies/X$ncells)
	fit1 <- lm(lnS ~ factor(Exp)*curves + dose + dose2, data=X)
	#print(summary(fit1))
	fit2 <- lm(lnS ~ (factor(Exp) + dose + dose2) %in% curves, data=X)
	anv <- anova(fit1,fit2,test="F")
	invisible(list(fit1=fit1, fit2=fit2, anv=anv))

}

# Maximum likelihood, ple fixed
fitLQ1.MLdiff <- function(X) {     #LQ model, ML method, ple fitted
	X1 <- subset(X, dose>0)
	X1$dose2 <- X1$dose^2
	X1$lcells <- log(X1$ncells)
	dose <- NULL   # to suppress R CMD check notes when looking for visible bindings
	lcells <- NULL   # dito
	if (0 %in% X$dose) {
	  uexp <- unique(X$curveExp)  #unique experiment numbers
	  X0 <- subset(X, dose==0)
ple <- sapply(1:length(uexp), function(i) {ind <- which(X0$curveExp==uexp[i]); sum(X0$ncolonies[ind])/(sum(X0$ncells[ind]))})  # pooled 0-dose ratio
	  pe <- sapply(1:nrow(X1), function(i) ple[which(uexp==X1$curveExp[i])]) #assign to experiment
	  X1$pe <- pe
	}
	#if (!(0 %in% X$dose)) X1$logPle <- log(X1$pe)

	fit1 <- glm(ncolonies ~ -1 + dose + dose2, offset=lcells+log(pe),
	  family=quasipoisson(link="log"), data=X1)
	#print(summary(fit1))
	fit2 <- glm(ncolonies ~ (-1 + dose + dose2) %in% curves, offset=lcells+log(pe),
	  family=quasipoisson(link="log"), data=X1)

	anv <- anova(fit1,fit2,test="F")
	invisible(list(fit1=fit1, fit2=fit2, anv=anv))
}

# Least squares, ple fixed
fitLQ1.LSdiff <- function(X) {     #LQ model, LS method, ple fitted
	X1 <- subset(X, dose>0)
	X1$dose2 <- X1$dose^2
    X1$lnS <- log(X1$ncolonies/X1$ncells)
	dose <- NULL   # to suppress R CMD check notes when looking for visible bindings
	if (0 %in% X$dose) {
	  uexp <- unique(X$curveExp)  #unique experiment numbers
	  X0 <- subset(X, dose==0)
ple <- sapply(1:length(uexp), function(i) {ind <- which(X0$curveExp==uexp[i]); sum(X0$ncolonies[ind])/(sum(X0$ncells[ind]))})  # pooled 0-dose ratio
	  pe <- sapply(1:nrow(X1), function(i) ple[which(uexp==X1$curveExp[i])]) #assign to experiment
	  X1$pe <- pe
	}
	#if (!(0 %in% X$dose)) X1$logPle <- log(X1$pe)

	fit1 <- lm(lnS ~ -1 + dose + dose2, offset=log(pe), data=X1)
	#print(summary(fit1))
	fit2 <- lm(lnS ~ (-1 + dose + dose2) %in% curves, offset=log(pe), data=X1)
	anv <- anova(fit1,fit2,test="F")
	invisible(list(fit1=fit1, fit2=fit2, anv=anv))
}

# Least squares, ple fixed, Franken weighted
fitLQ1.LSdiffFr <- function(X) {     #LQ model, LS method, ple fitted
	X1 <- subset(X, dose>0)
	X1$dose2 <- X1$dose^2
    X1$lnS <- log(X1$ncolonies/X1$ncells)
    X1$wt <- X1$ncolonies*X1$ncells/(X1$ncells - X1$ncolonies)
	dose <- NULL   # to suppress R CMD check notes when looking for visible bindings
	wt <- NULL   # dito
	if (0 %in% X$dose) {
	  uexp <- unique(X$curveExp)  #unique experiment numbers
	  X0 <- subset(X, dose==0)
ple <- sapply(1:length(uexp), function(i) {ind <- which(X0$curveExp==uexp[i]); sum(X0$ncolonies[ind])/(sum(X0$ncells[ind]))})  # pooled 0-dose ratio
	  pe <- sapply(1:nrow(X1), function(i) ple[which(uexp==X1$curveExp[i])]) #assign to experiment
	  X1$pe <- pe
	}
	#if (!(0 %in% X$dose)) X1$logPle <- log(X1$pe)

	fit1 <- lm(lnS ~ -1 + dose + dose2, offset=log(pe), weights=wt, data=X1)
	#print(summary(fit1))
	fit2 <- lm(lnS ~ (-1 + dose + dose2) %in% curves, offset=log(pe), weights=wt, data=X1)
	anv <- anova(fit1,fit2,test="F")
	invisible(list(fit1=fit1, fit2=fit2, anv=anv))
}


	if(!method %in% c("ml","ls","franken")) stop("Invalid method, use ml, ls or franken!")
	if(!PEmethod %in% c("fit","fix")) stop("Invalid PE handling string, use fit or fix!")
	if(!(0 %in% X$dose) & !("pe" %in% names(X))) stop("No 0-doses and no pe-column, see help!")
	
	if(!curvevar %in% names(X)) stop("curvevar not found in data frame!")
	X$curves <- X[,curvevar]
	X$curveExp <- paste(X$curves, X$Exp, sep="_")

	if(!(0 %in% X$dose)) PEmethod <- "fix"
	if(method=="ml" & PEmethod=="fit") fitcomp <- fitLQ.MLdiff(X)
	if(method=="ml" & PEmethod=="fix") fitcomp <- fitLQ1.MLdiff(X)
	if(method=="ls" & PEmethod=="fit") fitcomp <- fitLQ.LSdiff(X)
	if(method=="ls" & PEmethod=="fix") fitcomp <- fitLQ1.LSdiff(X)
	if(method=="franken") {fitcomp <- fitLQ1.LSdiffFr(X); PEmethod <- "fix"}

	catln("*** Overall comparison for two linear-quadratic cell survival curves ***")
	catln("Compared curves:", levels(X$curves))
	catln("method:", method)
	catln("PEmethod:", PEmethod)
	catln()
	catln("Test used: F-test")
	anvtab <- fitcomp$anv
	  rownames(anvtab) <- c("","values")  
	  print(anvtab[2,5:6])
	catln("Use 'print' to see detailed results")
	catln()
	#catln("Coefficients for Model 2: logarithmic PEs, dose and dose-squared coefficients")
	#print(summary(fitcomp$fit2))
	#print(fitcomp$anv)

	fitcomp$type <- c("cfa", method)
	fitcomp$PEmethod <- PEmethod
	class(fitcomp) <- "cellsurvLQdiff"
	invisible(fitcomp)
}

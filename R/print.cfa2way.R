print.cfa2way <- function(x, labels=c(A="A",B="B"), ...) {
	catln <- function(...) {cat(...,"\n")}

	print.fitval <- function(fit, type) {
	if (type == "ml") {
    ssqwresTot <- sum(fit$residuals^2*fit$weight)  #total sum of squared weighted residuals  
	catln("Goodness-of-fit values")
	catln("	Residual Deviance:", fit$deviance)
	catln("	Total sum of squared weighted residuals ssqwresTot:", ssqwresTot)
	catln("	Residual Degrees of Freedom:", fit$df.residual)
	catln("	Dispersion parameter:", summary(fit)$dispersion)
	}
	
	if (type == "ls") {
	#if (nrow(fit$data)!=nrow(fit$model)) stop("Error: more data rows than in the model!")
    ssqresTot <- sum(fit$residuals^2)  #total sum of squared residuals  
	catln()
	catln("Goodness-of-fit values")
	catln("	Total sum of squared residuals ssqresTot:", ssqresTot)
	catln("	Residual Degrees of Freedom:", fit$df.residual)
	catln("	Multiple R-squared:", summary(fit)$r.squared)
	}

	}

	fitcomp <- x
	if(!class(fitcomp)=="cfa2way") stop("Fit object not of class 'cfa2way'!")
	#if (!("cfa2way" %in% fitcomp$type)) stop("Error: fit object not of type cfa2way!")

	fit <- fitcomp$fit2
	sfit <- summary(fit)
	coef <- sfit$coef
	  coef <- cbind(coef, surv_percent=round(exp(coef[,1])*100,1))
	coef1 <- summary(fitcomp$fit1)$coef
	  coef1 <- cbind(coef1, surv_percent=round(exp(coef1[,1])*100,1))
	nr <- nrow(coef)
	nr1 <- nrow(coef1)
	#nam <- sub("factor\\(Exp\\)","PE",rownames(coef))
	#nam[(nr-1):nr] <- c("alpha","beta")
	#rownames(coef) <- nam

	catln("*** Logarithmic linear two-way ANOVA for factors A and B with interaction ***")
	catln("=============================================================================")
	catln("A=", labels["A"], ", B=", labels["B"])
	catln("Postscript digits for A or B: 0 inactive, 1 active")
	catln("surv_percent = exp(Estimate)*100")
	catln()

	catln("Null hypothesis (Model 1): no interaction")
	catln("-----------------------------------------")
	idx <- grep("Exp",rownames(coef1))
	print(coef1[-idx,])
	catln()
	print.fitval(fitcomp$fit1, fitcomp$type[3])
	catln()

	catln("Alternative hypothesis (Model 2): interaction")
	catln("---------------------------------------------")
	catln("parametrization:", fitcomp$type[2])
	idx <- grep("Exp",rownames(coef))
	print(coef[-idx,])
	catln()
	print.fitval(fit, fitcomp$type[3])
	catln()

	anv <- fitcomp$anv
	attr(anv, "heading") <- c("Analysis of Variance Table and F-test", "Model 2 versus Model 1")
	print(anv)

}

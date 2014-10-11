print.cellsurvLQdiff <- function(x, ...) {   #fitcomp contains objects fit1, fit2, anv
	
##### Function definition part #####
	catln <- function(...) {cat(...,"\n")}

	print.fitval <- function(fit, type) {
	if (type == "ml") {
    rsswTot <- sum(fit$residuals^2*fit$weight)  #total sum of squared weighted residuals  
	catln("Goodness-of-fit values")
	catln("	Residual Deviance:", fit$deviance)
	catln("	Total sum of squared weighted residuals rsswTot:", rsswTot)
	catln("	Residual Degrees of Freedom:", fit$df.residual)
	catln("	Dispersion parameter:", summary(fit)$dispersion)
	}
	
	if (type == "ls") {
	#if (nrow(fit$data)!=nrow(fit$model)) stop("Error: more data rows than in the model!")
    rssTot <- sum(fit$residuals^2)  #total sum of squared residuals  
	catln()
	catln("Goodness-of-fit values")
	catln("	Total sum of squared residuals rssTot:", rssTot)
	catln("	Residual Degrees of Freedom:", fit$df.residual)
	catln("	Multiple R-squared:", summary(fit)$r.squared)
	}

	if (type == "franken") {
    rsswTot <- sum(fit$residuals^2*fit$weight)  #total sum of squared weighted residuals  
	catln("Goodness-of-fit values")
	catln("	Total sum of squared weighted residuals rsswTot:", rsswTot)
	catln("	Residual Degrees of Freedom:", fit$df.residual)
	catln("	Multiple R-squared:", summary(fit)$r.squared)
	}
	}

##### Execution part #####
	fitcomp <- x
	if(!class(fitcomp)=="cellsurvLQdiff") stop("Fit object not of class 'cellsurvLQdiff'!")
	#if (!("cfa" %in% fitcomp$type)) stop("fit object not of type 'cfa'!")
	#if (fitcomp$PEmethod != "fit") stop("PEmethod not 'fit'. Use standard R, print(fit) or print(summary(fit))!")
	
	fit <- fitcomp$fit2
	sfit <- summary(fit)
	coef <- sfit$coef
	coef1 <- summary(fitcomp$fit1)$coef
	nr <- nrow(coef)
	nr1 <- nrow(coef1)
	nam <- sub("factor\\(Exp\\)","PE",rownames(coef))
	  nam <- sub("dose2","beta",nam)
	  nam <- sub("dose","alpha",nam)
	rownames(coef) <- nam
	idx <- grep("dose",rownames(coef1))
	nam1 <- sub("dose2","beta",rownames(coef1))
	  nam1 <- sub("dose","alpha",nam1)
	rownames(coef1) <- nam1
	catln("Overall comparison test for coefficients alpha and beta of LQ-models")
	catln("====================================================================")
	catln("method =", fitcomp$type[2])
	catln("PEmethod =", fitcomp$PEmethod)
	catln()

	if (fitcomp$PEmethod == "fit") catln(nr-4, "PEs fitted as intercepts. To look at, use simple R print function.")
	catln("Null hypothesis (Model 1): one set of shape parameters alpha and beta for all data")
	catln("----------------------------------------------------------------------------------")
	print(coef1[idx,])
	catln()
	print.fitval(fitcomp$fit1, fitcomp$type[2])
	catln()

	catln("Alternative hypothesis (Model 2): two sets of shape parameters alpha and beta")
	catln("-----------------------------------------------------------------------------")
	print(coef[(nr-3):nr,])
	catln()
	print.fitval(fit, fitcomp$type[2]) #fit is fitcomp$fit2, see lines before
	
	catln()
	anv <- fitcomp$anv
	attr(anv, "heading") <- c("Analysis of Variance Table and F-test", "Model 2 versus Model 1")
	print(anv)

	
}

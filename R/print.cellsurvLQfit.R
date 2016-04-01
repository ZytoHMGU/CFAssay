print.cellsurvLQfit <- function(x, ...) {
	catln <- function(...) {cat(...,"\n")}
	as.lm <- function(x) {class(x) <- "lm"; x}
	as.glm <- function(x) {class(x) <- "glm"; x}

	dose <- NULL   # to suppress R CMD check notes when looking for visible bindings

	fit <- x
	if(!class(fit)=="cellsurvLQfit") stop("Fit object not of class 'cellsurvLQfit'!")
	#if (!("cfa" %in% fit$type)) stop("fit object not of type 'cfa'!")
	#if (fit$PEmethod != "fit") stop("PEmethod not 'fit'. Use standard R, print(fit) or print(summary(fit))!")
	
  data0 <- subset(fit$data, dose==0)
	if("ml" %in% fit$type) sfit <- summary(as.glm(fit))
	if("ls" %in% fit$type | "franken" %in% fit$type) sfit <- summary(as.lm(fit))
	coef <- sfit$coef
	nr <- nrow(coef)
	uexp <- unique(fit$data$Exp)  #unique experiment numbers
  
## Coeffients
	catln("*** Coefficients of LQ-model for cell survival ***")
	catln("method =", fit$type[2])
	catln("PEmethod =", fit$PEmethod)
	catln()

	if(fit$PEmethod == "fit") {
	uexpfit<- sapply(1:(nr-2), function(i) strsplit(rownames(coef)[i], ")")[[1]][2])
	nam <- sub("factor\\(Exp\\)","PE",rownames(coef))
	nam[(nr-1):nr] <- c("alpha","beta")
	rownames(coef) <- nam
	#catln("Coefficients alpha and beta of LQ-model and fitted logarithmic plating efficiencies PE")
	catln("Logarithmic plating efficiencies PE fitted as intercepts")
	catln("see remark in the manual, 1.2")
	print(coef[1:(nr-2),])
	catln()
	catln("Shape parameters alpha and beta")
	print(coef[(nr-1):nr,])
	catln()
	catln("Observed and fitted plating efficiencies (%):")
	plefit <- exp(coef[1:(nr-2),1])
	ple <- sapply(1:length(uexpfit), function(i) {ind <- which(data0$Exp==uexpfit[i]); sum(data0$ncolonies[ind])/(sum(data0$ncells[ind]))})
	if(!all(sort(uexp) == sort(uexpfit))) stop("Mismatch of experiment names and coefficient names!")
	print(data.frame(Experiment=uexpfit, PE=round(ple*100,1), PEfitted=round(plefit*100,1)))
	}

	if(fit$PEmethod == "fix") {
	  rownames(coef) <- c("alpha", "beta")
	  catln("Shape parameters alpha and beta")
	  print(coef)
	}

## Quality parameters ml or franken
	if ("ml" %in% fit$type | "franken" %in% fit$type) {
    rsswTot <- sum(fit$residuals^2*fit$weight)  #total sum of squared weighted residuals  
	catln()
	if("ml" %in% fit$type) catln("Residual Deviance:", fit$deviance)
	catln("Total residual sum of weighted squares rsswTot:", rsswTot)
	catln("Residual Degrees of Freedom:", fit$df.residual)
	if("ml" %in% fit$type) catln("Dispersion parameter:", sfit$dispersion)
	if("franken" %in% fit$type) catln("Multiple R-squared:", sfit$r.squared)
	rssw <- sapply(1:length(uexp), function(i) {ind <- which(fit$data$Exp==uexp[i]); sum(fit$residuals[ind]^2*fit$weight[ind])})
	catln()
	catln("Fraction rssw of rsswTot per Experiment")
	print(data.frame(Experiment=uexp, rssw=round(rssw,2), perCent=round(rssw/rsswTot*100,1)))
	}

## Quality parameters ls	
	if ("ls" %in% fit$type) {
	if (nrow(fit$data)!=nrow(fit$model)) stop("More data rows than in the model!")
    rssTot <- sum(fit$residuals^2)  #total sum of squared residuals  
	catln()
	catln("Total residual sum of squares rssTot:", rssTot)
	catln("Residual Degrees of Freedom:", fit$df.residual)
	catln("Multiple R-squared:", sfit$r.squared)
	rss <- sapply(1:length(uexp), function(i) {ind <- which(fit$data$Exp==uexp[i]); sum(fit$residuals[ind]^2)})
	catln()
	catln("Fraction rss of rssTot per Experiment")
	print(data.frame(Experiment=uexp, rss=round(rss,2), perCent=round(rss/rssTot*100,1)))
	}
	
}

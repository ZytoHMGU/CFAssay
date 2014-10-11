plotExp.cellsurvLQfit <- function(x, xlim=NULL, ylim=c(0.001, 1.5), xlab="Dose (Gy)", ylab="Survival (1 = 100%)", ...) {
	fit <- x
	if (!"cfa" %in% fit$type) stop("Fit object not of type cfa!")
	#if (fit$PEmethod != "fit") stop("PEmethod not 'fit'. Experimental not yet implemented.")

	X <- fit$data
	doses <- unique(X$dose)
	maxd <-max(doses)
	b <- fit$coef
	npar <- length(b)
	uexp <- unique(X$Exp)  #unique experiment numbers
	if(fit$PEmethod=="fit") names(b)[1:length(uexp)] <- uexp

	par(lwd=2)

	lcells <- NULL   # to suppress R CMD check notes when looking for visible bindings
	Exp <- NULL   # dito
	logPle <- NULL   # dito
	wt <- NULL   # dito


##### ML method #####
	if("ml" %in% fit$type) {
    rsswTot <- sum(fit$residuals^2*fit$weight)  #total sum of squared weighted residuals
    #main loop
	for (i in 1:length(uexp)) {
	  ind <- which(X$Exp==uexp[i])            #indices of data rows in X for ith experiment
	  if(fit$PEmethod=="fit") {
	    fiti <- glm(ncolonies ~ dose + dose2, offset=lcells, family=quasipoisson(link="log"), data=subset(X, Exp==uexp[i]))
	    pe <- exp(b[i])     #fitted plating efficiency of ith experiment
	  }
	  if(fit$PEmethod=="fix") {
	    fiti <- glm(ncolonies ~ -1 + dose + dose2, offset=lcells+logPle, family=quasipoisson(link="log"), data=subset(X, Exp==uexp[i]))
	    pe <- exp(X$logPle[X$Exp==uexp[i]][1])
	  }
	  
	  bi <- fiti$coef
	  dp <- summary(fiti)$dispersion

	  sf <- X$ncolonies[ind]/X$ncells[ind]/pe #survival fractions of ith experiment
	  sde.poisson <- sqrt(fit$fitted.values[ind])/X$ncells[ind]/pe  #fitted Poisson s.d.

	  rssw <- sum(fit$residuals[ind]^2*fit$weight[ind])
	  do <- X$dose[ind]
	  sf01 <- sf
	    if(fit$PEmethod=="fix") {do <- c(0, do); sf01 <- c(1, sf)} #formal, otherwise error (!?)

	  plot(do, sf01, log="y", ylim=ylim, col=i, xlab=xlab, ylab=ylab)  #points
	  segments(X$dose[ind], sf+sde.poisson, X$dose[ind], sf-sde.poisson, col=i) #error bars
	  
	  #curve i
	  if(fit$PEmethod=="fit")
	    curve(exp(bi[1] + bi[2]*x + bi[3]*x^2)/pe, from=0, to=max(X$dose), col=i, add=TRUE)
	  if(fit$PEmethod=="fix")
	    curve(exp(bi[1]*x + bi[2]*x^2), from=0, to=max(X$dose), col=i, add=TRUE)
	  #mean curve  
	  curve(exp(b[npar-1]*x + b[npar]*x^2), from=0, to=max(X$dose), add=TRUE)

	  legend(0.62*maxd, ylim[2], c(paste("Replicate curve Exp",uexp[i]), "Mean curve"), text.col=c(i,1))
	  text(0, ylim[1], pos=4, paste("rssw", round(rssw,1), "of", round(rsswTot,1), " ", round(rssw/rsswTot*100,1), "%\n Dispersion parameter of experiment:", round(dp,2)))
	  title(main=paste("ML-fit, mean curve versus Exp", uexp[i]))
	  #cat("plot experiment",uexp[i], ", rssw", rssw, "\n")
	} #end of for-loop
	} #end of ml-method


##### LS method or Franken et al. (Nature Prot. 2006) #####
	if("ls" %in% fit$type | "franken" %in% fit$type) {

    if("ls" %in% fit$type)
      rssTot <- sum(fit$residuals^2)  #total sum of squared residuals  
    if("franken" %in% fit$type)
      rsswTot <- sum(fit$residuals^2*fit$weight)  #total sum of squared weighted residuals  

	for (i in 1:length(uexp)) {
	  ind <- which(X$Exp==uexp[i])      #indices of data rows in X for ith experiment

	  if(fit$PEmethod=="fit") {
	    fiti <- lm(lnS ~ dose + dose2, data=subset(X,Exp==uexp[i]))
	    pe <- exp(b[i])     #fitted plating efficiency of ith experiment
		rss <- sum(fit$residuals[ind]^2)
	  }

	  if(fit$PEmethod=="fix") {
	  	if("ls" %in% fit$type) {
	      fiti <- lm(lnS ~ -1 + dose + dose2, offset=logPle, data=subset(X,Exp==uexp[i]))
		  rss <- sum(fit$residuals[ind]^2)
	  	}

	    if("franken" %in% fit$type) {
	      fiti <- lm(lnS ~ -1 + dose + dose2, offset=logPle, weights=wt, data=subset(X,Exp==uexp[i]))
		  rssw <- sum(fit$residuals[ind]^2*fit$weight[ind])
	    }  #wt in fit$data

	    pe <- exp(X$logPle[X$Exp==uexp[i]][1])
	  }

	  bi <- fiti$coef
	  sf <- X$ncolonies[ind]/X$ncells[ind]/pe #survival fractions of ith experiment
	  do <- X$dose[ind]
	  sf01 <- sf
	  if(fit$PEmethod=="fix") {do <- c(0, do); sf01 <- c(1, sf)}
	  plot(do, sf01, log="y", ylim=ylim, col=i, xlab=xlab, ylab=ylab)  #points
	  
	  #curve i
	  if(fit$PEmethod=="fit")
	    curve(exp(bi[1] + bi[2]*x + bi[3]*x^2)/pe, from=0, to=max(X$dose), col=i, add=TRUE)
	  if(fit$PEmethod=="fix")
	    curve(exp(bi[1]*x + bi[2]*x^2), from=0, to=max(X$dose), col=i, add=TRUE) #curve i
	  #mean curve
	  curve(exp(b[npar-1]*x + b[npar]*x^2), from=0, to=max(X$dose), add=TRUE)

	  legend(0.62*maxd, ylim[2], c(paste("Replicate curve Exp",uexp[i]), "Mean curve"), text.col=c(i,1))
	  if("ls" %in% fit$type) {ltext <- "rss"; rs <- rss; rsTot <- rssTot}
	  if("franken" %in% fit$type) {ltext <- "rssw"; rs <- rssw; rsTot <- rsswTot}
	  text(0, ylim[1], pos=4, paste(ltext, round(rs,2), "of", round(rsTot,2), " ", round(rs/rsTot*100,3),"%"))
	  if("ls" %in% fit$type) title(main=paste("LS-fit, mean curve versus Exp", uexp[i]))
	  if("franken" %in% fit$type) title(main=paste("WLS-fit, mean curve versus Exp", uexp[i]))
	  #cat("plot experiment",uexp[i], ", rss", rss, "\n")
	} #end of for-loop
	} #end of ls-method
	
}

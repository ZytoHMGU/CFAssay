plotExp.cfa2way <- function(x, labels=c(A="A",B="B"), ...) {
	fitcomp <- x
	if(!class(fitcomp)=="cfa2way") stop("Fit object not of class 'cfa2way'!")

	Exp <- NULL   # to suppress R CMD check notes when looking for visible bindings

	fit <- fitcomp$fit2
	data <-cbind(fit$data, fitted=fit$fitted)
    ssqwresTot <- sum(fit$residuals^2*fit$weight)  #total sum of squared weighted residuals  
	if(fitcomp$type[3]=="ls")
      ssqresTot <- sum(fit$residuals^2)  #total sum of squared residuals  
	Exps <- as.vector(unique(data$Exp))
	par(lwd=2)
	for (i in Exps) {
	  ind <- which(data$Exp == i)
	  dsub <- subset(data, Exp == i)
	  # sum because potentially more than one observation per "dose"
	  y <- as.vector(by(dsub, list(dsub$B,dsub$A), function(x) sum(x$ncolonies)/sum(x$ncells)))
	  y1 <- as.vector(by(dsub, list(dsub$B,dsub$A), function(x) sum(x$fitted)/sum(x$ncells)))
	  if(fitcomp$type[3]=="ls")
	    y1 <- as.vector(by(dsub, list(dsub$B,dsub$A), function(x) exp(mean(x$fitted))))
	  if(fitcomp$type[3]=="ml") 	  	
	    sde.poisson <- as.vector(by(dsub, list(dsub$B,dsub$A), function(x) sqrt(sum(x$fitted))/sum(x$ncells)))
	  

	  ssqwres <- sum(fit$residuals[ind]^2*fit$weight[ind])
	  if(fitcomp$type[3]=="ls")
	    ssqres <- sum(fit$residuals[ind]^2)

	  barplot(y, ylim=c(0,1), width=2/3, space=0.5, names.arg=c("Control",labels["B"],labels["A"],paste(labels["A"],"+",labels["B"])), ylab="number of colonies per cells")
	  points((0:3)+2/3, y1)
	  if(fitcomp$type[3]=="ml")
	    segments((0:3)+2/3, y1+sde.poisson, (0:3)+2/3, y1-sde.poisson, col=1) #error bars
	  text(2.7,0.92, paste("bars: observed","points: modelled means",sep="\n"), pos=4)
	  title(main=paste("Observed vs modelled means, Experiment", i))
	  if(fitcomp$type[3]=="ml")
	    text(4, 0.85, pos=2, paste("ssqwres", round(ssqwres,1), "of", round(ssqwresTot,1), " ", round(ssqwres/ssqwresTot*100,1),"%"))
	  if(fitcomp$type[3]=="ls")
	    text(4, 0.85, pos=2, paste("ssqres", round(ssqres,3), "of", round(ssqresTot,3), " ", round(ssqres/ssqresTot*100,1),"%"))
	    

	  box()
	}

}

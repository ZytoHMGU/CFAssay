plot.cellsurvLQfit <- function(x, xlim=NULL, ylim=c(0.008, 1.0), xlab="Dose (Gy)", ylab="Survival (1 = 100%)", col=1, pch=1, add=FALSE, ...) {

	fit <- x
	if(!class(fit)=="cellsurvLQfit") stop("Fit object not of class 'cellsurvLQfit'!")
	#if (!"cfa" %in% fit$type) stop("Error: fit object not of type cfa!")
	data <- fit$data
	data$dose2 <- data$dose^2
	data$lcells <- log(data$ncells)  #falls LS
	uexp <- unique(data$Exp)
	doses <- unique(data$dose)
	maxd <-max(doses)
	b <- fit$coef[c("dose","dose2")]
	if(is.null(xlim)) xlim <- c(0,maxd)
	par(lwd=2)
	if(!add)
	  curve(exp(b[1]*x + b[2]*x^2), from=0, to=maxd, log="y", col=col, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab)
	if(add)
	  curve(exp(b[1]*x + b[2]*x^2), from=0, to=maxd, col=col, add=TRUE)
	if(0 %in% doses) {
    #b <- fit$coef
	#S0 <- exp(b[1:length(unique(fit$data$Exp))]); S0
    #nc <- length(unique(data$Exp))
    # Plating efficiencies from seperate experiment fits (ML, possibly without 0-dose data):
	S0 <- pes(data)$S0  #CFAssay function pes
  names(S0) <- rownames(pes(data))
	meanSF <- sfpmean(data, S0)  #CFAssay function sfpmean
	}
	if(!(0 %in% doses)) {
	  data$pe <- exp(data$logPle)
	  meanSF <- sfpmean(data)
	}
	pts <- meanSF[1,]
	sems <- meanSF[2,]
	points(doses, pts, col=col, pch=pch)
	segments(x0=doses, y0=pts-sems, x1=doses, y1=pts+sems, col=col)
}

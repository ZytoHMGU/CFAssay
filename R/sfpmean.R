
sfpmean <- function(X, S0=NULL)
## Pointwise mean survival fractions for curves with several experimental replicates
{

sf.mean <- function(X)
## Point mean for replicates with different plating efficiencies
# X must have a column "S0" !
{
	v <- X$ncells*X$S0
	fit <- glm(X$ncolonies ~  v -1, family=quasipoisson(link="identity"),data=X)
	if (summary(fit)$dispersion<1)
	  fit <- glm(X$ncolonies ~  v -1, family=poisson(link="identity"),data=X)
	summary(fit)$coef
}

	dose <- NULL   # to suppress R CMD check notes when looking for visible bindings

	if(is.null(S0) & !("pe" %in% names(X)))
	  stop("In function sfpmean: 'S0' has to be given in argument list or 'pe' has to be a column in the data frame! See help document.")
	  
	doses <- unique(X$dose)
	
	if(is.null(S0)) {
	  #S0 <- X$pe
	#pleExp <- sapply(1:nrow(X1), function(i) ple[which(uexp==X1$Exp[i])]) #assign to experiment
	  pmean <- sapply(1:length(doses), function(i) {X1 <- subset(data.frame(X, S0=X$pe), dose==doses[i]); sf.mean(X1)})
	}
	if(!is.null(S0)) { #attention for change of S0!
    if(is.null(names(S0))) stop("S0 is not a named vector!")
    if(length(grep("(Exp)", names(S0))) == length(unique(X$Exp)))
      {ExpNames <- sapply(1:length(S0), function(i) strsplit(names(S0)[i], ")")[[1]][2])} else
      {ExpNames <- names(S0)}
  if(!all(sort(unique(X$Exp)) == sort(ExpNames))) stop("Mismatch of experiment names in data frame and names of S0!")
  S01 <- sapply(1:nrow(X), function(i) S0[X$Exp[i]])
  X$S0 <- S01  # attention: X$S0 not S0
  pmean <- sapply(1:length(doses), function(i) {X1 <- subset(X, dose==doses[i]); sf.mean(X1)})
	}
	colnames(pmean) <- paste("dose_",doses,sep="")
	rownames(pmean) <- 1:nrow(pmean)
	rownames(pmean)[1:2] <- c("SF", "stdev")
	pmean[1:2,]
}

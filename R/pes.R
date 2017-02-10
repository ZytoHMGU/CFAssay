
pes <- function(X)
## Calculation of plating efficiencies for curve data with several experimental replicates
{
	uexp <- unique(X$Exp)
	nexp<- length(uexp)
	X$dose2 <- X$dose^2
	X$lcells <- log(X$ncells)
	lcells <- NULL   # to suppress R CMD check notes when looking for visible bindings
	Exp <- NULL   # dito
	#repeated measurements for dose==0 are pooled (sum)
	if (0 %in% X$dose) {
	  pe <- sapply(1:nexp, function(i) {ind <- which(X$dose==0 & X$Exp== uexp[i]); sum(X$ncolonies[ind])/(sum(X$ncells[ind]))})
	#plating efficiencies S0 fitted by curves of repeated experiments
	  S0 <- sapply(1: nexp, function(i) glm(ncolonies ~ dose + dose2, offset=lcells, family=quasipoisson(link="log"), data=subset(X, Exp==uexp[i]))$coef[1])
	  S0 <- exp(S0)
	}
	if (!(0 %in% X$dose) & "pe" %in% names(X)) {
	  pe <- X$pe
	  X$logPle <- log(pe)
	  X$lnS <- log(X$ncolonies/X$ncells)
	  #plating efficiencies S0 fitted by curves of repeated experiments
	  S0 <- sapply(1: nexp, function(i) lm(lnS ~ -1 + dose + dose2, offset = logPle, family=quasipoisson(link="log"), data=subset(X, Exp==uexp[i]))$coef[1])
	  S0 <- exp(S0)
	}
	indu <- which(!duplicated(X$Exp))   #included 2017/02/09
	sf0 <- data.frame(Exp=uexp, pe=pe[indu], S0=S0)
	rownames(sf0) <- uexp
	sf0
}

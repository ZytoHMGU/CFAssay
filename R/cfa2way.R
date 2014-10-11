cfa2way <- function(X, A, B, param="A/B", method="ml") {
	catln <- function(...) {cat(...,"\n")}
	if(!param %in% c("A/B","B/A","A*B")) stop("param should be one of A/B, B/A or A*B")
	if(!method %in% c("ml","ls")) stop("Invalid method, use ml or ls!")

	nameA <- A
	nameB <- B
	X$A <- factor(X[,nameA])
	X$B <- factor(X[,nameB])
	X$Exp <- factor(X$Exp)  #to ensure that Exp is factor
	
	if(method=="ml") {
	  X$lcells <- log(X$ncells)
	  fo <- formula(paste("ncolonies ~ Exp -1 +", param))
	  lcells <- NULL   # to suppress R CMD check notes when looking for visible bindings
	  fit1 <- glm(ncolonies ~ Exp - 1 + A + B, family=quasipoisson(link="log"), offset=lcells, data=X)
	  fit2 <- glm(fo, family=quasipoisson(link="log"), offset=lcells, data=X)
	  anv <- anova(fit1,fit2,test="F")
	}

	if(method=="ls") {
      X$lnS <- log(X$ncolonies/X$ncells)
	  fo <- formula(paste("lnS ~ Exp -1 +", param))
	  fit1 <- lm(lnS ~ Exp - 1 + A + B, data=X)
	  fit2 <- lm(fo, data=X)
	  anv <- anova(fit1,fit2,test="F")
	  fit2$data <- X
	}

	catln("*** Two-way ANOVA for factors A and B with interaction ***")
	catln("A=", nameA, ", B=", nameB)
	catln("Test for interaction: F-test")
	anvtab <- anv
	  rownames(anvtab) <- c("","values")  
	  print(anvtab[2,5:6])
	catln("Use 'print' to see detailed results")
	catln()

	#catln("A=", nameA, ", B=", nameB)
	#coef <- summary(fit2)$coef
	#  idx <- grep("Exp",rownames(coef))
	#print(coef[idx,])
	#catln()
	#print(coef[-idx,])
	#catln("ANOVA F-Test\n0-hypothesis (Model 1) no interaction\nalternative hypothesis (Model 2) interaction")
	#catln()
	#attr(anv, "heading") <- c("Analysis of Variance Table", "Model 2 versus Model 1")
	#print(anv)
	fitcomp <- list(fit1=fit1, fit2=fit2, anv=anv)
	
	fitcomp$type <- c("cfa2way", param, method)
	class(fitcomp) <- "cfa2way"
    invisible(fitcomp)
}

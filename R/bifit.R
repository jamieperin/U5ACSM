
bifit <- function(x, y, c=6, tol=1e-6, intercept = TRUE, tolerance = 1e-07, yname = NULL) {
	# Uses lsfit() to compute WLS using Tukey's biweight function
	# c is the tuning constant, typically around 6 to 9
	# tol is for convergence of the biweight fit; tolerance is for the matrix decomposition in lsfit
	# Like lsfit, input x matrix should not include a column of ones
	# Like lsfit, input y can be matrix if there are multiple left-hand sides

	if (is.vector(y)) {
		coef.old <- iter <- 0; coef.new <- 1; wt <- NULL
		while (max(abs(coef.new-coef.old)) > tol) {
			iter <- iter+1; #print(iter)
			z <- lsfit(x, y, wt, intercept=intercept, tolerance=tolerance, yname=yname)
			S <- max(median(abs(z$resid)), .1)
               			u <- z$resid / (30*S)
			coef.old <- coef.new
			coef.new <- z$coef }

		resid <- z$resid
		names(resid) <-  names(u) <- names(y)
		names(coef.new) <- paste("b", 0:(ncol(cbind(1,x))-1), sep="")
		return(list(coef=coef.new, residuals=z$resid, intercept=intercept, wt=wt, u=u)) }

	if (is.matrix(y)) {
		resid <- wt <- u <- coef <- NULL
		for (j in 1:ncol(y)) {
			print(paste("Y",j,":",sep=""))
			z <- bifit(x, y[,j], c=c, tol=tol, intercept=intercept, tolerance=tolerance, yname=yname)

                        #################################################################
			resid <- cbind(resid, z$resid)
			wt <- cbind(wt, z$wt)
			u <- cbind(u, z$u)
			coef <- cbind(coef, z$coef) }
		dimnames(resid) <- dimnames(u) <- dimnames(y) # <- dimnames(wt)
		dimnames(coef) <- list(paste("b",0:(ncol(cbind(1,x))-1),sep=""),dimnames(y)[[2]])
		return(list(coef=coef, residuals=resid, intercept=intercept, wt=wt, u=u)) } }

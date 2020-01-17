mxhat.logquad <- function(coefs, sex, Q5, k=rep(0,length(Q5))) {

	# Inputs:   Q5, sex, k (optional; default is k=0), and model coefficients
	#           Q5 and k can be scalar or vector, but must have same length
	#           coefs contains coefficients ax, bx, cx, and vx
	# Outputs:  Vector (or matrix) of predicted sex-specific mx values
	#           with age groups 0, 1-4, 5-9, ..., 110+
	#ages <- c("0", "1-4", paste(seq(5,105,5), seq(9,109,5), sep="-"), "110+")
#  cat("LOGQUAD\n")
  ages<- c("0-6day","7-27day","1-5m","6-11m","12-23m","24-59m")

	if (length(Q5)!=length(k)) { print("error: Q5 and k input vectors must have same length"); break }
	#if (!is.array(coefs)) { print("Error: table of coefficients must be an array"); break }
	ax <- coefs[, "ax"]; bx <- coefs[, "bx"]; cx <- coefs[, "cx"]; vx <- coefs[, "vx"]
	h <- log(Q5)
	if (length(Q5)==1) {
		# Compute age-specific mx from h, h^2, and coefficients
		mx <- exp(ax + bx*h + cx*h^2 + vx*k)
#		cat("\nmx\n")
#		print(mx)
		# Force 4q1 (and thus 4m1) to be consistent with 1q0 and 5q0
#		a1 <- coale.demeny.a0 (mx[1], sex)
#		a4 <- coale.demeny.4a1(mx[1], sex)
#		Q1 <- mx[1] / ( 1 + (1-a4)*mx[1] )
#		Q4 <- 1 - (1-Q5)/(1-Q1)
		#mx[2] <- Q4 / ( 4 - (4-a4)*Q4 )

            #Force cumulative qx to be consistent
            for(ii in length(mx):2){
                if(mx[ii] < mx[ii - 1]) mx[ii-1] <- mx[ii]
                }

        }
	else {
		mx <- matrix(NA, length(ages), length(Q5))
		dimnames(mx) <- list(ages, names(Q5))
		for (j in 1:length(Q5)) mx[,j] <- mxhat.logquad(coefs, sex, Q5[j], k[j]) }
	return(mx) }


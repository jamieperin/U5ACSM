#To fit cause specific under 5 Q, when qx's are predicted directly
#To match cause specific neonatal mortality
kcalc <- function(coefs, Q5, match) {
    ages<- c("0-6day","7-27day","1-5m","6-11m","12-23m","24-59m")
    vx <- coefs[, "vx"]
	if (length(Q5)==1) {
		mxhat1 <- mxhat.logquad(coefs, sex="Total", Q5, k=0)
		tmp.fcn <- function(k, match, vx, mx)  abs(-match + mx*exp(k*vx))
		tmp<- optimize(tmp.fcn, lower=-100, upper=100, match= match, mx=mxhat1[2], vx=vx[2], maximum=FALSE, tol=0.000000001)
		#print(tmp)
		k<- tmp$minimum

        }
	else {
		k <- Q5
		for (j in 1:length(Q5))  k[j] <- kcalc( coefs, Q5[j], match[j])
                }
	return(k) }

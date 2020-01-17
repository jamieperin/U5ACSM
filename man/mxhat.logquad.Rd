\name{mxhat.logquad}
\alias{mxhat.logquad}

\title{Estimating specific life tables of fixed child ages}
\usage{
mxhat.logquad(coefs, sex="Total", Q5, k)
}
\description{
	 Inputs:   Q5, sex, k (optional; default is k=0), and model coefficients
	           Q5 and k can be scalar or vector, but must have same length
	           coefs contains coefficients ax, bx, cx, and vx \\
	           
	 Outputs:  Vector (or matrix) of predicted qx (cummulative) values
	           with age groups 0-6 days, 7-27 days, 1-5 months, 6-11 months, 
	           12-23 months, and 24-59 months
	}
\examples{
data(camr)
camr<- camr[order(camr$year, camr$reg, camr$res, camr$age),] #By year, strata, and age group
rate<- "qc_nCH9"

#format mortality rates
qx.all<- array(camr[,rate], dim = c(6, length(camr[,rate])/6))

#Estimate log quadratic model with cause specific u5mr
qxc<- camr[which(camr$agegp=="24-59m"),c("reg","res","year",rate)]
bifit.f<-bifit( cbind(log(qxc[,rate]), log(qxc[,rate])^2), t(log(qx.all+0.000000001)), c=6)
beta.f<-bifit.f$coef
beta.f

#Determine error of quadratic fit
yhat.1<- cbind(1, log(qxc[,rate]), log(qxc[,rate])^2) %*% bifit.f$coef
head(yhat.1)
resid.1<- yhat.1 - t(log(qx.all+0.000000001))

#Determine variability around average pattern
svd.f<- svd(resid.1[,],1,1)
vx.f<- svd.f$v
vx.f[1]<-0 

#Estimate k for specific life table,
#given u5 pneumonia mortality of 0.005756 and 
#neonatal pneumonia mortality of 0.002924
ax.f <- beta.f[1,]; bx.f <- beta.f[2,]; cx.f <- beta.f[3,]
coefs<- data.frame(ax = ax.f, bx = bx.f, cx=cx.f, vx = vx.f)
coefs
kest <- kcalc(coefs, Q5 = 0.005756, match = 0.002924)
kest

#At typical pattern at U5 pneumonia mortality of 0.005756
mxhat.logquad(coefs, sex="Total", Q5 = 0.005756, k= 0)

#For specific life table with U5 pneumonia mortality of 0.005756
#and neonatal pneumonia mortality of 0.002924
mxhat.logquad(coefs, sex="Total", Q5 = 0.005756, k= kest)

}
\references{
Wilmoth J, Zureick S, Canudas-Romo V, Inoue M, Sawyer C. A flexible two-dimensional mortality model for use in indirect estimation. \emph{Population studies}. 2012 Mar 1;66(1):1-28.
}
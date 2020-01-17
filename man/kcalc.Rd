\name{kcalc}
\alias{kcalc}

\title{Estimating specific life tables of fixed child ages}
\usage{
kcalc(coefs, U5, match = NMR)
}
\description{
	To fit cause specific under 5 Q, when qxs are predicted directly
  To match cause specific neonatal mortality}
\examples{
data(camr)
camr<- camr[order(camr$year, camr$reg, camr$res, camr$age),] #By year, strata, and age group
rate<- "qc_nCH9"

#format mortality rates
qx.all<- array(camr[,rate], dim = c(6, length(camr[,rate])/6))

#Estimate log quadratic model with cause specific u5mr
qxc<- camr[which(camr$agegp=="24-59m"),c("reg","res","year",rate)]
bifit.f<-bifit( cbind(log(qxc[,rate]), log(qxc[,rate])^2), t(log(qx.all+0.00000001)), c=6)
beta.f<-bifit.f$coef
beta.f

#Determine error of quadratic fit
yhat.1<- cbind(1, log(qxc[,rate]), log(qxc[,rate])^2) %*% bifit.f$coef
resid.1<- yhat.1 - t(qx.all)

#Determine variability around average pattern
svd.f<- svd(resid.1[,],1,1)
vx.f<- svd.f$v
vx.f[1]<-0 

#Estimate k for specific life table,
#given u5 pneumonia mortality of 0.005756 and 
#neonatal pneumonia mortality of 0.002924
ax.f <- beta.f[1,]; bx.f <- beta.f[2,]; cx.f <- beta.f[3,]
coefs<- data.frame(ax = ax.f, bx = bx.f, cx=cx.f, vx = vx.f)
kest <- kcalc(coefs, Q5 = 0.005756, match = 0.002924)
kest
}
\references{
Wilmoth J, Zureick S, Canudas-Romo V, Inoue M, Sawyer C. A flexible two-dimensional mortality model for use in indirect estimation. \emph{Population studies}. 2012 Mar 1;66(1):1-28.
}
\name{bifit}
\alias{bifit}

\title{WLS}
\usage{
bifit(x.f, y.f, c=6)
}
\description{
	Uses lsfit() to compute WLS using Tukey's biweight function
	c is the tuning constant, typically around 6 to 9
	tol is for convergence of the biweight fit; tolerance is for the matrix decomposition in lsfit
	Like lsfit, input x matrix should not include a column of ones
	Like lsfit, input y can be matrix if there are multiple left-hand sides
}
\examples{
data(u5mr)
data(camr)
rate<- "qc_nCH9"
qx.all<- array(camr[,rate], dim = c(6, length(camr[,rate])/6))
bifit.f<-bifit( cbind(log(u5mr$u5mr), log(u5mr$u5mr)^2), t(qx.all), c=6)
bifit.f$coef

yhat.1<- cbind(1, log(u5mr$u5mr), log(u5mr$u5mr)^2) %*% bifit.f$coef
resid.1<- yhat.1 - t(qx.all)
}
\references{
Wilmoth J, Zureick S, Canudas-Romo V, Inoue M, Sawyer C. A flexible two-dimensional mortality model for use in indirect estimation. \emph{Population studies}. 2012 Mar 1;66(1):1-28.
}
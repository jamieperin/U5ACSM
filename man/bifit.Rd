\name{bifit}
\alias{bifit}

\title{WLS}
\usage{
bifit()
}
\description{
	Uses lsfit() to compute WLS using Tukey's biweight function
	c is the tuning constant, typically around 6 to 9
	tol is for convergence of the biweight fit; tolerance is for the matrix decomposition in lsfit
	Like lsfit, input x matrix should not include a column of ones
	Like lsfit, input y can be matrix if there are multiple left-hand sides
}
\examples{
bifit()
}
\references{
Wilmoth J, Zureick S, Canudas-Romo V, Inoue M, Sawyer C. A flexible two-dimensional mortality model for use in indirect estimation. \emph{Population studies}. 2012 Mar 1;66(1):1-28.
}
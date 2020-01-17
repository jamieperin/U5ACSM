\name{mxhat.logquad}
\alias{mxhat.logquad}

\title{Estimating specific life tables of fixed child ages}
\usage{
kcalc()
}
\description{
	 Inputs:   Q5, sex, k (optional; default is k=0), and model coefficients
	           Q5 and k can be scalar or vector, but must have same length
	           coefs contains coefficients ax, bx, cx, and vx \\
	           
	 Outputs:  Vector (or matrix) of predicted sex-specific mx values
	           with age groups 0, 1-4, 5-9, ..., 110+
	}
\examples{
mxhat.logquad()
}
\references{
Wilmoth J, Zureick S, Canudas-Romo V, Inoue M, Sawyer C. A flexible two-dimensional mortality model for use in indirect estimation. \emph{Population studies}. 2012 Mar 1;66(1):1-28.
}
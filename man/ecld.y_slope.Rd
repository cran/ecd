% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ecld-y-slope-method.R
\name{ecld.y_slope}
\alias{ecld.y_slope}
\alias{ecld.y_slope_trunc}
\title{Analytic solution for the slope of \eqn{y(x)} in lambda distribution}
\usage{
ecld.y_slope(object, x)

ecld.y_slope_trunc(object, t = 1)
}
\arguments{
\item{object}{an object of ecld class}

\item{x}{a vector of \eqn{x} values}

\item{t}{numeric, for MGF truncation}
}
\value{
numeric
}
\description{
Analytic solution for the slope of \eqn{y(x)} if available.
\eqn{ecld.y_slope_trunc} calculates the MGF truncation point where
\eqn{dy/dx+t=1}. SGED is supported.
}
\examples{
ld <- ecld(sigma=0.01*ecd.mp1)
x <- seq(-0.1, 0.1, by=0.01)
ecld.y_slope(ld,x)
ecld.y_slope_trunc(ld)
}
\author{
Stephen H-T. Lihn
}
\keyword{y_slope}

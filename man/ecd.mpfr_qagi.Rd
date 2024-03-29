% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ecd-mpfr-qagi-method.R
\name{ecd.mpfr_qagi}
\alias{ecd.mpfr_qagi}
\title{Utility to integrate mpfr with infinity via qagi}
\usage{
ecd.mpfr_qagi(
  object,
  f,
  lower,
  upper,
  ...,
  abs.tol = .Machine$double.eps^0.25,
  show.warning = TRUE
)
}
\arguments{
\item{object}{an object of ecd class}

\item{f}{an R function taking a numeric first argument and
returning a numeric vector of the same length.
Returning a non-finite element will generate an error.}

\item{lower}{numeric, the lower limit of integration. Can be infinite.}

\item{upper}{numeric, the upper limit of integration. Can be infinite.}

\item{...}{addtional arguments for \code{f}.}

\item{abs.tol}{numeric, the suggested absolute tolerance.}

\item{show.warning}{logical, to suppress warnings or not.}
}
\value{
The \code{integrate} object
}
\description{
This utility supplements \code{Rmpfr::integrateR} with the quadpack qagi method to handle
integration involving infinity. Qagi is a transformation of \eqn{x/sigma=(1-t)/t} for positive x, 
and \eqn{x/sigma=(t-1)/t} for negative x. \eqn{t=0} is represented by \code{.Machine$double.eps}.
This utility requires (a) \code{lower} or \code{upper} is \code{+/-Inf}; 
(b) \code{lower} and \code{upper} are of the same sign.
}
\author{
Stephen H. Lihn
}
\keyword{utility}

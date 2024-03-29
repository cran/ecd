% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ecdb-read-generic.R
\name{read.ecdb}
\alias{read.ecdb}
\alias{read}
\alias{read,ecdb-method}
\title{Read API for the ecdb}
\usage{
\method{read}{ecdb}(
  object,
  alpha = NULL,
  gamma = NULL,
  alim = NULL,
  glim = NULL,
  cusp = 0,
  polar_ext = FALSE
)

read(
  object,
  alpha = NULL,
  gamma = NULL,
  alim = NULL,
  glim = NULL,
  cusp = 0,
  polar_ext = FALSE
)

\S4method{read}{ecdb}(
  object,
  alpha = NULL,
  gamma = NULL,
  alim = NULL,
  glim = NULL,
  cusp = 0,
  polar_ext = FALSE
)
}
\arguments{
\item{object}{an object of ecdb class}

\item{alpha, gamma}{numeric vectors of points for cartesian product}

\item{alim, glim}{length-two numeric vectors of min and max range}

\item{cusp}{numeric. Type of cusp. Only 0 and 1 are allowed.
If cusp=1, read cusp data on the critical line.
Reading cusp data must be done from the \code{alpha} side.
Default: 0.}

\item{polar_ext}{logical, for polar coordinate extension: 
\code{R, theta, angle}.
Default: \code{FALSE}.}
}
\value{
The data.frame from ECDATTR table.
}
\description{
Read ecdb into data.frame. This can be accomplished by either specifying
the range of \code{alpha,gamma} or the cartesian product of \code{alpha,gamma}
point by point, or both. If both are specified, it follows a similar logic as plot 
how \code{x,y} is scoped by \code{xlim,ylim}.
}
\keyword{ecdb}

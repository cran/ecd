% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ecd-class.R
\docType{class}
\name{ecd-class}
\alias{ecd-class}
\title{The ecd class}
\description{
This S4 class is the major object class for elliptic lambda distribution.
It stores the ecd parameters, numerical constants that facilitates
quadpack integration, statistical attributes, and optionally,
an internal structure for the quantile function.
}
\section{Slots}{

\describe{
\item{\code{call}}{The match.call slot}

\item{\code{alpha,gamma,sigma,beta,mu}}{a length-one numeric.
These are core ecd parameters.}

\item{\code{cusp}}{a length-one numeric as cusp indicator. 
0: not a cusp; 
1: cusp specified by \code{alpha}; 
2: cusp specified by \code{gamma}.}

\item{\code{lambda}}{a length-one numeric, the leading exponent for the special model, default is 3.}

\item{\code{R,theta}}{a length-one numeric. These are derived ecd parameters in polar coordinate.}

\item{\code{use.mpfr}}{logical, internal flag indicating whether to use mpfr.}

\item{\code{const}}{A length-one numeric as the integral of \eqn{exp(y(x))} that normalizes the PDF.}

\item{\code{const_left_x}}{A length-one numeric marking the left point of PDF integration.}

\item{\code{const_right_x}}{A length-one numeric marking the right point of PDF integration.}

\item{\code{stats}}{A list of statistics, see \code{ecd.stats} for more information.}

\item{\code{quantile}}{An object of ecdq class, for quantile calculation.}

\item{\code{model}}{A vector of four strings representing internal classification:
\code{long_name.skew}, code{long_name},
\code{short_name.skew}, \code{short_name}.
This slot doesn't have formal use yet.}
}}

\section{Details}{

  The elliptic lambda distribution is the early research target of what becomes the stable 
  lambda distribution. It is inspired from elliptic curve, and is defined by a depressed cubic polynomial, 
  \deqn{
    -y(z)^3 - \left(\gamma+\beta z\right) y(z) + \alpha = z^2,
  }
  where \eqn{y(z)} must approach minus infinity as \eqn{z} approaches plus or minus infinity.
  The density function is defined as \deqn{
    P\left(x; \alpha, \gamma, \sigma, \beta, \mu\right)
    \equiv\, \frac{1}{C\,\sigma} 
    e^{y\left(\left|\frac{x-\mu}{\sigma}\right|\right)},
  }
  and \eqn{C} is the normalization constant,\deqn{
    C = \int_{-\infty}^{\infty}e^{y(z)}\,dz,
  }
  where 
  \eqn{\alpha} and \eqn{\gamma} are the shape parameters, 
  \eqn{\sigma} is the scale parameter, 
  \eqn{\beta} is the asymmetric parameter, 
  \eqn{\mu} is the location parameter.
}

\references{
For elliptic lambda distribution, see 
  Stephen Lihn (2015). 
  \emph{On the Elliptic Distribution and Its Application to Leptokurtic Financial Data}.
  SSRN: \url{https://www.ssrn.com/abstract=2697911}
}
\author{
Stephen H. Lihn
}
\keyword{class}
\keyword{constructor}

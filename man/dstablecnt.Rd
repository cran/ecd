% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lamp-stable-cnt-distribution-method.R
\name{dstablecnt}
\alias{dstablecnt}
\alias{pstablecnt}
\alias{rstablecnt}
\alias{qstablecnt}
\alias{cfstablecnt}
\alias{kstablecnt}
\title{Stable Count distribution}
\usage{
dstablecnt(x, alpha = NULL, nu0 = 0, theta = 1, lambda = NULL)

pstablecnt(x, alpha = NULL, nu0 = 0, theta = 1, lambda = NULL)

rstablecnt(n, alpha = NULL, nu0 = 0, theta = 1, lambda = NULL)

qstablecnt(q, alpha = NULL, nu0 = 0, theta = 1, lambda = NULL)

cfstablecnt(s, alpha = NULL, nu0 = 0, theta = 1, lambda = NULL)

kstablecnt(alpha = NULL, nu0 = 0, theta = 1, lambda = NULL)
}
\arguments{
\item{x}{numeric, vector of responses.}

\item{alpha}{numeric, the shape parameter, default is NULL. User must provide
either alpha or lambda.}

\item{nu0}{numeric, the location parameter, default is 0.}

\item{theta}{numeric, the scale parameter, default is 1.}

\item{lambda}{numeric, alternative shape parameter, default is NULL.}

\item{n}{numeric, number of observations.}

\item{q}{numeric, vector of quantiles.}

\item{s}{numeric, vector of responses for characteristic function.}
}
\value{
numeric, standard convention is followed:
        d* returns the density,
        p* returns the distribution function,
        q* returns the quantile function, and
        r* generates random deviates.
        The following are our extensions:
        k* returns the first 4 cumulants, skewness, and kurtosis,
        cf* returns the characteristic function.
}
\description{
Implements the stable count distribution
(based on stabledist package) for stable random walk simulation.
Quartic stable distribution is implemented through gamma distribution.
}
\section{Details}{

  The stable count distribution is the conjugate prior of the stable distribution.
  The density function is defined as \deqn{
    \mathit{N}_{\alpha}\left(\nu;\nu_{0},\theta\right) 
    \equiv\frac{\alpha}{\mathit{\Gamma}\left(\frac{1}{\alpha}\right)}\,
    \frac{1}{\nu-\nu_{0}}\,L_{\alpha}\left(\frac{\theta}{\nu-\nu_{0}}\right),
    \:\mathrm{where}\,\nu>\nu_{0}.
  }{
    N_\alpha(\nu; \nu_0, \theta) = \alpha/\Gamma(1/\alpha) * 
    1/(\nu-\nu_0) * L_\alpha(1/(\nu-\nu_0))
  }
  where \eqn{\nu>\nu_0}. \eqn{\alpha} is the stability index, 
  \eqn{\nu_0} is the location parameter, and \eqn{\theta} is the scale parameter.
  \cr
  At \eqn{\alpha=0.5} aka \eqn{\lambda=4}, it is called "quartic stable count distribution",
  which is a gamma distribution with shape of 3/2. It has the closed form of \deqn{
    \mathit{N}_{\frac{1}{2}}\left(\nu;\nu_{0},\theta\right)
    \equiv\frac{1}{4\sqrt{\pi}\,\theta^{3/2}}
    \left(\nu-\nu_{0}\right)^{\frac{1}{2}}
    e^{-\frac{\nu-\nu_{0}}{4\theta}}
  }{
    N_\alpha(\nu; \nu_0, \theta) = 1/(4 sqrt(\pi) \theta^1.5)
    (\nu-\nu_0)^0.5 exp(-(\nu-\nu_0)/(4\theta))
  }
}

\examples{
  # generate the pdf of the VIX distribution
  x <- c(0, 100, by=0.1)
  pdf <- dstablecnt(x, nu0=10.4, theta=1.6, lambda=4)
  
}
\references{
For more detail, see Section 2.4 and Section 3.3 of 
  Stephen Lihn (2017). 
  \emph{A Theory of Asset Return and Volatility 
  under Stable Law and Stable Lambda Distribution}.
  SSRN: 3046732, \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3046732}.
  This distribution is also documented formally 
  in Wikipedia: \url{https://en.wikipedia.org/wiki/Stable_count_distribution}.
}
\author{
Stephen H-T. Lihn
}
\keyword{Stable}

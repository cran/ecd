% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lamp-stdlap-distribution-method.R
\name{dstdlap}
\alias{dstdlap}
\alias{pstdlap}
\alias{qstdlap}
\alias{rstdlap}
\alias{cfstdlap}
\alias{kstdlap}
\alias{dstdlap_poly}
\title{Standardized Laplace process and distribution}
\usage{
dstdlap(x, t = 1, convo = 1, beta = 0, mu = 0)

pstdlap(x, t = 1, convo = 1, beta = 0, mu = 0)

qstdlap(q, t = 1, convo = 1, beta = 0, mu = 0)

rstdlap(n, t = 1, convo = 1, beta = 0, mu = 0)

cfstdlap(s, t = 1, convo = 1, beta = 0, mu = 0)

kstdlap(t = 1, convo = 1, beta = 0, mu = 0)

dstdlap_poly(x, t = 1, convo = 1, beta = 0, mu = 0)
}
\arguments{
\item{x}{numeric, vector of responses.}

\item{t}{numeric, the time parameter, of which the variance is t.}

\item{convo}{numeric, the convolution number, default is 1, which is Laplace without convolution.
There is a special provision in \code{rstdlap}, where it will simulate the Wiener process
if \code{convo=Inf} and \code{beta=0}.}

\item{beta}{numeric, skewness parameter according to skewed lambda distribution, default is 0.}

\item{mu}{numeric, location parameter, default is 0.}

\item{q}{numeric, vector of quantiles.}

\item{n}{numeric, number of observations.}

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
Implements the standardized Laplace process and distribution.
Be aware of the performance concerns:
(a) The cumulative density function is implemented by direct integration over the density.
(b) The quantile function is implemented by root finding on cumulative density function.
}
\section{Details}{

  The Lihn-Laplace distribution is the stationary distribution of Lihn-Laplace process.
  The density function is defined as \deqn{
    f_{\mathit{L}}^{\left(m\right)}\left(x;t,\beta,\mu\right)
    \equiv\frac{1}{\sqrt{\pi}\Gamma\left(m\right)\sigma_{m}}\,
    \left|\frac{x-\mu}{2B_{0}\sigma_{m}}\right|^{m-\frac{1}{2}}
    K_{m-\frac{1}{2}}\left(\left|\frac{B_{0}(x-\mu)}{\sigma_{m}}\right|\right)
    e^{\frac{\beta (x-\mu)}{2\sigma_{m}}}
  }{
    f_L^(m)(x;t,\beta,\mu) =
    1/(sqrt(\pi)\Gamma(m)\sigma_m)
    |(x-\mu)/2B0 \sigma_m|^(m-1/2)
    K_(m-1/2)(|B0 (x-\mu)/\sigma_m|)
    exp(\beta (x-\mu)/2\sigma_m)
  }
  where \deqn{
    \sigma_{m}\equiv\sqrt{\frac{t}{m\left(2+\beta^{2}\right)}},
    \:
    B_{0}\equiv\sqrt{1+\frac{1}{4}\beta^{2}}.
  }{
    \sigma_m = sqrt(t/m/(2+\beta^2)),
    B0 = sqrt(1+\beta^2/4).
  }
  \eqn{K_n(x)} is the modified Bessel function of the second kind.
  \eqn{t} is the time or sampling period,
  \eqn{\beta} is the asymmetric parameter,
  \eqn{\mu} is the location parameter.
}

\examples{
  # generate the pdf at time t=1 for the second convolution and beta = 0.1 for skewness
  x <- c(-10, 10, by=0.1)
  pdf <- dstdlap(x, t=1, convo=2, beta=0.1)

}
\references{
For more detail, see Section 5.4 of
  Stephen Lihn (2017).
  \emph{A Theory of Asset Return and Volatility
  under Stable Law and Stable Lambda Distribution}.
  SSRN: 3046732, \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3046732}.
}
\author{
Stephen H-T. Lihn
}
\keyword{Laplace}

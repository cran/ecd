% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ecd-cusp-constructor.R
\name{ecd.cusp}
\alias{ecd.cusp}
\title{Cusp constructor of ecd class}
\usage{
ecd.cusp(
  alpha = NaN,
  gamma = NaN,
  sigma = 1,
  mu = 0,
  with.stats = TRUE,
  with.quantile = FALSE,
  bare.bone = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{alpha}{numeric, the flatness parameter. Default: NaN.}

\item{gamma}{numeric, the sharpness parameter. Default: NaN.}

\item{sigma}{numeric, the scale parameter. Must be positive. Default 1.}

\item{mu}{numeric, the location parameter. Default: 0.}

\item{with.stats}{logical, also calculate statistics, default is \code{TRUE}.}

\item{with.quantile}{logical, also calculate quantile data, default is \code{FALSE}.}

\item{bare.bone}{logical, skip both const and stats calculation, default is \code{FALSE}.
This for debug purpose for issues on integrating \eqn{e^y(x)}.}

\item{verbose}{logical, display timing information, for debugging purpose, default is \code{FALSE}.}
}
\value{
The ecd class
}
\description{
Construct an ecd class for cusp distribution by specifying either alpha or gamma,
but not both. At the moment, it can't handle beta.
}
\examples{
d <- ecd.cusp(alpha=1)
d <- ecd.cusp(gamma=-1)
}
\author{
Stephen H. Lihn
}
\keyword{constructor}
\keyword{cusp}
\keyword{ecd}

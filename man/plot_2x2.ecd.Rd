% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ecd-plot-2x2-generic.R
\name{plot_2x2.ecd}
\alias{plot_2x2.ecd}
\alias{plot_2x2}
\alias{plot_2x2,ecd-method}
\title{Standard 2x2 plot for sample data}
\usage{
plot_2x2.ecd(object, ts, EPS = FALSE, eps_file = NA)

plot_2x2(object, ts, EPS = FALSE, eps_file = NA)

\S4method{plot_2x2}{ecd}(object, ts, EPS = FALSE, eps_file = NA)
}
\arguments{
\item{object}{An object of ecd class.}

\item{ts}{The xts object for the timeseries.}

\item{EPS}{Logical, indicating whether to save the plot to EPS, default = FALSE}

\item{eps_file}{File name for eps output}
}
\description{
Standard 2x2 plot for sample data
}
\examples{
\dontrun{
plot_2x2(d, ts)
}
}
\keyword{plot}

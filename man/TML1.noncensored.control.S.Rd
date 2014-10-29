\name{TML1.noncensored.control.S}
\alias{TML1.noncensored.control.S}
\title{
Control parameters for S-estimate of location and scale
}
\description{
      Auxiliary function for \code{\link{TML1.noncensored}}. Typically only used
      internally by \code{TML1.noncensored}, but may be used to construct a control argument.
      This function provides default values.
}
\usage{
TML1.noncensored.control.S(tlo = 1e-04, mxf = 50, mxs = 50, ntm = 50,
                           tls = 1e-06, h = 100)
}
\arguments{
  \item{tlo}{Relative tolerance in the iterative algorithms.}
  \item{mxf}{Maximum number of iterations in computing the location estimate.}
  \item{mxs}{Maximum number of iterations in computing the scale estimate.}
  \item{ntm}{Parameter used in iteration monitoring. When the number of iterations
  is a multiple of ntm, the current parameter values are printed.}
  \item{tls}{Tolerance for denominators. If a scale estimate is less than tls, the scale
  estimate is set equal to tls.}
  \item{h}{The number of subdivisions of the interval (min(yi), max(yi)) used is the
computation of the estimate \eqn{\lambda^(0)}.}
}

\value{
List containing the desired values for each of the control parameters, plus the value
Beta0 of \eqn{\beta}.
}

\keyword{misc }


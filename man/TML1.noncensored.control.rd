\name{TML1.noncensored.control}
\alias{TML1.noncensored.control}

\title{Control Parameters forTruncated Maximum Likelihood Estimation
       of Location and Scale }

\description{
      Auxiliary function for \code{\link{TML1.noncensored}}. Typically only used
      internally by \code{TML1.noncensored}, but may be used to construct a control argument.
      This function provides default values.
}

\usage{
TML1.noncensored.control(iv = 1, gam = 0.1, maxit = 200, tol = 1e-04)}

\arguments{
  \item{iv}{
     \itemize{
     \item 0: use and do not change the initial estimate of scale.
     \item 1: compute a truncated maximum likelihood estimate of scale.}}
  \item{gam}{ Relaxation factor for the IRLS algorithm for the final estimate. Set 0 < gam <= 1.}
  \item{maxit}{Maximum number of iterations in the IRLS algorithm for the final estimate. }
  \item{tol}{ Relative tolerance in the IRLS algorithm for the final estimate. }   }

\value{
A list with components named as the arguments.}

\seealso{ \code{\link{TML1.noncensored}}}


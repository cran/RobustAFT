\name{TML1.noncensored}
\alias{TML1.noncensored}
\title{ Truncated Maximum Likelihood Estimates of Location and Scale}
\description{
      This functions computes the truncated maximum likelihood estimates of location 
      and scale described in Marazzi and Yohai (2004). 
      It assumes that the error distribution is approximately Gaussian or log-Weibull. 
      The cut-off values for outlier rejection are fixed or adaptive.
      This function is a simplified version of \code{\link{TML.noncensored}} for the 
      case without covariates.
}
\usage{
TML1.noncensored(y, errors= c("Gaussian", "log-Weibull"), cu = NULL, 
     initial = c("S", "input"), otp = c("adaptive", "fixed"), 
     cov = c("no", "parametric", "nonparametric"), input = NULL, 
     control = list(), ...)
}

\arguments{
  \item{y}{ Observation vector }
  \item{errors}{
      \itemize{
      \item "Gaussian": the error distribution is assumed to be approximately Gaussian.
      \item "log-Weibull" : the error distribution is assumed to be approximately 
        log-Weibull. }}
  \item{cu}{Preliminary minimal upper cut-off. The default is 2.5 in the Gaussian case 
        and 1.855356 in the log-Weibull case. }
  \item{initial}{ 
      \itemize{
      \item "S"     : initial S-estimate.
      \item "input" : the initial estimate is given on input.}}
  \item{otp}{
      \itemize{
      \item "adaptive": adaptive cut-off.
      \item "fixed"   : non adaptive cut-off.}}

  \item{cov}{
      \itemize{
      \item "no": no estimate of the covariance matrix of the estimates is provided on 
        output.
      \item "parametric": a parametric estimate of the covariance matrix of the 
                    location-scale estimates is provided (to be used when n is small).
      \item "nonparametric": a nonparametric estimate of the covariance 
                    matrix of the location-scale estimates is provided.}}

  \item{input}{
      Initial input estimates of location and scale.\cr 
      Required when initial="input". 
      \itemize{
      \item "Gaussian case"   : list(theta=...,sigma=...) initial input estimates. theta: 
       location; sigma: scale.
      \item "log-Weibull case" : list(tau=...,v=...) initial input estimates of location 
       (tau) and scale (v).}}

  \item{control}{ Control parameters. For the default values, see the function 
       \code{\link{TML1.noncensored.control}.}}
  \item{\dots}{ If initial="S", parameters for the computation of the initial S estimates.
        See the function \code{\link{TML1.noncensored.control.S}} for the default values.}
}
\value{
  A list with the following components:
  \item{th0 }{Initial location estimate (S or input).}
  \item{v0 }{Initial scale estimate (S or input).}
  \item{nit0 }{Reached number of iteration if initial="S"}
  \item{th1 }{Final location estimate.}
  \item{v1 }{Final scale estimate.}
  \item{nit1 }{Reached iteration number in IRLS algorithm for final estimate (only for the 
    log-Weibull case).}
  \item{tu, tl }{Final cut-off values.}
  \item{alpha }{Estimated proportion of retained observations.}
  \item{tn }{Number of retained observations.}
  \item{beta }{Consistency constant for scale.}
  \item{wi }{Vector of weights (0 for rejected observations, 1 for retained observations).}
  \item{CV0 }{Covariance matrix of the initial estimates (th0,v0). }
  \item{CV1 }{Covariance matrix of the final   estimates (th1,v1). }
}

\references{ 
  Marazzi A., Yohai V. (2004). Adaptively truncated maximum likelihood regression with 
  asymmetric errors. \emph{Journal of Statistical Planning and Inference}, 122, 271-291.
 }

\seealso{\code{\link{TML.noncensored}}, \code{\link{TML1.noncensored.control}}, 
   \code{\link{TML1.noncensored.control.S}} }

\keyword{ robust }
\keyword{ univar }

\examples{
      data(Z243)
      Cost <- Z243$CouTot                         
      y    <- log(Cost)
      ctrl <- TML1.noncensored.control(iv=1,tol=1e-3)
 \dontrun{    
      z    <- TML1.noncensored(y,errors="log-Weibull", initial="S",otp="adaptive",
              cov="no",control=ctrl)
}
}


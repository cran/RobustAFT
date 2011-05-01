\name{TML.censored}
\alias{TML.censored}
\title{Truncated Maximum Likelihood Regression With Censored Observations}
\description{
      This function computes the truncated maximum likelihood estimates of accelerated failure time regression 
      described in Locatelli et al. (2010). 
      The error distribution is assumed to follow approximately a Gaussian or a log-Weibull distribution.
      The cut-off values for outlier rejection are fixed or adaptive.
}
\usage{
TML.censored(formula, delta, data, errors = "Gaussian",  initial = "S", 
             input = NULL, otp = "fixed", cu = NULL, control.S=list(), 
             control.ref=list(), control.tml=list())
}
\arguments{
  \item{formula}{ A \code{\link[stats]{formula}}, i.e., a symbolic description of the model
              to be adjusted (cf. \code{\link[stats]{glm}} or \code{\link[stats]{lm}}). }
  \item{data}{ An optional data frame containing the variables in the model. If not
              found in \code{data}, the variables are taken from \code{environment(formula)},
              typically the environment from which \code{robaft} is called. }
  \item{delta}{ Vector of 0 and 1. 
      \itemize{
          \item 0: censored observation.
          \item 1: complete observation.}}
  \item{errors}{
    \itemize{
        \item "Gaussian": the error distribution is assumed to be Gaussian.
        \item "logWeibull" : the error distribution is assumed to be log-Weibull.}}  
  \item{initial}{
      \itemize{
        \item "S": initial S-estimate.
        \item "input": the initial estimate is given on input.}}
  \item{input}{ A list(theta=c(...),sigma=...): initial input estimates where
      theta is a vector of p coefficients and sigma a scalar scale.\cr
      Required when initial="input".}
  \item{otp}{
      \itemize{
        \item "adaptive": adaptive cut-off.
        \item "fixed": non adaptive cut-off.}}
  \item{cu}{ Preliminary minimal upper cut-off. 
             The default is 2.5 in the Gaussian case and 1.855356 in the log-Weibull case. }  
  \item{control.S}{ A list of control parameters for the computation of the initial S estimates. 
           See the function \code{\link{TML.censored.control.S}} for the default values.}                                     
  \item{control.ref}{ A list of control parameters for the refinement algorithm of the initial S estimates. 
           See the function \code{\link{TML.censored.control.ref}} for the default values.}
  \item{control.tml}{ AA list of control parameters for the computation of the final estimates.
           See the function \code{\link{TML.censored.control.tml}} for the default values. }
}

\value{
  
  \code{TML.censored} returns an object of class "TML". 
  The function \code{\link[base]{summary}} can be used to obtain or print a summary of the results. 
  The generic extractor functions \code{\link[stats]{fitted}}, \code{\link[stats]{residuals}} and 
  \code{\link[stats]{weights}} can be used to extract various elements of the object returned 
  by \code{TML.censored}. The function \code{\link[stats]{update}} can be used to update the model.

  An object of class "TML" is a list with at least the following components: 
  \item{th0 }{Initial coefficient estimates.}
  \item{v0 }{Initial scale estimate.}
  \item{nit.ref }{Reached number of iteration in the refinement step for the initial estimates.}
  \item{th1 }{Final coefficient estimates.}
  \item{v1 }{Final scale estimate.}
  \item{nit.tml }{Number of iterations reached in IRLS algorithm for the final estimates.}
  \item{tu,tl }{Final cut-off values.}
  \item{alpha }{Estimated proportion of retained observations.}
  \item{tn }{Number of retained observations.}
  \item{weights }{Vector of weights (0 for rejected observations, 1 for retained observations).}
  \item{COV }{Covariance matrix of the final estimates (th1[1],...,th1[p],v1) (where p=ncol(X)).}
  \item{residuals }{ Residuals of noncensored observations are calculated as response minus fitted values. 
        For censored observations, the the expected residuals given that 
        the response is larger than the recorded censored value are provided.}
  \item{fitted.values }{The fitted mean values.}
  \item{call }{The matched call.}
  \item{formula }{The formula supplied.}
  \item{terms }{The \code{\link[stats]{terms}} object used.}
  \item{data }{The \code{data argument}.}
}
\references{
  Locatelli I., Marazzi A., Yohai V. (2010). Robust accelerated failure time regression.
  \emph{ Computational Statistics and Data Analysis}, 55, 874-887.
 }

\seealso{ \code{\link{TML.censored.control.ref}},
          \code{\link{TML.censored.control.tml}},
          \code{\link{TML.censored.control.S}},
          \code{\link{TML.noncensored}}
}

\keyword{ Regression }
\keyword{ Robust}
\keyword{ Accelerated Failure Time}

\examples{
     # This is the example described in Locatelli et al. (2010). 
     # The estimates are slighty different than those of the paper due to changes 
     # in the algorithm for the final estimate.
     #
     data(MCI)
     attach(MCI)
     
     # Exploratory Analysis
     plot(Age,log(LOS),type= "n",cex=0.7)

     # (1) filled square : regular,   complete
     # (2) empty  square : regular,   censored
     # (3) filled triangle : emergency, complete
     # (4) empty  triangle : emergency, censored

     points(Age[Dest==1 & TypAdm==0], log(LOS)[Dest==1 & TypAdm==0], pch=15,cex=0.7) # (1)
     points(Age[Dest==0 & TypAdm==0], log(LOS)[Dest==0 & TypAdm==0], pch=0, cex=0.7) # (2) 
     points(Age[Dest==1 & TypAdm==1], log(LOS)[Dest==1 & TypAdm==1], pch=17,cex=0.7) # (3) 
     points(Age[Dest==0 & TypAdm==1], log(LOS)[Dest==0 & TypAdm==1], pch=2, cex=0.7) # (4) 

     # Maximum Likelihood
     ML   <- survreg(Surv(log(LOS), Dest) ~ TypAdm*Age, dist="gaussian")
     summary(ML)
     B.ML <- ML$coef
     S.ML <- ML$scale
     
     abline(c(B.ML[1]        ,B.ML[3]        ),lwd=1,col="grey",lty=1)
     abline(c(B.ML[1]+B.ML[2],B.ML[3]+B.ML[4]),lwd=1,col="grey",lty=1)
     
     # Robust Accelerated Failure Time Regression with Gaussian errors
     ctrol.S   <- list(N=150, q=5, sigma0=1, MAXIT=100, TOL=0.001,seed=123)

     ctrol.ref <- list(maxit.sigma=2,tol.sigma=0.0001,maxit.Beta=2,tol.Beta=0.0001,
           Maxit.S=50, tol.S.sigma=0.001, tol.S.Beta=0.001,alg.sigma=1,nitmon=FALSE)

     ctrol.tml <- list(maxit.sigma=50,tol.sigma=0.0001,maxit.Beta=50,tol.Beta=0.0001,
      Maxit.TML=50, tol.TML.sigma=0.001, tol.TML.Beta=0.001, alg.sigma=1,nitmon=FALSE)
     
     WML<-TML.censored(log(LOS)~TypAdm*Age,data=MCI,delta=Dest,otp="adaptive",
          control.S=ctrol.S,control.ref=ctrol.ref,control.tml=ctrol.tml)

     summary(WML)
     
     B.WML<-coef(WML)
     abline(c(B.WML[1]         ,B.WML[3]         ),lty=1, col="red")
     abline(c(B.WML[1]+B.WML[2],B.WML[3]+B.WML[4]),lty=1, col="red")
}

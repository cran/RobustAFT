\name{TML.censored.control.tml}
\alias{TML.censored.control.tml}

\title{Control parameters for the IRLS algorithm of the final TML.censored estimates }

\description{
      Auxiliary function for \code{\link{TML.censored}}. 
      Typically only used internally by \code{TML.censored}, but may be used to provide
  	  a control argument. This function provides default values.}

\usage{
TML.censored.control.tml(maxit.sigma=20, tol.sigma=0.0001, maxit.Beta=20, 
    tol.Beta=0.0001,Maxit.TML=50, tol.TML.sigma=0.001, tol.TML.Beta=0.001, 
    alg.sigma=1, nitmon = FALSE)}

\arguments{
  \item{maxit.sigma}{ Maximum number of iterations in scale step.}
  \item{tol.sigma}{ Tolerance for sigma in scale step.}
  \item{maxit.Beta}{ Maximum number of iterations in coefficient step.}
  \item{tol.Beta}{ Tolerance for coefficients in coefficient step.}
  \item{Maxit.TML}{ Maximum number of iterations for global cycle.}
  \item{tol.TML.sigma}{ Tolerance for sigma in global cycle.}
  \item{tol.TML.Beta}{ Tolerance for coefficients in global cycle.}
  \item{alg.sigma}{ Type of algorithm in scale step:
    \itemize{ 
    \item 1: fixed point algorithm.
    \item 2: regula falsi.}}
  \item{nitmon}{ Set to TRUE if iteration monitoring is desired. Default=FALSE.}}

\value{
  A list with components named as the arguments.}

\seealso{ \code{\link{TML.censored}},
          \code{\link{TML.censored.control.S}}, 
          \code{\link{TML.censored.control.ref}}}

\examples{
     ### In the example(TML.censored), the control argument for the final estimates 
	 ### can be built using this function:
	 
	 \dontrun{
     data(MCI)
     attach(MCI)
     
     # Robust Accelerated Failure Time Regression with Gaussian errors
     ctrol.ref <- list(maxit.sigma=2,tol.sigma=0.0001,maxit.Beta=2,tol.Beta=0.0001,
           Maxit.S=50, tol.S.sigma=0.001, tol.S.Beta=0.001,alg.sigma=1,nitmon=FALSE)

     ctrol.tml <- TML.censored.control.tml(maxit.sigma=50,tol.sigma=0.0001,
           maxit.Beta=50,tol.Beta=0.0001, Maxit.TML=50, tol.TML.sigma=0.001, 
           tol.TML.Beta=0.001, alg.sigma=1,nitmon=FALSE)
     
     WML   <- TML.censored(log(LOS)~TypAdm*Age,data=MCI,delta=Dest,otp="adaptive",
           control.ref=ctrol.ref,control.tml=ctrol.tml)

     summary(WML)
}
}

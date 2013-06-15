\name{TML.censored.control.S}
\alias{TML.censored.control.S}

\title{Control parameters for the computation of the initial S estimates in TML.censored}

\description{
      Auxiliary function for \code{\link{TML.censored}}. 
      Typically only used internally by \code{TML.censored}, but may be used to provide a
      control argument. This function provides default values.}

\usage{
TML.censored.control.S(N=100, q=6, sigma0=1, MAXIT=100, TOL=0.01, seed=153)}

\arguments{
  \item{N}{ Number of subsamples.}
  \item{q}{ Subsample size.}
  \item{sigma0}{ Initial value of scale.}
  \item{MAXIT}{ Maximum number of iterations for solving the equation for scale.}
  \item{TOL}{ Relative tolerance for scale.}
  \item{seed}{ Seed for the random number generator.}  }

\value{
  A list with components named as the arguments.}

\seealso{ \code{\link{TML.censored}},
          \code{\link{TML.censored.control.ref}},
          \code{\link{TML.censored.control.tml}}}

\examples{
     ### In the example(TML.censored), the control argument for the refinement 
     ### algorithm can be built using this function:

     data(MCI)
     attach(MCI)
     
     # Robust Accelerated Failure Time Regression with Gaussian errors
     ctrol.S   <- list(N=150, q=5, sigma0=1, MAXIT=100, TOL=0.001,seed=123)

     ctrol.ref <- TML.censored.control.ref(maxit.sigma=2,tol.sigma=0.0001,
               maxit.Beta=2,tol.Beta=0.0001, Maxit.S=50, tol.S.sigma=0.001, 
               tol.S.Beta=0.001,alg.sigma=1,nitmon=FALSE)

     ctrol.tml <- list(maxit.sigma=50,tol.sigma=0.0001,maxit.Beta=50,
               tol.Beta=0.0001, Maxit.TML=50, tol.TML.sigma=0.001, 
               tol.TML.Beta=0.001, alg.sigma=1,nitmon=FALSE)
 \dontrun{        
     WML       <- TML.censored(log(LOS)~TypAdm*Age,data=MCI,delta=Dest,
               otp="adaptive",control.S=ctrol.S,control.ref=ctrol.ref,
               control.tml=ctrol.tml)

     summary(WML)
}
}

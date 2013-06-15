\name{TML.noncensored.control}
\alias{TML.noncensored.control}

\title{Control Parameters for Truncated Maximum Likelihood Regression 
       Without Censored Observations}

\description{
      Control parameters for \code{\link{TML.noncensored}}. 
      Typically only used internally by \code{TML.noncensored}, but may be used to 
      construct a control argument. This function provides default values.
}
\usage{
TML.noncensored.control(iv = 1, nrep = 0, gam = 0.1, nitmon = FALSE, 
                maxit = 200, tol = 1e-04, fastS = FALSE, seed=1313)
}
\arguments{

  \item{iv}{
     \itemize{
     \item 0: use and do not change the initial estimate of scale.
     \item 1: compute a truncated maximum likelihood estimate of scale.}}
  \item{nrep}{
      \itemize{
      \item Number of subsamples to be used in the computation of the S-estimate.
      \item 0: exhaustive sampling if the observation number is not too large.}}
  \item{gam}{ Relaxation factor for the IRLS algorithm of final estimate. 
       Set 0 < gam <= 1.}
  \item{nitmon}{ Set to TRUE if iteration monitoring in IRLS algorithm for the final 
       estimate is desired. Default=FALSE. }
  \item{maxit}{ Maximum number of iterations in IRLS algorithm for the final estimate. }
  \item{tol}{ Relative tolerance in IRLS algorithm. }
  \item{fastS}{
      \itemize{
      \item "TRUE"  : the initial S-estimate is computed using 
           \code{\link[robustbase]{lmrob.S}} from the \pkg{robustbase} package. 
           The control parameters are taken from \code{lmrob.control}.
      \item "FALSE" : the initial S-estimate is computed using \code{hysest} from the
                \pkg{robeth} package.}}
  \item{seed}{ Seed for the random number generator in the resampling algorithm for 
           the initial S-estimate. }
}

\value{
  A list with components named as the arguments.}

\seealso{ \code{\link{TML.noncensored}}}

\examples{
     ### In the example(TML.noncensored), the control argument can be built 
     ### using this function:

     data(D243)
     Cost <- D243$Cost                             # Cost (Swiss francs)
     LOS  <- D243$LOS                              # Length of stay (days)
     Adm  <- D243$Typadm; Adm <- (Adm==" Urg")*1   # Type of admission 
                                                   # (0=on notification, 1=Emergency)
     Ass  <- D243$Typass; Ass <- (Ass=="P"   )*1   # Type of insurance 
                                                   # (0=usual, 1=private)
     Age  <- D243$age                              # Age (years)
     Dst  <- D243$dest;   Dst <- (Dst=="DOMI")*1   # Destination 
                                                   # (1=Home, 0=another hospital)
     Sex  <- D243$Sexe;   Sex <- (Sex=="M"   )*1   # Sex (1=Male, 0=Female)

     # Truncated maximum likelihood regression with Gaussian errors

     ctrol <- TML.noncensored.control(iv=1, nrep=0, gam=0.2, fastS=TRUE, nitmon=FALSE)
     z     <- TML.noncensored(log(Cost)~log(LOS)+Adm+Ass+Age+Dst+Sex, otp="adaptive")
     summary(z)
}

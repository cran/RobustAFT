\name{TML.noncensored}
\alias{TML.noncensored}

\title{Truncated Maximum Likelihood Regression Without Censored Observations}

\description{
      This function computes the truncated maximum likelihood regression estimate
      described in Marazzi and Yohai (2004). 
      The error distribution is assumed to follow approximately a Gaussian or a log-Weibull distribution. 
      The cut-off values for outlier rejection are fixed or adaptive.
}

\usage{
TML.noncensored(formula, data, errors = "Gaussian", cu = NULL, 
                initial = "S",otp = "fixed", cov = "parametric", 
                input = NULL, control = list(), ...)}

\arguments{
  \item{formula}{ A \code{\link[stats]{formula}}, i.e., a symbolic description of the model
              to be fit (cf. \code{\link[stats]{glm}} or \code{\link[stats]{lm}}). }
  \item{data}{ An optional data frame containing the variables in the model. If not
              found in \code{data}, the variables are taken from \code{environment(formula)},
              typically the environment from which \code{TML.noncensored} is called. }
  \item{errors}{
      \itemize{
        \item "Gaussian": the error distribution is assumed to be Gaussian.
        \item "logWeibull" : the error distribution is assumed to be log-Weibull. }}
  \item{cu}{ Preliminary minimal upper cut-off. 
             The default is 2.5 in the Gaussian case and 1.855356 in the log-Weibull case. }
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
      \item "no": no estimate of the covariance matrix of the coefficients is provided on output.
      \item "parametric": a parametric estimate of the covariance matrix of the
                    coefficients is provided (to be used when n is small).
      \item "nonparametric": a nonparametric estimate of the covariance
                    matrix of the coefficients is provided.}}
  \item{input}{
      Initial input estimates of location and scale.\cr 
      Required when initial="input". 
      \itemize{
      \item "Gaussian case"   : list(theta=...,sigma=...) initial input estimates. theta: location; sigma: scale.
      \item "logWeibull case" : list(tau=...,v=...) initial input estimates of location (tau) and scale (v).}}

  \item{control}{ Control parameters. For the default values, see the function \code{\link{TML.noncensored.control}}.}
  \item{\dots}{ If fastS=TRUE, parameters for \code{lmrob.S}. See the function \code{\link[robustbase]{lmrob.control}}
                (from the robustbase package) for the default values. }
}

\value{ 
  \code{TML.noncensored} returns an object of class "TML". 
  The function \code{\link[base]{summary}} can be used to obtain or print a summary of the results. 
  The generic extractor functions \code{\link[stats]{fitted}}, \code{\link[stats]{residuals}} and 
  \code{\link[stats]{weights}} can be used to extract various elements of the value returned 
  by \code{TML.noncensored}. The function \code{\link[stats]{update}} can be used to update the model.

  An object of class "TML" is a list with the following components: 
  \item{th0 }{Initial coefficient estimates (S or input).}
  \item{v0 }{Initial scale (S or input).}
  \item{nit0 }{Reached number of iteration in \code{lmrob.S} (available only if fastS is TRUE).}
  \item{th1 }{Final coefficient estimates.}
  \item{v1 }{Final scale (S or input).}
  \item{nit1 }{Number of iterations reached by the IRLS algorithm for the final estimates.}
  \item{tu,tl }{Final cut-off values.}
  \item{alpha }{Estimated proportion of retained observations.}
  \item{tn }{Number of retained observations.}
  \item{beta }{Consistency constant for scale.}
  \item{weights }{Vector of weights (0 for rejected observations, 1 for retained observations).}
  \item{COV}{Covariance matrix of the final estimates (th1[1],...,th1[p],v1) (where p=ncol(X)).}
  \item{residuals }{The residuals, that is response minus fitted values.}
  \item{fitted.values }{The fitted mean values.}
  \item{call }{The matched call.}
  \item{formula }{The formula supplied.}
  \item{terms }{The \code{\link[stats]{terms}} object used.}
  \item{data }{The \code{data argument}.}}

\references{
  Marazzi A., Yohai V. (2004). Adaptively truncated maximum likelihood regression with asymmetric errors.
  \emph{Journal of Statistical Planning and Inference}, 122, 271-291.
}

\seealso{ \code{\link{TML.noncensored.control}},
          \code{\link{TML1.noncensored}}, \code{\link{TML1.noncensored.control}}, 
          \code{\link{TML.censored}}
}

\keyword{ Regression }
\keyword{ Robust}
\keyword{ Accelerated Failure Time}

\examples{
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

     z    <- TML.noncensored(log(Cost)~log(LOS)+Adm+Ass+Age+Dst+Sex,
              otp="adaptive",control=list(fastS=TRUE))

     summary(z)
     
     # Truncated maximum likelihood regression with log-Weibull errors

     w    <- TML.noncensored(log(Cost)~log(LOS)+Adm+Ass+Age+Dst+Sex,
             errors="logWeibull",otp="adaptive",control=list(fastS=TRUE))

     summary(w)

}


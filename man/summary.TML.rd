\name{summary.TML}
\alias{summary.TML}
\alias{print.TML}
\alias{print.summary.TML}
\alias{coef.TML}
\alias{vcov.TML}

\title{Summarizing Truncated Maximum Likelihood regression }

\description{
      Summary and print \code{\link[utils]{methods}} for \code{R} object of class "TML" and \code{\link[base]{print}} 
      method for the summary object.
      Further, methods \code{\link[stats]{fitted}}(), \code{\link[stats]{residuals}}(), \code{\link[stats]{weights}}() 
      or \code{\link[stats]{update}}() work (via the default methods), and \code{\link[stats]{coef}}(), \code{\link[stats]{vcov}}() 
      have explicitly defined TML methods. }

\usage{
\method{summary}{TML}(object, ...)
\method{print}{TML}(x, digits = max(3, getOption("digits") - 3), ...)
\method{coef}{TML}(object, ...)
\method{vcov}{TML}(object, ...)

\method{print}{summary.TML}(x, digits = max(3, getOption("digits") - 3),
  signif.stars = getOption("show.signif.stars"), ...)
}

\arguments{
  \item{object}{ An object of class "TML", usually, a result of a call to \code{\link{TML.noncensored}} or \code{\link{TML.censored}}. }
  \item{\dots}{ Potentially more arguments passed to methods.}
  \item{digits}{ Number of digits for printing, see \code{digits} in \code{\link[base]{options}}.}
  \item{x}{ An object of class "TML" or "summary.TML". }
  \item{signif.stars}{ Logical indicating if the P-values should be visualized by so called "significance stars". }
}

\details{
      \code{summary.TML} returns an object of \code{\link[base]{class}} "summary.TML".
       
      \code{print.TML} returns a printed summary of object of class "TML".
      
      \code{print.summary.TML} tries to be smart about formatting the coefficients, standard errors, etc, and gives "significance stars" if signif.stars is TRUE (as per default when \code{\link[base]{options}} where not changed).
      
      \code{coef.TML} returns the final coefficient estimates (value \code{th1} of a "TML" object), and \code{vcov.TML} returns the covariance matrix of the final estimates (value \code{CV1} of a "TML" object).
}
\value{          
  An object of class "summary.TML" is a list with the following components:

  \item{call }{The component from \code{object}.}
  \item{terms }{The component from \code{object}.}
  \item{residuals }{The component from \code{object}.}
  \item{fitted.values }{The component from \code{object}.}
  \item{tn }{The component from \code{object}.}
  \item{coefficients }{The matrix of coefficients, standard errors, t-values and p-values. Aliased coefficients are omitted.}
  \item{aliased }{Named logical vector showing if the original coefficients are aliased.}
  \item{df }{Degrees of freedom, a 3-vector (p, n-p, p*), the last being the number of non-aliased coefficients.} 
  \item{sigma }{The final scale estimate from \code{object}.}
  \item{cutoff.values }{A vector of the final lower and upper cut-off values from \code{object}.}
}
\seealso{ \code{\link{TML.noncensored}}, \code{\link{TML.censored}}, \code{\link[base]{summary}}, \code{\link[base]{print}}
}
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
     z    <- TML.noncensored(log(Cost)~log(LOS)+Adm+Ass+Age+Dst+Sex, otp="adaptive", 
             cov="nonparametric", control=list(fastS=TRUE))

     z                  # -> print.TML(....)
     sumz <- summary(z) # -> summary.TML(....)
     sumz               # -> print.summary.TML(....)
     coef(z)            # -> coef.TML(....)
     vcov(z)            # -> vcov.TML(....)
}

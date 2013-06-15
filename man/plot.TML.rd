\name{plot.TML}
\alias{plot.TML}

\title{Plot Method for "TML" objects }

\description{
      Diagnostic plots for elements of class "TML". Three plots (selectable by which) 
      are currently available:
      a residual Q-Q plot, a plot of response against fitted values and a plot of 
      standardized residuals against fitted values.
}

\usage{
\method{plot}{TML}(x, which = 1:3, caption = c("Residual QQ-plot",
  "Response vs. Fitted Values", "Standardized Residuals vs. Fitted Values"),
  panel = points, sub.caption = deparse(x$call$formula), main = "",
  ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...)
}

\arguments{
  \item{x}{ An object of class "TML", usually, a result of a call to 
            \code{\link{TML.noncensored}} or \code{\link{TML.censored}}.}
  \item{which}{ If a subset of the plots is required, specify a subset of the numbers 
            \code{1:3}.}
  \item{caption}{ Caption for the different plots.}
  \item{panel}{ Panel.}
  \item{sub.caption}{ Sub titles.}
  \item{main}{ Main title.}
  \item{ask}{ If ask=TRUE, plot.TML() operates in interactive mode.}
  \item{\dots}{ Optional arguments for \code{\link[graphics]{par}}.}
}
\details{
      The residual Q-Q plot is build with respect to the \code{errors} argument of the 
       object. 
      This means that the expected order statistics are calculated either for a Gaussian 
      or a log-Weibull distribution.      
      The two horizontal dotted lines on the first and the third plots represent the upper 
      and lower cut-off values for outlier rejection.      
      Observations that were not retained for the estimation (outliers) are identified on 
      the third plot.
}

\seealso{ \code{\link{TML.noncensored}}, \code{\link{TML.censored}}, 
          \code{\link[graphics]{plot.default}}}

\keyword{ plot }

\examples{
     data(D243)
     Cost <- D243$Cost                             # Cost (Swiss francs)
     LOS  <- D243$LOS                              # Length of stay (days)
     Adm  <- D243$Typadm; Adm <- (Adm==" Urg")*1   # Type of admission 
                                                   # (0=on notification, 1=Emergency)
\dontrun{     
     # Truncated maximum likelihood regression with log-Weibull errors
     w  <- TML.noncensored(log(Cost)~log(LOS)+Adm, errors="log-Weibull", 
           otp="adaptive", control=list(fastS=TRUE))
     plot(w)
     plot(w, which = 1)
     plot(w, which = 2)
     plot(w, which = 3)
}
}

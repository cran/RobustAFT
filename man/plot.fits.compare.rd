\name{plot.fits.compare}
\alias{plot.fits.compare}

\title{Plot Method for "fits.compare" objects}

\description{
      Comparative plots for objects of class "fits.compare".
}

\usage{
\method{plot}{fits.compare}(x, xplots = FALSE, ask = TRUE, which = 1:4, 
             leg.position = c("topleft", "topleft", "topleft"), ...)
}

\arguments{
  \item{x}{ An object of class "fits.compare", usually, a result of a call to \code{\link{fits.compare}}.}
  \item{xplots}{ If xplots=TRUE, plots of the independent variables versus the residuals are produced. }
  \item{ask}{ If ask=TRUE, plot.fits.compare() operates in interactive mode. }
  \item{which}{ If a subset of the plots is required, specify a subset of the numbers \code{1:4}.}
  \item{leg.position}{ A vector of character string specifying the legend position of the second, third and fourth plots.}
  \item{\dots}{ Optional arguments for \code{\link[graphics]{par}}.}
}

\details{
    For clarity reasons, at most three models should be compared.
    Four default plots (selectable by which) are produced: histograms of the residuals of each model, 
    a residual Q-Q plot, response against fitted values and residuals against fitted values.
    Additional plots are produced if \code{xplots=TRUE}.
}

\seealso{ \code{\link{fits.compare}}, \code{\link[graphics]{plot.default}}, \code{\link{plot.TML}}
}

\examples{
\dontrun{
     data(D243)
     Cost <- D243$Cost                             # Cost (Swiss francs)
     LOS  <- D243$LOS                              # Length of stay (days)
     Adm  <- D243$Typadm; Adm <- (Adm==" Urg")*1   # Type of admission 
                                                   # (0=on notification, 1=Emergency)

     lwrob <- TML.noncensored(log(Cost)~log(LOS)+Adm, errors="logWeibull")
     reg   <- lm(log(Cost)~log(LOS)+Adm)

     comp  <- fits.compare(least.squares=reg, TML.logWeibull=lwrob)
     plot(comp, leg.position=c("topleft", "topleft", "bottomleft"), xplots=TRUE)
}
}

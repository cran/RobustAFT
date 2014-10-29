\name{fits.compare}
\alias{fits.compare}
\alias{print.fits.compare}
\title{Numerical comparison of several fits }
\description{
      Creates a class "fits.compare" object allowing the user to display model summary 
      statistics in a form allowing easy comparison of models.
}
\usage{
fits.compare(...)
}
\arguments{
  \item{\dots}{ one or more class "lm", class "lm.robust" or class "TML" objects. Names given to objects in the list are used as labeling information in the printed output.}
}
\details{
      The \code{fits.compare} function processes its arguments one at a time to create a named list of objects. 
      The object returned is a member of class "fits.compare". 
      Because of differences in the computed summary statistics, the list of input objects is currently limited to class "lm", 
      class "lm.robust" and class "TML" objects.
      The \code{print.fits.compare} function displays a textual comparison of the input models, 
      and the \code{plot.fits.compare} function provides comparative plots.
}
\value{
An object of class "fits.compare" containing the list of input models to be compared.
}
\seealso{ \code{\link{TML.noncensored}}, \code{\link{plot.fits.compare}}
}
\examples{

\dontrun{
     data(D243)
     Cost <- D243$Cost                             # Cost (Swiss francs)
     LOS  <- D243$LOS                              # Length of stay (days)
     Adm  <- D243$Typadm; Adm <- (Adm==" Urg")*1   # Type of admission 
                                                   # (0=on notification, 1=Emergency)

     lwrob <- TML.noncensored(log(Cost)~log(LOS)+Adm, errors="logWeibull")
     grob  <- TML.noncensored(log(Cost)~log(LOS)+Adm)
     reg   <- lm(log(Cost)~log(LOS)+Adm)

     fits.compare(least.squares=reg, TML.logWeibull=lwrob, TML.Gaussian=grob)
}
}

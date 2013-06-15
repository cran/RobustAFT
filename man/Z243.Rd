\name{Z243}
\alias{Z243}
\docType{data}
\title{ Sample of 100 hospital stays for medical back problems}
\description{
  Sample of 100 patients hospitalized for medical back problems in Switzerland
}
\usage{data(Z243)}
\format{
  A data frame with 100 observations on the following 14 variables.
  \describe{
    \item{\code{NoAdm}}{Admission number}
    \item{\code{APDRG}}{DRG: Always 243}
    \item{\code{Sex}}{Gender: 1=Male, 0=Female}
    \item{\code{Age}}{Age in years}
    \item{\code{LOS}}{Length of stay (days)}
    \item{\code{CouTot}}{Total Cost (Swiss francs)}
    \item{\code{CsansInv}}{Cost (Swiss francs)}
    \item{\code{Adm}}{Type of admission (0=on notification, 1=Emergency)}
    \item{\code{Ass}}{Type of insurance  (0=usual, 1=private)}
    \item{\code{Death}}{0=No, 1=Yes}
    \item{\code{BBD}}{A numeric vector}
    \item{\code{BBDaggr}}{A numeric vector}
    \item{\code{Charls}}{A numeric vector}
    \item{\code{LOSF}}{Adjusted length of stay}
  }
}

\examples{
data(Z243)
}
\keyword{datasets}

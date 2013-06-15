\name{D243}
\alias{D243}
\docType{data}
\title{ Sample of 100 hospital stays for medical back problems}
\description{
  Sample of 100 patients hospitalized for medical back problems in Switzerland
}
\usage{data(D243)}
\format{
  A data frame with 100 observations on the following 11 variables.
  \describe{
    \item{\code{Sexe}}{Gender: M=Male, F=Female}
    \item{\code{age}}{Age in years}
    \item{\code{dest}}{Destination: DOMI=Home else=another hospital}
    \item{\code{Typadm}}{Type of admission: Urg=Emergency else=on notification}
    \item{\code{Typass}}{Type of insurance: P=Private else=usual}
    \item{\code{LOS}}{Length of stay (days)}
    \item{\code{APDRG}}{DRG code: Always 243}
    \item{\code{Cost}}{Cost (Swiss francs)}
    \item{\code{CSansInv}}{Intermediate cost}
    \item{\code{BBDaggr}}{a numeric vector}
    \item{\code{BBD}}{a numeric vector}
  }
}

\examples{
data(D243)
}
\keyword{datasets}

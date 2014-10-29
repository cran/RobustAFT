\name{predict.TML}
\alias{predict.TML}

\title{Predict method for "TML" objects }

\description{
      Obtains predictions from a fitted Truncated Maximum Likelihood (TML) object.
}

\usage{
\method{predict}{TML}(object, newdata = NULL, ...)}

\arguments{
  \item{object}{ An object of class "TML", usually, a result of a call to \code{\link{TML.noncensored}} or \code{\link{TML.censored}}. }
  \item{newdata}{ Optionally, a vector, a matrix or a data frame containing the variables with which to predict. 
                  If omitted, the fitted values of \code{object} are returned.}
  \item{\dots}{ Additional arguments affecting the predictions produced.}
}

\details{
      \code{newdata} must have the same number of variables (that is of columns) as the model.
      If \code{object} is a model with an intercept, \code{newdata} must have a first column of 1.
}

\value{
  Returns a vector of predictions.
}

\seealso{ \code{\link{TML.noncensored}}, \code{\link{TML.censored}}, \code{\link[stats]{predict}}
}

\examples{
\dontrun{
     data(D243)
     Cost <- D243$Cost                             # Cost (Swiss francs)
     LOS  <- D243$LOS                              # Length of stay (days)
     Adm  <- D243$Typadm; Adm <- (Adm==" Urg")*1   # Type of admission 
                                                   # (0=on notification, 1=Emergency)

     # Fitting the model
     z    <- TML.noncensored(log(Cost)~log(LOS)+Adm, errors="logWeibull")

     # With a vector of data
     vec  <- c(1, 2.4, 1)
     predict(object = z, newdata = vec)
     # With a matrix of data
     mat  <- matrix(c(1,1,2.4,2.7,1,0), ncol=3)
     predict(z, mat)
     # With a data frame
     dat  <- as.data.frame(cbind("intercept"=c(1,1,1), "log(LOS)"=c(2.4,2.7,2.2), 
             "Adm"=c(1,0,1)))
     predict(z, dat)
}
}


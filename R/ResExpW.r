ResExpW <- function(r)
{
  nr <- length(r)
  num <- rep(0,nr)
  for (i in 1:nr) num[i] <- intxlwbl(r[i])
  den <- 1-plweibul(r)
  ind <- (den < 1e-5)
  val <- num/den
  val[ind] <- r[ind]
  val
}
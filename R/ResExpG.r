ResExpG <- function(r)
{
  num <- dnorm(r)
  den <- 1-pnorm(r)
  val <- num/den
  ind <- (den < 1e-5)
  val[ind] <- r[ind]
  val
}
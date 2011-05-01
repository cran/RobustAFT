TMLjac11.G <- function(d.Beta,d.sigma,rs0,delta,X,cu)
{
# Jacobian of Beta-equation wrt Beta.hat
  n     <- length(rs0); cl <- -cu
  D1 <- D2 <- D3 <- D <- rep(0,n)
  rsd   <- (rs0-X%*%d.Beta)/d.sigma
  Fo    <- pnorm(rsd)
  fo    <- dnorm(rsd)
  ai    <- (pmax(rs0,cl)-X%*%d.Beta)/d.sigma
  bi    <- (cu  -        X%*%d.Beta)/d.sigma
  foai  <- dnorm(ai)
  fobi  <- dnorm(bi)
  fopai <- foai*(-ai)
  fopbi <- fobi*(-bi)
  D1    <- - delta*ww(rs0,cl,cu)
  D2    <- - (1-delta)*fo/(1-Fo)^2*(foai-fobi )
  D3    <- - (1-delta)/(1-Fo)*( fopai-fopbi )
  D     <- D1 + D2 + D3
  Jac   <- t(X) %*% (as.vector(D)*X) /d.sigma/n
  Jac
}
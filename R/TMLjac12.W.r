TMLjac12.W <- function(d.Beta,d.sigma,rs0,delta,X,cl,cu)
{
# Jacobian of Beta-equation wrt sigma.hat
  n     <- length(rs0)
  D1 <- D2 <- D3 <- D <- rep(0,n)
  rsd   <- (rs0-X%*%d.Beta)/d.sigma
  Fo    <- plweibul(rsd)
  fo    <- dlweibul(rsd)
  ai    <- (pmax(rs0,cl)-X%*%d.Beta)/d.sigma
  bi    <- (cu  -        X%*%d.Beta)/d.sigma
  foai  <- dlweibul(ai)
  fobi  <- dlweibul(bi)
  fopai <- foai*(-ps0W(ai))
  fopbi <- fobi*(-ps0W(bi))
  D1    <- -  delta*ww(rs0,cl,cu)*psp0W(rsd)*rsd
  D2    <- - (1-delta)*fo/(1-Fo)^2*(foai-fobi)*rsd
  D3    <- - (1-delta)/(1-Fo)*(fopai*ai-fopbi*bi )
  D     <- D1 + D2 + D3
  Jac   <- t(X)%*%(as.vector(D))/d.sigma/n
  Jac
}
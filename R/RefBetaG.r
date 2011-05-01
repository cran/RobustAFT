RefBetaG <- function(X,y,delta,Beta,sigma,maxit,tol,nitmon)
{
# iteratively rewighting algorithm for Beta refinement, Gaussian case
  p <- length(Beta); n <- length(y); nu <- sum(delta); nc <- n-nu
  nit <- 1; Beta1 <- rep(100,p)
  I1 <- I0 <- I2 <- rep(0,nc)
  D1 <- D2 <- ym <- rep(0,n)
  indu <- (1:n)[delta==1]
  indc <- (1:n)[delta==0]
  while ( max(abs(Beta1-Beta)) > tol & (nit < maxit) ) {
    nit <- nit+1; Beta1 <- Beta
    mui <- X %*% as.matrix(Beta1)
    rs  <- (y-mui)/sigma
    if (nu > 0) {
      ru      <- rs[delta==1]
      wi      <- PspSG(ru)
      cnd     <- ru != 0
      wi[cnd] <- (PsiSG(ru[cnd])/ru[cnd])
    }
    if (nu < n) {
      rc <-  rs[delta==0]
      for (i in 1:nc) {
        I0[i] <- integrate(intg0,lower=rc[i],upper=6)$val
        I2[i] <- integrate(intg2,lower=rc[i],upper=6)$val
      }
      I1 <- sigma*I2+mui[delta==0]*I0
      vi <- I1/( 1-pnorm(rc) )
      ui <- I0/( 1-pnorm(rc) )
    }
    D1[indu] <- wi
    D1[indc] <- vi
    D2[indu] <- wi
    D2[indc] <- ui
    ym[indu] <- y[indu]
    ym[indc] <- rep(1,nc)
    A <- t(X)%*%(D1*ym)
    B <- t(X)%*%(D2*X)
    Beta <- solve(B)%*%A
    if(nitmon) cat(nit, Beta, Beta1, "\n")
  }
  list(Beta=Beta,nit=nit)
}
QMatrix.G <- function(d.Beta,d.sigma,rs0,delta,X,cu,const)
{
  cl <- -cu; p <- ncol(X);  n  <- length(rs0)
  Hu <- Hc <- matrix(0,nrow=(p+1),ncol=(p+1))
  ku <- kc <- 0
  iu  <- delta==1; ic <- delta==0;    nc <- sum(ic); nu <- n-nc
  if(nc < n) {
    Xu  <- X[iu,]; r0u <- rs0[iu]
    ru  <- (r0u-Xu%*%d.Beta)/d.sigma
    hi  <- ww(r0u,cl,cu)*ru
    ki  <- ww(r0u,cl,cu)*ru^2 - const
    Hu[1:p,1:p]      <- t(Xu)%*% (as.vector(hi^2) * Xu)
    Hu[(p+1),(p+1)]  <- sum(ki^2)
    Hu[1:p,(p+1)]    <- t(Xu)%*%as.matrix(hi*ki)
    Hu[(p+1),1:p]    <- Hu[1:p,(p+1)]
  }
  if(nc > 0) {
    Xc  <- X[ic,]; r0c <- rs0[ic]
    rc  <- (r0c-Xc%*%d.Beta)/d.sigma
    Fo  <- pnorm(rc)
    ai  <- (pmax(r0c,cl)- Xc%*%d.Beta)/d.sigma
    bi  <- (         cu - Xc%*%d.Beta)/d.sigma
    Ihi <- dnorm(ai)-dnorm(bi)
    Iki <- ai * dnorm(ai) - bi * dnorm(bi) + pnorm(bi) - pnorm(ai)
    ehi <- eki <- rep(0,nc)
    ehi[bi>ai]  <- Ihi[bi>ai]/(1-Fo[bi>ai])
    eki[bi>ai]  <- Iki[bi>ai]/(1-Fo[bi>ai]) - const
    Hc[1:p,1:p]      <- t(Xc)%*% (as.vector(ehi^2) * Xc)
    Hc[(p+1),(p+1)]  <- sum(eki^2)                              # check
    Hc[1:p,(p+1)]    <- t(Xc)%*%as.matrix(ehi*eki)              # check
    Hc[(p+1),1:p]    <- Hc[1:p,(p+1)]
  }
  (Hu+Hc)/n
}
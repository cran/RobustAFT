"invM2.w" <-
function(l,u,theta,sigma,rs,wi,XtX,xbar,estim=c("SA","TMLA","SI","TMLI")) {
n <- length(rs); p <- ncol(XtX); xk <- 1.717817
if (estim=="SA")  {zpsp <- RobustAFT:::psp.weight(rs,ips=2,xk=1.717817)
                   a1c  <- sum(zpsp)/(n*sigma)
                   b1c  <- sum(zpsp*rs)/(n*sigma)
                   zpsi <- RobustAFT:::psi.weight(rs,ips=2,xk=1.717817)
                   a2c  <- sum(zpsi)/(n*sigma)
                   b2c  <- sum(zpsi*rs)/(n*sigma)}
if (estim=="TMLA"){sp1 <- exp(rs); sp2 <- (exp(rs)*(1+rs)-1)
                   a1c <- sum(wi*sp1)/(n*sigma)
                   b1c <- sum(wi*sp1*rs)/(n*sigma)
                   a2c <- sum(wi*sp2)/(n*sigma)
                   b2c <- sum(wi*sp2*rs)/(n*sigma)}
if (estim=="SI")  {a1c <- integrate(RobustAFT:::Pspphi.w, lower=-xk,upper=xk)$value/sigma
                   a2c <- integrate(RobustAFT:::Psiphi.w, lower=-xk,upper=xk)$value/sigma
                   b1c <- integrate(RobustAFT:::Pspzphi.w,lower=-xk,upper=xk)$value/sigma
                   b2c <- integrate(RobustAFT:::Psizphi.w,lower=-xk,upper=xk)$value/sigma}
if (estim=="TMLI"){a1c <- integrate(RobustAFT:::s1pphi.w, lower=l,  upper=u)$value/sigma
                   a2c <- integrate(RobustAFT:::s2pphi.w, lower=l,  upper=u)$value/sigma
                   b1c <- integrate(RobustAFT:::s1pzphi.w,lower=l,  upper=u)$value/sigma
                   b2c <- integrate(RobustAFT:::s2pzphi.w,lower=l,  upper=u)$value/sigma}
A1 <- a1c*XtX;    b1 <- b1c*xbar
a2 <- a2c*xbar;   b2 <- b2c
M  <- matrix(0,ncol=p+1,nrow=p+1)
M[1:p,1:p] <- A1; M[1:p,p+1] <- as.matrix(b1)
M[p+1,1:p] <- a2; M[p+1,p+1] <- b2
Minv <- solve(M)
list(Minv=Minv,XtX=XtX,xbar=xbar)}


\name{RobustAFT-package}
\alias{RobustAFT-package}
\alias{RobustAFT}
\docType{package}
\title{
Robust Accelerated Failure Time Model Fitting
}
\description{
This package computes the truncated maximum likelihood regression estimates
described in Marazzi and Yohai (2004) and Locatelli et al. (2010). 
The error distribution is assumed to follow approximately a Gaussian or a log-Weibull distribution. 
The cut-off values for outlier rejection are fixed or adaptive.
The main functions of this package are \code{\link{TML.noncensored}} and \code{\link{TML.censored}}.
}
\details{
\tabular{ll}{
Package: \tab RobustAFT\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-05-01\cr
License: \tab GPL-2 or later\cr
LazyLoad: \tab yes\cr
}
}
\author{
Alfio Marazzi <Alfio.Marazzi@chuv.ch>, Jean-Luc Muralti

Maintainer: A. Randriamiharisoa <Alex.Randriamiharisoa@chuv.ch>
}
\references{
  Marazzi A., Yohai V. (2004). Adaptively truncated maximum likelihood regression with asymmetric 
  errors. \emph{Journal of Statistical Planning and Inference}, 122, 271-291.

  Locatelli I., Marazzi A., Yohai V. (2010). Robust accelerated failure time regression.
  \emph{Computational Statistics and Data Analysis}, 55, 874-887.

  Marazzi A. (1993) \emph{Algorithm, Routines, and S functions for Robust Statistics}.  
  Wadsworth & Brooks/cole, Pacific Grove, California.
}
\keyword{ regression }
\keyword{ robust}
\keyword{ package }
\examples{
# Example 1. This is the example described in Marazzi and Yohai (2004).
# ---------
# The two following auxiliary functions, not included in the library, 
#must be loaded.
# ------------- Auxiliary functions -------------------------------------

SDmux.lw <- function(x,theta,sigma,COV){
# Standard deviation of the conditional mean estimate: log-Weibull case
np <- length(theta); nc <- ncol(COV); nr <- nrow(COV)
if (np!=length(x)) cat("length(x) must be the same as length(theta)")
if (nc!=nr)        cat("COV is not a square matrix")
if (nc!=(np+1))    cat("ncol(COV) must be the same as length(theta)+1")
log.mu.x <- t(x)\%*\%theta+lgamma(1+sigma) # log of conditional mean estimate
mu.x     <- exp(log.mu.x)                # conditional mean estimate
dg       <- digamma(1+sigma)
COV.TT   <- COV[1:np,1:np]
Var.S    <- COV[(np+1),(np+1)]
COV.TS   <- COV[1:np,(np+1)]
V.mu.x   <- t(x)\%*\%COV.TT\%*\%x + dg^2*Var.S + 2*dg*(t(x)\%*\%COV.TS)
SD.mu.x  <- as.numeric((sqrt(V.mu.x))*mu.x)
SD.mu.x}

plt <- function(LOS,Cost,Adm,theta.fr,sigma.fr,sd0.fr,sd1.fr,theta.ml,
                sigma.ml,sd0.ml,sd1.ml){
# Plot of the conditional mean and confidence intervals: log-Weibull case
par(mfrow=c(1,2),oma=c(0,0,2,0))
plot(LOS,Cost,type="n")
points(LOS[Adm==0],Cost[Adm==0],pch=1)
points(LOS[Adm==1],Cost[Adm==1],pch=16,col=2)
x0t   <- x0\%*\%theta.fr; x1t <- x1t <- x1\%*\%theta.fr
lines(l0,exp(x0t)*gamma(1+sigma.fr))
lines(l0,exp(x1t)*gamma(1+sigma.fr),col=2)
z0min <- exp(x0t)*gamma(1+sigma.fr)-2.576*sd0.fr
z0max <- exp(x0t)*gamma(1+sigma.fr)+2.576*sd0.fr
z1min <- exp(x1t)*gamma(1+sigma.fr)-2.576*sd1.fr
z1max <- exp(x1t)*gamma(1+sigma.fr)+2.576*sd1.fr
lines(l0,z0min,lty=2,col=1)
lines(l0,z0max,lty=2,col=1)
lines(l0,z1min,lty=2,col=1)
lines(l0,z1max,lty=2,col=1)
polygon(c(l0,rev(l0)), c(z0min,rev(z0max)), border=FALSE, density=10, angle=90)
polygon(c(l0,rev(l0)), c(z1min,rev(z1max)), border=FALSE, density=12, angle=90,col=2)
plot(LOS,Cost,type="n")
points(LOS[Adm==0],Cost[Adm==0],pch=1)
points(LOS[Adm==1],Cost[Adm==1],pch=16,col=2)
x0t   <- x0\%*\%theta.ml; x1t <- x1t <- x1\%*\%theta.ml
lines(l0,exp(x0t)*gamma(1+sigma.ml))
lines(l0,exp(x1t)*gamma(1+sigma.ml),col=2)
z0min <- exp(x0t)*gamma(1+sigma.ml)-2.576*sd0.ml
z0max <- exp(x0t)*gamma(1+sigma.ml)+2.576*sd0.ml
z1min <- exp(x1t)*gamma(1+sigma.ml)-2.576*sd1.ml
z1max <- exp(x1t)*gamma(1+sigma.ml)+2.576*sd1.ml
lines(l0,z0min,lty=2,col=1)
lines(l0,z0max,lty=2,col=1)
lines(l0,z1min,lty=2,col=1)
lines(l0,z1max,lty=2,col=1)
polygon(c(l0,rev(l0)), c(z0min,rev(z0max)), border=FALSE, density=10, angle=90)
polygon(c(l0,rev(l0)), c(z1min,rev(z1max)), border=FALSE, density=12, angle=90,col=2)}

#------ End of auxiliary functions --------------------------------------------------

library(RobustAFT)

data(D243)
Cost <- D243$Cost                            # Cost (Swiss francs)
LOS <- D243$LOS                              # Length of stay (days)
Adm <- D243$Typadm; Adm <- (Adm==" Urg")*1   # Type of admission 
                                             # (0=on notification, 1=Emergency)
Ass <- D243$Typass; Ass <- (Ass=="P"   )*1   # Type of insurance (0=usual, 1=private)
Age <- D243$age                              # Age (years)
Dst <- D243$dest;   Dst <- (Dst=="DOMI")*1   # Destination (1=Home, 0=another hospital)
Sex <- D243$Sexe;   Sex <- (Sex=="M"   )*1   # Sex (1=Male, 0=Female)

\dontrun{
# Plot data
par(mfrow=c(1,2))
plot(LOS,Cost); plot(log(LOS),log(Cost))

# log-Weibull fits
# ----------------
# Full robust model
zwff     <- TML.noncensored(log(Cost)~log(LOS)+Adm+Ass+Age+Sex+Dst,
            errors="logWeibull")
summary(zwff)

# Reduced model
zwfr     <- update(zwff,log(Cost)~log(LOS)+Adm)
summary(zwfr)

# Residual plots
par(mfrow=c(1,2))
plot(zwfr,which=c(1,3))

# Plot robust predictions on log-log scale
par(mfrow=c(1,1))
l0       <- seq(from=2,to=60,by=0.5)
x0       <- as.matrix(cbind(1,log(l0),0))
x1       <- as.matrix(cbind(1,log(l0),1))
plot(log(LOS),log(Cost),type="n")
points(log(LOS[Adm==1]),log(Cost[Adm==1]),pch=16,col=2)
points(log(LOS[Adm==0]),log(Cost[Adm==0]),pch=1)
lines(log(l0),predict(zwfr,x0))
lines(log(l0),predict(zwfr,x1),col=2)

# Maximum likelihood : full model
zmlf     <- TML.noncensored(log(Cost)~log(LOS)+Adm+Ass+Age+Sex+Dst,
            errors="logWeibull",cu=100)
summary(zmlf)

# Maximum likelihood : reduced model
zmlr     <- update(zmlf,log(Cost)~log(LOS)+Adm)
summary(zmlr)

# Plot conditional means and confidence intervals
l0       <- seq(from=2,to=62,by=0.5)
x0       <- as.matrix(cbind(1,log(l0),0))
x1       <- as.matrix(cbind(1,log(l0),1))
theta.fr <- coef(zwfr)
sigma.fr <- zwfr$v1
COV.fr   <- vcov(zwfr)
sd0.fr   <- apply(x0,1,SDmux.lw,theta.fr,sigma.fr,COV.fr)
sd1.fr   <- apply(x1,1,SDmux.lw,theta.fr,sigma.fr,COV.fr)
theta.ml <- coef(zmlr)
sigma.ml <- zmlr$v1
COV.ml   <- zmlr$COV
sd0.ml   <- apply(x0,1,SDmux.lw,theta.ml,sigma.ml,COV.ml)
sd1.ml   <- apply(x1,1,SDmux.lw,theta.ml,sigma.ml,COV.ml)
plt(LOS,Cost,Adm,theta.fr,sigma.fr,sd0.fr,sd1.fr,theta.ml,sigma.ml,sd0.ml,sd1.ml)

# Gaussian fits (for comparison)
# ------------------------------
# Reduced model
zgfr <- TML.noncensored(log(Cost)~log(LOS)+Adm,errors="Gaussian")
summary(zgfr)

# Residual plots
par(mfrow=c(1,2))
plot(zgfr,which=c(1,3))

# Classical Gaussian fit
lr <- lm(log(Cost)~log(LOS)+Adm)
summary(lr)

# Compare several fits
# --------------------
comp <- fits.compare(TML.logWeibull=zwfr,ML.logWeibull=zmlr,least.squares=lr)
comp
plot(comp,leg.position=c("topleft","topleft","bottomleft")) # click on graphics
}
#
#------------------------------------------------------------------------------
#
# Example 2. This is the example described in Locatelli Marazzi and Yohai (2010).
# ---------
# This is the example described in Locatelli et al. (2010). 
# The estimates are slighty different due to changes in the algorithm for the 
# final estimate.
\dontrun{
# Remove data of Example 1
rm(Cost,LOS,Adm,Ass,Age,Dst,Sex)

data(MCI)
attach(MCI)
     
# Exploratory Analysis

par(mfrow=c(1,1)) 

plot(Age,log(LOS),type= "n",cex=0.7)

# (1) filled square   : regular,   complete   
# (2) empty  square   : regular,   censored
# (3) filled triangle : emergency, complete
# (4) empty  triangle : emergency, censored

points(Age[Dest==1 & TypAdm==0], log(LOS)[Dest==1 & TypAdm==0], pch=15,cex=0.7)  # (1)
points(Age[Dest==0 & TypAdm==0], log(LOS)[Dest==0 & TypAdm==0], pch=0, cex=0.7)  # (2)
points(Age[Dest==1 & TypAdm==1], log(LOS)[Dest==1 & TypAdm==1], pch=17,cex=0.7)  # (3)
points(Age[Dest==0 & TypAdm==1], log(LOS)[Dest==0 & TypAdm==1], pch=2, cex=0.7)  # (4)

# Maximum Likelihood
ML   <- survreg(Surv(log(LOS), Dest) ~ TypAdm*Age, dist="gaussian")
summary(ML)
B.ML <- ML$coef; S.ML <- ML$scale
abline(c(B.ML[1]        ,B.ML[3]        ),lwd=1,col="grey",lty=1)
abline(c(B.ML[1]+B.ML[2],B.ML[3]+B.ML[4]),lwd=1,col="grey",lty=1)
     
# Robust Accelerated Failure Time Regression with Gaussian errors
ctrol.S   <- list(N=150, q=5, sigma0=1, MAXIT=100, TOL=0.001,seed=123)

ctrol.ref <- list(maxit.sigma=2,tol.sigma=0.0001,maxit.Beta=2,tol.Beta=0.0001,
          Maxit.S=50, tol.S.sigma=0.001, tol.S.Beta=0.001,alg.sigma=1,nitmon=FALSE)

ctrol.tml <- list(maxit.sigma=50,tol.sigma=0.0001,maxit.Beta=50,tol.Beta=0.0001,
    Maxit.TML=50, tol.TML.sigma=0.001, tol.TML.Beta=0.001, alg.sigma=1,nitmon=FALSE)

WML       <- TML.censored(log(LOS)~TypAdm*Age,data=MCI,delta=Dest,otp="adaptive",
             control.S=ctrol.S,control.ref=ctrol.ref,control.tml=ctrol.tml)

summary(WML)

B.WML     <-coef(WML)
abline(c(B.WML[1]         ,B.WML[3]         ),lty=1, col="red")
abline(c(B.WML[1]+B.WML[2],B.WML[3]+B.WML[4]),lty=1, col="red")
#
detach(MCI)
}
}

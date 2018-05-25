"liepsu" <-
function(upper=.dFvGet()$upr) {
 epsi2 <- integrate(Psi2phi, -upper,upper)$value
 epsp  <- integrate(Pspphi,  -upper,upper)$value
list(epsi2=epsi2,epsp=epsp)
}

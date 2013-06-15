"E0.vect" <-
function(xbar,kl,ku,l,u) {
i1  <- integrate(RobustAFT:::Psiphi.w, lower=kl,upper=ku)$value
i2  <- integrate(RobustAFT:::Chiphi.w, lower=l,upper=u)$value
E0l <- matrix(c(i1*xbar,i2),ncol=1)
E0l}


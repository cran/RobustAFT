"E0.vect.n" <-
function(xbar,ku,u) {
i1  <- 0 #integrate(RobustAFT:::Psiphi.n, lower=-ku,upper=ku)$value
i2  <- integrate(RobustAFT:::Chiphi.n, lower=-u,upper=u)$value
E0l <- matrix(c(i1*xbar,i2),ncol=1)
E0l}


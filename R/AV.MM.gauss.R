"AV.MM.gauss" <- 
function(k0=1.5477, k1=4.6873, sigma=1)
{old <- comval()
 if (abs(k1 - 4.6873) < 0.0001) {
  Q1 <- 0.04509882; M1 <- 0.2069976 }
 else {
  dfcomn2(ipsi=4, xk=k1)
  Q1 <- integrate(Psi2phi, -10,10)$value
  M1 <- integrate(Pspphi,  -10,10)$value}
 if (abs(k0 - 1.5477) < 0.0001) {
  Q2 <- 0.1406683; M2 <- 0.3893525 }
 else {
 dfcomn2(ipsi=4, xk=k0)
 Beta <- integrate(Chiphi,  -10,10)$value
 Q2   <- integrate(Chi2phi, -10,10)$value-Beta^2
 M2   <- integrate(Chipzphi,-k0,k0)$value}
v.lambda <- (sigma^2 * Q1)/M1^2
v.sigma <- (sigma^2 * Q2)/M2^2
dfcomn2(ipsi=old$ipsi, xk=old$xk)
list(V.lambda = v.lambda, V.sigma = v.sigma)}


\name{dfcomn2}
\alias{dfcomn2}
\title{ Assigns values to the ROBETH parameters included in common blocks}
\description{
 See Marazzi A. (1993), p.405}
\usage{
dfcomn2(ipsi = -9, c = -1.345, h1 = -1.7, h2 = -3.4, h3 = -8.5, 
       xk = -1.548, d = -1.345, beta = -0.5, bet0 = -1, iucv = -1, 
       a2 = 0, b2 = -3, chk = -9, ckw = -2, bb = -1, bt = -1, 
       cw = -1, em = -1.345, cr = -2, vk = -1, np = -2, nu = -1, 
       v7 = -1, iwww = -1)
}
\arguments{
  \item{ipsi}{ Option parameter for the choice of  \eqn{\psi}{psi}. Set -4 <= ipsi <= 4 }
  \item{c}{ Parameter \code{c} of the Huber function }
  \item{h1}{ Parameter \eqn{h1} of the Hampel function }
  \item{h2}{ Parameter \eqn{h2} of the Hampel function }
  \item{h3}{ Parameter \eqn{h3} of the Hampel function }
  \item{xk}{ Parameter \eqn{k} of the rescaled Tukey biweight }
  \item{d}{ See reference}
  \item{beta}{ Parameter \eqn{\beta}{beta} to make\eqn{\sigma}{sigma} estimate 
               asymptotically unbiased }
  \item{bet0}{ Parameter \eqn{\beta_0}{bet0} to make\eqn{\sigma}{sigma} estimate 
               asymptotically unbiased }
  \item{iucv}{ Option parameter for the choice of u(s), u'(s), v(s), v'(s), w(s) or w'(s)}
  \item{a2}{ Parameter a^2 of Huber's mimimax u-function }
  \item{b2}{ Parameter b^2 of Huber's mimimax u-function }
  \item{chk}{ Parameter c of the Hampel-Krasker u-function }
  \item{ckw}{ Parameter c of the Krasker-Welsch u-function }
  \item{bb}{ Parameter b of the Mallows-unstandard u-function }
  \item{bt}{ Option parameter for w(s) or w'(s) }
  \item{cw}{ Option parameter for w(s) or w'(s) }
  \item{em}{ Parameter em for unstandard u-function }
  \item{cr}{ Parameter cr for unstandard u-function }
  \item{vk}{ Parameter vk for unstandard u-function }
  \item{np}{ Parameter np for unstandard u-function }
  \item{nu}{ Parameter nu for unstandard u-function }
  \item{v7}{ Parameter v for unstandard u-function }
  \item{iwww}{ Option parameter for the choice of \eqn{\bar\omega}{bar w}. 
               Set 0 <= iwww <= 3 }
}
\value{
  See reference
}

\references{
  Marazzi A. (1993) \emph{Algorithm, Routines, and S functions 
for Robust Statistics}. Wadsworth & Brooks/cole, Pacific Grove,
California. p.405
}
\keyword{utilities}


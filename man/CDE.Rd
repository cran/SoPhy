\name{ADE}
\alias{ADE}
\alias{CDE}
\title{Advection-dispersion equation}
\description{
  Solution of a simple advection-dispersion equation
}
\usage{
  ADE(z, time, C0, dispersion, velocity) 
}

\arguments{
  \item{z}{vector of distances from surface}
  \item{time}{vector of time points}
  \item{C0}{constant concentration at the surface}
  \item{dispersion}{constant dispersion coefficient}
  \item{velocity}{constant velocity}
}
\details{
  \command{ADE} solves the PDE 
  \deqn{\frac{\partial C(z,t)}{\partial t} = D \frac{\partial^2
      C(z,t)}{\partial z^2}- v\, \frac{\partial C(z,t)}{\partial z}}{d
    C(z,t) / d t =  D * d^2 C(z,t) / d z^2  - v * d C(z,t) / d z}
  on the one-dimensional semi-finite column, i.e.
  \eqn{z\in[0,\infty]}{z in [0,infty]}.
  Here \eqn{D} is the dispersion coefficient and \eqn{v} the
  velocity. Further, \eqn{C(\infty, t)=0}{C(infty,t)=0} for all \eqn{t}
  and \eqn{C(z,0)=0} for all \eqn{z}.
  Further, \eqn{C(0, t)=C_0}{C(0,t)=\code{C0}}.
}
\value{
  \command{ADE} calculates the concentration on a grid given by
  \code{z} and \code{time}, i.e.,
  \command{ADE} returns a (\code{length(z)} x \code{length(time)})-matrix. 
}

\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}
\references{
  Schlather, M. and Huwe, B. (2004)
  The use of the language interface of R: two examples for modelling water
  flux and solute transport. \emph{Computers \& Geosciences} \bold{30}, 197-201.
  
  Tindall, J.A. and Kunkel, J.R (1999) \emph{Unsaturated Zone Hydrology
    for Scientists and Engineers.} New Jersey: Prentice-Hall. % p. 280
 }

\examples{
par(cex=1.5, mar=c(4.2,4,0.2,0.4))
z <- ADE(z=1:100, t=10^(0:3), C0=1, d=1, v=1)
matplot(z, outer(-1:-100, rep(1,4)), type='l',
        xlab='concentration', ylab='z', 
        lwd=3, cex.lab=1.5, cex.axis=1.5, col=1)
legend(1.0, -100, legend=paste('t=', 10^(0:3), sep=''), xj=1, yj=0,
        lwd=3, lty=1:4, cex=1.5, col = 1)
}
\keyword{spatial}

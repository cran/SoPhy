\name{my.legend}
\alias{my.legend}
\title{legend for a spectrum of colours}
\description{
  the function plots a legend for a spectrum of colours
}
\usage{
my.legend(lb.x=0, lb.y=0, zlim, col, cex=0.8, y.intersp=0.02, bg='white') 
}
\arguments{
  \item{lb.x}{x-coordinate of the left bottom corner of the legend box}
  \item{lb.y}{y-coordinate of the left bottom corner of the legend box}
  \item{zlim}{value range for the colours}
  \item{col}{spectrum of colours, obtained by
    \code{\link{rainbow}} or \code{\link{grey}}, for example.}
  \item{cex}{letter size increse}
  \item{y.intersp}{the space between the small colour lines; this value
    might be changed if the size of the legend box is not appropriate.
  }
  \item{bg}{background colour}
}
\value{
  \code{NULL}
}
%%% unerklaerlicher fehler bei R CMD check, falls nicht auskommentiert.
%%% spaeter wieder probieren.
%\examples{
% library(SoPhy, lib="~/TMP")
%  x <- seq(0, 100, 1)
%  y <- seq(0, 100, 1)
%  xy <- expand.grid(x, y)
%  z <- matrix(nrow=length(x), sin(xy[,1]/3) + cos(xy[,2]/2 + xy[,1]/4))
%  zlim <- range(z)
%  col <- rainbow(100)
%  image(x, y, z, zlim=zlim, col=col)
%  my.legend(zlim=zlim, col=col, y.i=0.03)
% }
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}
\examples{
 ## see \link{read.swms2d.table}, for example
}
\keyword{hplot}

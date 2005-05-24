\name{calculate.horizons}
\alias{calculate.horizons}
\title{Update of the horizon structure}
\description{
  The function updates the positions of the boundary points, the
  affiliation of the pixels to the horizons and
  the location and the size of the smallest rectangle
  that include a horizon.
}
\usage{
calculate.horizons(h)
}
\arguments{
  \item{h}{a list of the same format as the output of 
    \command{\link{xswms2d}}}
}
\value{
  a list of the same format as the output of \command{\link{xswms2d}}.
  The values of the following list elements are changed:
  \item{h[[i]]$points}{points defining the vertices of
    the boundary line for horizon \eqn{i};
    changed if x-coordinates of a boundary definition are
    all genuinely within the grid size} 
 
  \item{idx.rf}{affiliation number of each grid point} 
 
  \item{h[[i]]$border}{grid points defining the boundary of horizon \eqn{i}} 
   
  \item{h[[i]]$cut.x, h[[i]]$cut.y}{clipping rectangle that includes the
    horizon \eqn{i}} 
 
  \item{h[[i]]$idx}{indication of the points within the clipping area that
     belong to horizon \eqn{i}} 
    
  \item{RF}{\code{NULL}} 
 
  \item{Stone.R}{\code{NULL}} 
 
  \item{Root.RF}{\code{NULL}} 
}
\note{
  The function must be called after the definition of the horizons has
  been changed by the user.
  If a first or a last defining point for a boundary is truely within
  the range of the x-coordinates of the grid then a point with the
  respective x coordinate on the boundary of the grid is added. 
}
\seealso{
  \command{\link{draw.horizons}},
  \command{\link{modify.horizons}},
  \command{\link{plotRF}},
  \command{\link{plotWater}},
  \command{\link{simulate.horizons}},
  \code{\link{SoPhy}},
  \command{\link{swms2d}} 
  \command{\link{xswms2d}} 
}
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}
\examples{
%  library(SoPhy, lib="~/TMP"); library(RandomFields)
% source("/home/schlather/R/SOPHY/SoPhy/R/Sophy.R")
% source("/home/schlather/R/SOPHY/SoPhy/R/swms2d.R")
% source("/home/schlather/R/SOPHY/SoPhy/R/simu.R")

############ fetch the standard definition of a profile ###########
h <- xswms2d(xlim=c(1, 100), ylim=c(1, 100), step=1, new=NULL)
##        new=NULL: xswms2d returns the standard definition for h
##                       without entering the interactive surface

############ modify the profile definition ########################
h$n <- 2         ## define a second horizon
h$H2 <- h$H1    
h$H2$type <- "H" ## genuine horizon, not polygon
h$H2$model$model[[1]]$scale <- 30 ## different
## coordinates of the boundary segments between the horizons: 
h$H2$points <- list(x=seq(1, 100, 5))             ## x coordinates
h$H2$points$y <- 80 + 20 * cos(h$H2$points$x / 10) ## y coordinates
h <- calculate.horizons(h)       ## update the internal parameters
                         ## after modifying the profile definition
draw.horizons(h)                 ## plot the profile

############ simulate water flux for the specified profile ########
h <- simulate.horizons(h)        ## simulate stochastic components
plotRF(h)                        ##   plot the simulated random field
swms2d.out <- swms2d(h, iter.print=1) ## numerical simulation
if (is.character(swms2d.out)) {
  cat(swms2d.out, "\n")          ## an error has occured 
} else {                         
  plotWater(swms2d.out, what="theta") ## plot the water contents 
}

}
\keyword{spatial}

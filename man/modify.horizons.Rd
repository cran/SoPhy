\name{modify.horizons}
\alias{modify.horizons}
\title{
  Stochastic variation of the definitions of the horizons
}
\description{
  The function modies the geometric definition of the horizons stochastically
}
\usage{
modify.horizons(h, percent=5, level.percent=5, rdistr=rnorm)
}

\arguments{
  \item{h}{list of the same format as the output of
    \command{\link{xswms2d}}}
  \item{percent}{magnitude of variation of a single point; see Details.}
  \item{level.percent}{magnitude of the length of the translation vector
    for the whole object; see Details.}
  \item{rdistr}{distribution for the translation.  The function must
    have two parameters:  the first parameter gives the number of random
    variables to generate;  the other parameter is called by name, \sQuote{s},
    and should determine (or be related to) the standard deviation of
    the distribution.}
}
\details{
  The function \command{modify.horizons} changes stochastically the
  location of the vertices of the defining segments for the
  boundary of a horizon. To this end, 
  the standard deviations \eqn{s_x}{s.x} and \eqn{s_y}{s.y}
  for the x coordinates and the y coordinates of the vertices are calculated.
  Then the \code{s} parameter for the distribution of the translation distance
  equals \eqn{s_x \cdot \code{percent} / 100}{s.x * \code{percent} /
    100}
  for the x coordinates (analogously for the y coordinates).
  \itemize{
    \item
    For genuine horizons there is also a shift of all the points in
    vertical direction by the distibution \code{rdistr(,s)} where
    is the product of
    the total length of the grid in y direction and \code{level.percent}
    devided by 100 times the number of genuine horizons.
    (Without the devision by the number of genuine horizons the user would
    have to pay much more attention that the horizons do not awefully overlap.)
    Further, if the first point is right to the second or the last point
    is left to the last but one point, the x coordinate of the
    first and/or the second point is replaced by the value of the
    neighbouring point minus or plus a random value that is uniformly 
    distributed on \eqn{[0, s_x]}{[0, s.x]}. Finally, if the x coordinate of
    the first or the last point is within the grid than points on 
    the border of the grid are added.

    \item
    If it is a polygon, then the whole figure is additionally
    shifted to the right and top according to the given distribution where
    \code{s} equals (the range of the
    respective coordinate) * \code{level.percent} / (100 * (the number
    of polygons)). 
}
  
  \command{modify.horizons} calls \command{\link{calculate.horizons}}
  before returning the list.
}
\value{
  List of the same format as the output of \command{\link{xswms2d}}.
  
  \command{modify.horizons} calls \command{\link{calculate.horizons}}
  before returning the list.
}

\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/~schlather}}

\seealso{\code{\link{SoPhy}}, \command{\link{xswms2d}}}
\examples{
%  library(SoPhy, lib="~/TMP"); library(RandomFields)
% source("/home/schlather/R/SOPHY/SoPhy/R/Sophy.R")
h <- xswms2d(xlim=c(1, 100), ylim=c(1, 100), step=1, new=NULL)
##        new=NULL: xswms2d returns the standard definition for h
##                          without entering the interactive surface
h$n <- 2         ## define a second horizon
h$H2 <- h$H1    
h$H2$type <- "H" ## genuine horizon, not polygon

## coordinates of the boundary segments between the horizons 
h$H2$points <- list(x=seq(1, 100, 5))             ## x coordinates
h$H2$points$y <- 40 + 20 * cos(h$H2$points$x / 5) ## y coordinates

h <- calculate.horizons(h) ## update the internal parameters
draw.horizons(h)           ## plot the horizons
for (i in 1:20) {
  readline("Press return")
  h <- modify.horizons(h)  ## stochastic modification of the horizons
  draw.horizons(h)        
  }
}
\keyword{spatial}


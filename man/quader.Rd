\name{Quader}
\alias{belichtet}
\alias{pos3D}
\alias{quader}
\title{Perspective plot}
\description{Simple functions that plot points in a perspective view
}
\usage{
   belichtet(x, radius=1, fix.scale=TRUE, correction=1.2)

   pos3D(x, inf)

   quader(size, inf, sun, bottomleftfront= c(0,0,0), 
          col=grey(seq(1,0,-0.001)), col.frame=c("grey", "black"), lty=c(2,1),
          cex.axis=1.5,  reverse=TRUE, unit = "cm", add=FALSE, plot=TRUE,
                   font = c("sans serif", "bold"), srt=NULL)

}
\arguments{
  \item{x}{matrix of three rows; the three dimensional coordinates}
  \item{radius}{radius of the dyed tube around each path;
     the tube is approximated by a sequence of balls}
   \item{fix.scale}{logical. If \code{FALSE}
     the plotted radius is modified according to the perspective distance
   }
   \item{correction}{positive real number. Do not change its value.
     This constant takes into account that the plotting procedure used
     in the function does not plot a circle of exactly the given radius.
   }
  \item{inf}{Vector of three components. The first two give the position
    of the infinitely far point; the third component gives the extension
    factor for the second dimension that is orthogonal to the screen;
    \code{inf} is of order \code{1/100}}
  \item{size}{dimension of the quader}
  \item{sun}{the three dimensional coordinates of the fictive sun}
  \item{bottomleftfront}{location of the quader; the standard is that
    the quader is in the first quadrant with one vertice in the origin}
  \item{col}{vector of colours; the point that is nearest to the sun
    gets the colour \code{col[1]}, the one that is the furthest away
    \code{rev(col)[1]}}
  \item{col.frame}{vector of two components. 
    The edges that are visible get colour
    \code{col.frame[2]} the others \code{col.frame[1]}} 
  \item{lty}{vector of two components.  the edges that are visible have
    line type \code{lty[2]} the  others \code{lty[1]}}
  \item{cex.axis}{The magnification to be used for axis annotation
    relative to the current, see \code{\link{par}}}
  \item{reverse}{logical. If \code{TRUE} the z-axis is numbered
    upside down}
  \item{unit}{}
  \item{add}{logical.  If \code{TRUE}, a quader is added to current figure}
  \item{plot}{logical. If \code{FALSE} only the values for the
    coordinate system are calculated and stored in the background,
    which will be used by \command{belichtet}, for example.
  }
  \item{font}{font used for the labelling of the axes, see also
    \code{\link[graphics]{font}}}
  \item{srt}{angle for the character rotation when labelling the y axis;
    if \code{NULL} the angle is calculated, but might be not always the
    best choice.
  }
}
\value{
  The functions \code{belichtet} and \code{quader} return \code{NULL},
  the function \code{pos3D} returns a matrix of 2 rows and
  \code{ncol(x)} coloumns.
}
\details{
  The function \code{quader} must be called first with \code{add=FALSE}.
  \code{quader} plots a perspective quader and provides global variables
  whose name starts with \code{.p3d}.
  If, for the second time, \code{quader} is called solely by the argument
  \code{add=TRUE}, the original quader is redrawn.

  The auxiliary function \code{pos3D} calculates the two-dimensional
  projected positions of the three dimensional points. If \code{inf} is
  not given, the one given in \code{quader} is used.

  \code{belichtet} plots the three dimensional points in perspective view.
}
\examples{
  % library(SoPhy, lib="~/TMP")
  f <- function(x) 
     cos((x[2,] + x[1, ]) / 5) + cos((x[2, ]-x[3,]) / pi + 3) +
     cos(x[3, ] / 2 + 2)
  x <- t(expand.grid(1:50, 1:100, 1:50))
  y <- x[, f(x)>0.5]
  quader(size=c(50, 100, 50), inf=c(150, -400, 1/500),
         sun=c(200, -30000, 50))
  belichtet(y)
  quader(add=TRUE, col.frame=c("transparent", "black"))
}
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}
\keyword{hplot}


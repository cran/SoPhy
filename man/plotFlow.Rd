\name{plotFlow}
\alias{plotFlow2d}
\alias{plotFlow3d}
\title{plotting of the simulated infiltration pattern}
\description{
  plotFlow2d plots two dimensional profiles of simulated stochastic flow
  patterns;
  plotFlow3d plots a 3d, perspective figure of simulated stochastic flow
  patterns.
}
\usage{
plotFlow2d(coord, pointradius=1, slice=2 * pointradius,
           full.size=TRUE, Profiles=1, dev=1, ps="",
           height=4, unit="cm", cex=2, correction=1.2,
           rl = function(x) readline(paste(x, ": press return")))

plotFlow3d(paths, horizons=c("no", "absorbing", "breakthrough"),
           drop.distr, n.balls=1, pointradius=1,
           dev=1, ps="3d.dye.pattern", ps.background=FALSE,
           profileheight=4, unit="cm", unit.scale=1, inf, sun,
           rl = function(x) readline(paste(x, ": press return")),
           low.resolution = TRUE)
}
\arguments{
  \item{coord}{output of \command{plotFlow3d}: \eqn{n x 3}{n x 3} matrix;
     3d coordinates for the centre of the balls that approximate the tube around
     each path}
   \item{pointradius}{radius of the dyed tube around each path;
     the tube is approximated by a sequence of balls}
   \item{slice}{thickness of the slice for which the paths are visible}
   \item{full.size}{logical. If \code{FALSE} then the projected radii
     of the sliced balls are plotted. If \code{TRUE}
     then instead of the projected radius the \code{pointradius}
     is used.
   }
   \item{Profiles}{number of the in the
     y-direction equally spaced profiles of the
     3d flow pattern}
   \item{dev}{plotting device, see \code{\link[RandomFields]{Dev}}
   }
   \item{ps}{name of the postscript file, see \code{\link[RandomFields]{Dev}}
   }
   \item{height}{height of the plotted figure, see
     \code{\link[RandomFields]{Dev}}
   }
   \item{unit}{string; if empty then labels are not given. Otherwise
     \code{unit} denotes the units of the plot}
   \item{cex}{used for \code{cex.axis} and \code{cex.lab}, see
     \code{\link[graphics]{par}}
   }
   \item{correction}{positive real number. Do not change its value.
     This constant takes into account that the plotting procedure used
     in the function does not plot a circle of exactly the given radius.
   }
   \item{paths}{list as returned by \command{\link{flowpattern}}}
   \item{horizons}{
     Temptative parameter.
     If \code{"absorbing"} the path will be \sQuote{infinitely} long
     in a layer at depth 50 about, so that fluxes will stop in this layer
     If \code{"breakthrough"} the path will be very long
     in a layer at depth 50 about, so that only a few paths are
     percolated.
     
     The values \code{"absorbing"} and
     \code{"breakthrough"} are used to create the figures in the
     Discussion of Schlather and Huwe (2004a).
   }
   \item{drop.distr}{function or missing.
     The stained path length reached by a single drop is
     calculated as follows.  First a Gaussian random field is simulated
     with the parameters given below. Then the Gaussian random field is
     marginally transformed to uniform distribution by the
     Gaussian distribution function. Then the inverse of the required
     distribution function is applied.  \code{drop.distr} gives this
     inverse.  It should allow for vectors as input variables.
     The \code{drop.distr} used to create \code{paths} is
     overwritten if \code{drop.distr} is not missing, but
     the realisation of the paths is always kept.
   }
   \item{n.balls}{the tube for each path segment is approximated by
     \code{n.balls} balls}
   \item{ps.background}{logical. If \code{dev=TRUE}
     then a postscript file is created, which is, in general, huge.
     Therefore, any postscript file created on a unix system
     is transformed into a tiff-file by `convert'
     and back again to a postscript file. These transformation are
     time consuming. If \code{ps.background=TRUE} then
     these transformations are started in the background.
   }
   \item{profileheight}{essentially the height of the figure; the value is
     used to calculated the width of the plot. The final height of the
     plot is higher then \code{profileheight} because additional space
     needed for the perspective plot}
   \item{unit.scale}{scaling factor for the numerical values of
     all three-dimensional coordinates}
   \item{inf}{Vector of three components. The first two give the position
    of the infinitely far point; the third component gives the extension
    factor for the second dimension that is orthogonal to the screen;
    \code{inf} is of order \code{1/100}}
  \item{sun}{the three dimensional coordinates of the fictive sun}
  \item{rl}{function. Called, after a figure, but the last,
    is plot on the screen.}
  \item{low.resolution}{logical. If \code{TRUE} and \code{dev=FALSE},
    i.e. postscript file will be produced, then on unix systems the
    postscript file will be rewritten to have much smaller size, but
    also lower resolution.
  }
}
\value{
  \command{plotFlow2d} returns \code{NULL}.
  \command{plotFlow3d} returns invisibly a matrix for the centers of the
  approximating balls
}
%\details{
%}
\seealso{\command{\link{flowpattern}}}

\references{
  \itemize{
    \item
    Schlather, M. and Huwe, B. (2004a)
    A stochastic model for 3-dimensional flow patterns in dye tracer
    experiments. \emph{Submitted to J. Hydrol. }% \bold{}, .
    
    \item
    Schlather, M. and Huwe, B. (2004b)
    A risk index for characterising flow pattern in soils
    using dye tracer distributions
    \emph{Submitted to J. Contam. Hydrol.} % \bold{}, .
   
  }
}
\examples{
 ## see \link{flowpattern}
}
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}
\keyword{hplot}

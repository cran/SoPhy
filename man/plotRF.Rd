
\name{plotRF}
\alias{plotRF}
\title{Plotting random fields}
\description{
  the function plots random fields including stones and roots using
  the definition of the horizons in \command{\link{xswms2d}}.
}
\usage{
  plotRF(h, col.txt='black', col.stones='white',
         col.rf=if (is.null(h$col.rf)) rainbow(100) else h$col.rf, 
         titl=TRUE, line = 0.7, lim=1, quadratic=TRUE,
         cex=1, cex.leg=0.8, pch=5,
         cex.pch = max(0.05, 3/min(length(h$grid.x), length(h$grid.y))),
         root.col=grey(1 - (i-1) / (1.6 * length(h$plants))),
         what=c('Root.RF', 'Stone.RF', 'RF'),
         transf=c("K", "H", "theta", "none"), legend = TRUE, ylim, zlim)
 }
\arguments{
  \item{h}{a list of the
    same format as the output of \command{\link{xswms2d}}}
  \item{col.txt}{colour of the title}
  \item{col.stones}{colour of the stones and the air}
  \item{col.rf}{colour spectrum of the random field;  plants are always grey}
  \item{titl}{logical. If \code{TRUE}  a title is plotted using the
    colour \code{col.txt}}
  \item{line}{place where the titel is plottes, see
    \command{\link[graphics]{title}}.}
  \item{lim}{value in \eqn{[0,1]} that is used only if
    \code{zlim} is missing; quantile for the upper bound 
    \code{zlim[2]}; the value of \code{lim} should be
    less than 1 if the random field contains some extreme values, since
    the colour scale is linear}
  \item{quadratic}{logical.  If \code{TRUE} the figure matrix is enlarged and
    filled with NA symmetrically in x-direction or on the bottom in
    y-direction such the matrix (and the figure) become quadratic}
  \item{cex}{A numerical value giving the amount by which the title
    text, the axes and the
    labels of the axes should be scaled relative to the default}
  \item{cex.leg}{A numerical value giving the amount by which the legend text
    should be scaled relative to the default}
  \item{pch}{the symbol for a root segment}
  \item{cex.pch}{size of the symbol}
  \item{root.col}{the colours for a root segments}
  \item{what}{chooses between the random field that includes modifications by
    roots and stones or by stones allone, or the random field without
    modifications}
  \item{transf}{Kind of random field:
    saturated hydraulic conductivity \code{"K"},
    pressure head \code{"H"},
    the water content \code{"theta"},
    or the original Gaussian random field (\code{"none"})}
  \item{legend}{logical. If (TRUE) a legend is added.}
  \item{ylim}{missing or vector of two components. If missing the range
    of the simulation is used.}
  \item{zlim}{missing or vector of two components. Values below
    \code{zlim[1]} are plotted in white,
    values above \code{zlim[2]} are plotted also in white.
    If missing then \code{zlim[1]} is the minimal value
    and \code{zlim[2]} is calculated by means
    of \code{lim}.
  }
}
\details{
  If \code{transf="K"} then
  the function plots the random field of hydraulic conductivity, i.e.,
  the Gaussian random field after
  \code{h$miller.link} and \code{h$millerK} have been applied.
  }
\value{
  the matrix of the
  random field values, always excluding the plants,
  if no error has occured
  and an error message otherwise.
}
\seealso{
  \command{\link{plotWater}}, \command{\link{xswms2d}}
}
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}
\examples{
   ## see \link{simulateHorizons}
  }
\keyword{hplot}


\name{draw.horizons}
\alias{draw.horizons}
\title{Draw of the profile of horizons}
\description{
  The function draws a profile, i.e. a sequence of horizons,
  by coloured areas or by boundary lines.
}
\usage{
 draw.horizons(h, areas=TRUE,
               col.hor=c('#000000', '#996600', '#660000', '#CC9933',
                         '#666600', '#CCCC99', '#CCCCCC', '#990000',
                         '#FFCC00', '#FFFFCC', rep('white', h$max.horizons)),
               border.col=NULL, picture=NULL,
               lwd = 2, quadratic=TRUE, all=TRUE, cex = 1, cex.leg = 1
               )
}

\arguments{
  \item{h}{a list of the same format as the output
    of \command{\link{xswms2d}}}
  \item{areas}{logical.  If \code{TRUE} the horizons are plotted in
    different colours.  If \code{FALSE} only the boundary lines are
    plotted.}
  \item{col.hor}{the colours of the horizons or the boundary lines;
    the vector must have 20 entries.  The first ten give the
    colours of the horizons, the second ten give the colour of the stones
    that belong to a certain horizon (here the position of the pixel not
    the center of the stone counts).
    }
  \item{border.col}{\code{NULL} or any vector of colours.  If not
    \code{NULL} the boundary pixels are plotted in the colour of the
    respective horizon.}
  \item{picture}{array of three dimensions where the third
    dimension has 3 or 4 components (RGB or RGBA coding).
    The pictures is given as a background figure if \code{areas=FALSE} and
    \code{all=TRUE}.
  }
  \item{lwd}{postive number, only used if \code{area=FALSE}.  Then it
    gives the width of the boundary lines.}
  \item{quadratic}{logical.  If \code{TRUE} the figure matrix is enlarged and
    filled with NA symmetrically in x-direction or on the bottom in
    y-direction such the matrix (and the figure) become quadratic.}
  \item{all}{If \code{FALSE} no axis and no background picture are
    plotted.}
  \item{cex}{type size for the axes and the labels, see
    \command{\link[graphics]{par}}}
  \item{cex.leg}{type size for the legend}
}
\value{
  invisible list of \code{xlim} and \code{ylim} of the figure.
}

\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/institute}}

\seealso{
  \command{\link{modify.horizons}},
  \code{\link{SoPhy}},
  \command{\link{xswms2d}}
}
\examples{
   ## see  \link{modify.horizons}  for an example
}
\keyword{spatial}












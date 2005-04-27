\name{plotWater}
\alias{plotWater}
\title{Plotting the result of SWMS2D}
\description{
  the function plots the H, Q, \eqn{\theta}{theta}, vx, vz, or conc,
  stored in a list as in the output of SWMS2D
}
\usage{
plotWater(h, instance,
          what=c("H", "Q", "theta", "vx", "vz", "Conc", "logH"),
          col.txt="black",
          col.simu=if (is.null(h$col.simu)) rainbow(100) else
                   if (what.nr == 7) rev(h$col.simu) else h$col.simu,
          col.exception = c("yellow", "darkred"), lim=1,
          titl=TRUE, line=0.2, quadratic=TRUE, cex=1, cex.leg=0.8,
          legend = TRUE, ylim, zlim)
}
\arguments{
  \item{h}{a list as returned by \code{\link{swms2d}} or a list as
    returned by  \code{\link{xswms2d}};
    usually \code{\link{create.waterflow}} has been called beforehand.
  }
  \item{instance}{
    the number of the time point for which the data should be plotted.
    If missing or \code{NULL} the last one is taken.
  }
  \item{what}{character.
    \code{what} gives the parameter to be plotted:
    \describe{
      \item{\code{H}}{water potential \eqn{H},}
      \item{\code{Q}}{discharge/recharge rates \eqn{Q} for internal
	sink/sources,} 
      \item{\code{theta}}{water contents \eqn{theta},}
      \item{\code{vx}}{x-components of the Darcian flux vector vx,}
      \item{\code{vz}}{z-components of the Darcian flux vector vz,}
      \item{\code{Conc}}{solute concentration}
    }
    The parameter \code{what} allows also the numbers \eqn{1,\ldots,6}{1,...,6}
    corresponding to \code{H},...,\code{Conc}.
  }
  \item{col.txt}{colour of the title}
  \item{col.simu}{colour spectrum of the simulated field of
    H, Q, \eqn{\theta}{theta}, vx, vz, or conc; areas of stones or air
    are white.}
  \item{col.exception}{vector of 2 components; colour plotted values
    that are below or above the given range, see \code{zlim}}
  \item{lim}{value in \eqn{[0,1]} that is used only if \code{zlim} is
    missing or \code{NULL};
    the range \eqn{\code{zlim}} of the plotted values depends on \code{what}.
    Let \eqn{q(p)} be the \eqn{p}-quantile of the values. Then 
    \itemize{
      \item \code{H} : \eqn{\code{zlim}=[q(1-\code{lim}); 0]}
      \item \code{theta} :
      \eqn{\code{zlim}=\command{range}(\{0,
	\theta_{s,i}\})}{\code{zlim}=\command{range}({0,
	  theta_{s,i}})}, where the \eqn{i} runs over all horizons
      and polygons. That is, the range is independent of the value of
      \code{lim}. See \code{\link{xswms2d}} for \eqn{\theta_s}{theta_s}.
      \item \code{Conc} : \eqn{\code{zlim}=[0, q(\code{lim})]}
      \item otherwise :
      \eqn{\code{zlim}=[q(1-\code{lim}), q(\code{lim})]}
    }
    quantile for the upper bound
    of the values of the random field; this variable should be
    less than 1 if the random field contains some extreme values, since
    the colour scale is linear}
  \item{titl}{logical or character. If \code{TRUE} a title is plotted using the
    colour \code{col.txt}}
  \item{line}{parameter of function \command{\link[graphics]{title}}}
  \item{quadratic}{logical.  If \code{TRUE} the figure matrix is enlarged and
    filled with NA symmetrically in x-direction or on the bottom in
    y-direction such the matrix (and the figure) become quadratic}
  \item{cex}{A numerical value giving the amount by which the title
    text, the points indicating the root segments, the axes and
    labels of the axes should be scaled relative to the default}
  \item{cex.leg}{A numerical value giving the amount by which the legend text
    should be scaled relative to the default}
  \item{legend}{logical. If (TRUE) a legend is added.}
  \item{ylim}{missing or vector of two components. If missing the range
    of the simulation is used.}
  \item{zlim}{missing or vector of two components. Values below
    \code{zlim[1]} are plotted in \code{col.exception[1]},
    values above \code{zlim[2]} are plotted in \code{col.exception[2]}.
    If missing or \code{NULL} the vector \code{zlim} is calculated by means
    of \code{lim}.
  }
}
\details{	
  For all variables the plot is over the whole range except for
  the following cases.
  If \eqn{theta} is plotted,
  the minimum is 0 and the maximum the maximal field
   capacity, If H is plotted, the maximum is at most 0.
}
\value{
  matrix of the plotted values
  if no error has occured and an error message otherwise.
}
\seealso{
  \code{\link{plotRF}}, \code{\link{xswms2d}}
}
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}
\examples{
   ## see  \link{create.waterflow}
  }
\keyword{hplot}

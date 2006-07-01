
\name{create.stones}
\alias{create.stones}
\title{Simulation of stones}
\description{
  This function simulates stones within a profile 
  according to a given soil definition.
  The algorithm returns an intermediate stage of the random sequential
  adsoption (RSA) model.
}
\usage{
create.stones(h, trials = 10)
}
\arguments{
  \item{h}{a list of the
    same format as the output of \command{\link{xswms2d}}}
  \item{trials}{parameter for the
    number of trials after which the simulation for a horizon terminates
    and a warning is given.
    The number of trials equals, for each horizon,
    \code{trials} times the Poisson random number of
    stones to simulate, see also Details.
    }
}
\details{
  In the simulation,
  each horizon is considered separately in the ordering they are
  defined.
  The ordering is important in case the stones are allowed to reach into
  neighbouring horizons.
  For each horizon:

  \enumerate{
    \item a Poisson number \eqn{p} of stones is drawn with mean
    \code{h$stone$lambda} times modulus of the area of the horizon [in
    the given units].
    
    \item while not all \eqn{p} points are simulated and the total number of
    trials for the horizon is less than \code{trials} times \eqn{p}:

    \enumerate{
      \item within the rectangle that is given by \code{cut.x} and
      \code{cut.y} and that encloses the area of the horizon a uniformly
      distributed point \eqn{P} is draw.  If the point is outside the horizon,
      the drawing is repeated, up to \code{trials} times.  If still not
      successful, the simulaton for a horizon terminates, a warning is
      given, the stones
      simulated up to now are kept and the next
      horizon is simulated.
      
      \item According to
      \cr
      \code{$stone$phi.distr}, \code{$stone$phi.mean},
      \code{$stone$phi.s},
      \cr
      \code{$stone$main.distr}, \code{$stone$main.mean},
      \code{$stone$main.s},
      \cr
      \code{$stone$sec.distr}, \code{$stone$sec.mean}, and
      \code{$stone$sec.s},
      \cr
      an ellipse is simulated with centre \eqn{P}.
      The length of the main axis is given by a random variable with
      distribution \code{stone$main.distr}, mean parameter
      \code{$stone$main.mean} and standard deviation
      \code{$stone$main.s}.
      The length of the secondary axis and
      the angle between the main axis and the horizon are given in
      a similar way. 
      
      \item The pixels that are covered by the ellipse are determined.

      \item The ellipse including \eqn{P} is rejected in the following
      cases

      \enumerate{
	\item x
	\code{$stone$no.overlap} is \code{TRUE} and any of the
	pixels belong already to a stone
	
	\item \code{$stone$no.lower} is \code{TRUE} and any pixel
	belongs to a horizon or polygon that precedes the current one.
	
	\item \code{$stone$no.upper} is \code{TRUE} and any pixel
	belongs to a horizon or polygon that follows the current one.
      }
      
    }
    
    \item if the total number of trials is reached without having
    \eqn{p} points, the algorithm fails.
  }
}
\value{
  The function returns the input argument \code{h}
  except that \code{h$idx.rf} and \code{h$Stone.RF} have been updated.
  
  \item{\code{$idx.rf}}{
    gives, at first instance, the association of a pixel
    to a horizon (0,...,\code{h$n}-1).
    If the pixel belongs to one or more stones the value
    of the maximum number of horizons is added.}

  \item{\code{$Stone.RF}}{
    \code{RF}, modified at the locations with stones
    according to the \code{h[[...]]$stone$value}s.
  }
}
\references{
  
  Schlather, M. and Huwe, B. (2006) Modelling inhomogeneous soils:
  Theory. \emph{In preparation.} %\bold{}, - .
  \cr

  Stoyan, D., Kendall, W.S., Mecke, J. (1995)
  \emph{Stochastic Geometry and its Applications} Chichester: Wiley, 2nd
  edition.
  \cr
  
  RSA model (also called SSI process in stochastic geometry)
  \itemize{
    \item
    Herczynski, R. (1975) Distribution function for random distribution
    of spheres. \emph{Nature}, \bold{255}, 540-541.
  }
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/institute}}
\seealso{\code{\link{SoPhy}}, \command{\link{xswms2d}}}

\examples{
%  library(SoPhy, lib="~/TMP"); library(RandomFields)
% source("/home/schlather/R/SOPHY/SoPhy/R/Sophy.R")
## get a standard definition first
h <- xswms2d(xlim=c(1, 100), ylim=c(1, 100), step=1, new=NULL)
h <- simulateHorizons(h)
h$H1$stone$lambda <- 0.002       ## in average 20 stones 
h$H1$stone$no.overlap  <- FALSE ## allow for overlapping stones
for (i in 1:5) {
  plotRF(create.stones(h), what="Stone") 
  readline("Press return")
}
}
\keyword{spatial}

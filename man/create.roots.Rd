\name{create.roots}
\alias{create.roots}

\title{Simulation of roots}
\description{
  This function generates root systems within a profile 
  according to a given specification of the soil and the root model.
}
\usage{
create.roots(h, trials = 10, PrintLevel=RFparameters()$PrintLevel,
                         message=NULL)
}

\arguments{
  \item{h}{a list of the
    same format as the output \code{\link{xswms2d}},
    see there for a description of the components.
  } 
 \item{trials}{number of trials to find a realisation of the plant
    positions (by a trial and error algorithm) where the distances have
    a distance of at least \code{$plants.mindist}.
  }
  \item{PrintLevel}{
     \itemize{
     \item \eqn{<}=0 : no messages are printed\cr
     \item 1 : failures of \code{optim} are announced
     \item 2 : occurences of partial failures
     \item 3 : the percentage of performed work is shown.
     \item \eqn{>}=4 : further tracing information
   }
   Note that if \code{PrintLevel>2} a running counter is shown that
   includes the printing of backspaces
   (\eqn{\mbox{\textasciicircum}}{^}H).
   The backspaces may have undesirable interactions with some few other
   R functions, e.g. 
   \command{\link[utils]{Sweave}}.
   See package \link[RandomFields]{RandomFields} for the default option
   \command{\link[RandomFields]{RFparameters}}\code{()$Print}.
  }
  \item{message}{
    function with one parameter (string). The function is used to plot messages. 
  }
}
\details{
  The simulation algorithm works currently as follows.
  \enumerate{
    \item
    For each plant type and according to \code{$plants.lambda} a Poisson
    random variable is drawn which gives the number of plants of this
    type.
    The starting points of the plants are created successively where the
    horizontal positions are uniformly distributed over the profile.
    The vertical position of the starting point is 1
    plus the vertical position of the first
    pixel that is not air at the drawn horizontal pixel.
    If the simulated horizontal distance has distance less than
    \code{$plants.mindist} to any other plant of this type, the
    horizontal distance is redrawn up to \eqn{\code{trials} - 1} times.
    If still not successful the last drawn position is taken.

    \item In the second step the root growth is simulated, for each
    plant separately and independently.
    For each plant the following steps are iterates until or root segments have
    stopped growing or the aimed total root length has been reached.
    In each step a list of active root segments is considered that is
    ordered with respect to some priority value \eqn{v}
    \enumerate{
      \item The root with the highest priority stops growing with
      probability \code{$stop.probab}.  If this happens this root segment
      is dropped from the list. This step is repeated until
      the list is void (then the simulation of the root system
      is terminated) or the growth does not stop.

      \item The list \eqn{p} of neighbours of the selected root
      is determined.  Out of the 8 potential neighbours only those are
      selected which do not leave the profile and that do not entre
      free air or a stone.
      In case \code{no.own.root} is \code{TRUE} the neighbouring pixel
      may also not be a already simulated segment of the current root system.

      \item
      The priority number of the remaining potential new root segments
      is a sum of the following values

      \enumerate{
	\item the direction change from previous segment to current to
	next is calculated (separately for x and y coordinate);
 	this value is the mean and the standard deviation (up to
	\code{$dir.ch.s}) of a Gaussian random variable.

	\item \code{$depth.bonus} or \code{$side.bonus}
	if the neighbouring pixel will be in one of the three
	pixels in vertical direction or in one of the two
	horizontal directions.

	\item \code{$rf.link} function value for \code{rf}
	that can be interpreted as water availability of soil quality.
      }
      
      \item if the distance to the previous knot is greater than
      \code{knot.mindist} then the current root segment is a knot
      with probility \code{$knot.probab}.

      The number of children is determined both by the number of
      available neighbouring pixels and the values of
      \code{shoots.4} and \code{shoots.3}.

      \item The lengths of the root segments to be added are calculated. 
      If then the aimed total length is exceeded, (biased) coins are
      thrown which of the neighbouring pixels are kept
      
      \item The remaining children/neighbouring pixels are added to the
      list of active root segments; The current, selected root segment
      is deleted from the list. 
      
      \item 
      The priority of root segments that have already been created in
      former steps are multiplied by \code{$age.bonus}.
    }
  }
}
\value{
  The function returns the input argument \code{h}
  except that \code{h$plants}, \code{h$plant.idx} and \code{h$rf}
  have been updated.

  \item{$plants}{is a list of simulated plants
    where each components contains the
    positions of the roots.  The latter is given by a matrix of 8 columns:
    \enumerate{
      \item horizontal coordinate of 
      a root segment.
      The coordinate is given in
      pixels. 
      \item vertical coordinate of a root segment,
      For the first pixel it is 1 plus the vertical position of the first
      pixel that is not air at the drawn horizontal pixel.
      \item index for the parent root pixel
      \item index where the subsequent root pixel is found; \code{NA} if
      it is a terminating root segment.  
      If the pixel is a knot with k children then
      the value gives the position of the first child, the others
      follow immediately.
      \item number of children in case it is a knot, 1 otherwise
      \item the number of preceding knots until and including this position
      \item aimed random total length of the roots if it is the first
      pixel, and the current distance to the surface along the roots,
      otherwise.
      \item distance to the preceding knot
    }
  }
  
  \item{$plant.idx}{
    vector or \code{NULL}.  \code{h$plant.idx[i]} gives the
    plant type (\eqn{1,\ldots,\code{}}{1,...,\code{plant.types}})
    of the \code{i}th simulated plant.}

  \item{$Root.RF}{
    \code{$Stone.RF} modified as follows:
    if \code{h$root[[i]]$rf.Kf} is \code{TRUE} then
    the locations of the roots of plant type \eqn{i} get
    the value of the function \code{h$Kf}.
  }
}

\references{  
  Schlather, M. and Huwe, B. (2005) Modelling inhomogeneous soils:
  Theory. \emph{In preparation.} %\bold{}, - .
  \cr
}

\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}

\seealso{ \code{\link{SoPhy}}, \code{\link{xswms2d}}}

\examples{
%  library(SoPhy, lib="~/TMP"); library(RandomFields)
% source("/home/schlather/article/R/SOPHY/SoPhy/R/simu.R")
## get a standard definition first
h <- xswms2d(xlim=c(1, 100), ylim=c(1, 100), step=1, new=NULL)
h$root[[1]]$plants.lambda <- 0.02
h <- simulate.horizons(h)
for (i in 1:5) {
  plotRF(create.roots(h), cex.pch=0.5, pch=16, root.col=1) 
  readline("Press return")
}
}
\keyword{spatial}

\name{create.waterflow}
\alias{create.waterflow}
\title{Preparation for swms2d}
\description{
 The function creates the input list for \command{\link{swms2d}}
}
\usage{
create.waterflow(h, Kat = 2, MaxIt = 20, hTab = c(0.001, 200), 
                 %Units = c(LUnit = 'cm', TUnit = 'sec',  
                 %             MUnit = '-', BUnit = '-'), 
                 dtMinMax = h$water$dtMinMax, DMul = c(1.1, 0.33),
                 atmadd=FALSE
                 )
}
\arguments{
  \item{h}{a list of the
    same format as the output \command{\link{xswms2d}}}
 % \item{Units}{vector of character.  Units for length, time, mass, and a
 %   further quantity, which is unclear}
  \item{Kat}{type of flow system:
    \itemize{
      \item 0 :  horizontal
      \item 1 :  axisymmetric (well or single drain)
      \item 2 :  vertical (profile)
    }
  }
  \item{MaxIt}{maximum number of iteration during any time step.}
  \item{hTab}{\code{c(hTab1, hTabN)}, interval of pressure heads
    within which a table of hydraulic properties is generated.}
  \item{dtMinMax}{\code{(dtmin, dtmax)}, minimum and maximum permitted
    time increment.}
  \item{DMul}{\code{c(dMul, dMul2)},
    \eqn{\code{dMul}\ge 1}{\code{dMul >= 1}},
    \eqn{\code{dMul2}\le 1}{\code{dMul2 <= 1}};
    if number of required iterations is less than 4 or greater than 6
    then the next time
    step is multiplied by \code{dMul} and \code{dMul2}, respectively.}
  \item{atmadd}{logical. If a segment contains several roots, the
    potential water uptake for all the roots is added if
    \code{atmadd=TRUE}. Otherwise the maximum potential water uptake
    is assigned to the segment.
  }
}
\details{
 The function is usually called internally; see Examples 
 for an explicite call.  The function performes several tasks
  \itemize{
    \item
    According to \code{h$water$red} the grid is thinned.

    \item
    Several checks are performed to insure that the users input parameters are
    sound.

    \item
    the boundary points and the numerical conditions are determined.

    \item 
    the root uptake is determined. This functionality
     is still in an experimental
    stage. 

    If several root segments belong to a
    single pixel then the boundary condition is taken that appears
    most. 
    If equality, then in the order of atmospheric condition, Neumann
    condition, and Dirichlet condition, the first one is taken as
    condition for the pixel. 

    In case of the atmospheric condition, the root
    uptake function is the minimal function of the standard form (given
    by P0, P2H, P2L, and P3), such
    that the values equal the maximal root uptake function
    at the considered pixel. 

    For the Neumann condition, the sum of the Neumann values is taken
    over all root segment at the pixel with Neumann condition. 
    
    For the Dirichlet condition, the minimum over all Dirichlet values
    is taken.

    \item the finite element mesh is created, where points that belong to
    the atmosphere or to a non-penetrable stone (\code{$stone$value=NA})
    are excluded.
  }
}
\value{
  The output is a list suitable as input list for \command{\link{swms2d}}.
  Additionally, the list has the following elements:

  \item{flux}{
    a \code{length(water.x)}\eqn{\cdot}{*}\code{length(water.y)}
    vector of logical values indicating
    which grid elements of the simulation
    grid are part of the finite element mesh.
  }
  \item{water.x}{the x coordinates of the thinned grid}
  \item{water.y}{the y coordinates of the thinned grid}}

\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/~schlather}}

\seealso{
         \command{\link{modify.horizons}},
         \command{\link{plotRF}},
         \command{\link{plotWater}},
         \command{\link{simulateHorizons}},
         \code{\link{SoPhy}},
         \command{\link{swms2d}},
	 \command{\link{xswms2d}}
       }

\examples{
%  library(SoPhy, lib="~/TMP"); library(RandomFields)
% source("/home/schlather/R/SOPHY/SoPhy/R/Sophy.R")
% source("/home/schlather/R/SOPHY/SoPhy/R/simu.R")
% source("/home/schlather/R/SOPHY/SoPhy/R/swms2d.R")
h <- xswms2d(xlim=c(1, 50), ylim=c(1, 30), step=1, new=NULL)
h$water$red <- 1
h$water$top.value <- -50   ## pressure head at surface
h$water$TPrint <- 10
h <- simulateHorizons(h)  ## simulation of the stochastic components
plotRF(h)
swms2d.in <- create.waterflow(h) ## rewrite h as a list readable by
##        swms2d; the following function can also directly be called 
##        with the argument h; then  create.waterflow is called
##        internally
str(swms2d.in)
swms2d.out <- swms2d(swms2d.in, iter.print=1) ## numerical simulation
if (is.character(swms2d.out)) cat(swms2d.out, "\n") else
plotWater(swms2d.out, what="H")
}

\keyword{spatial}% at least one, from doc/KEYWORDS

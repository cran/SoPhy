\name{flowpattern}
\alias{flowpattern}
\title{Simulation of flow pattern}
\description{
  Purely stochastic simulation of a water flow pattern as they appear in
  dye tracer experiments
}
\usage{
flowpattern(type=c('identical', 'unif', 'independent', 
                               'dependent', 'all'),
      length.profile=200, depth=100,  width.slice=1,     
      delta.x = depth * sqrt(x.var), delta.y = depth * sqrt(y.var),
      lambda.path=0.1, len.x, len.y,  grid=FALSE,

      x.name='whittle', x.var=1, x.v.scale=10, x.kappa=2, x.h.scale=1,
      y.name=x.name, y.var=x.var, y.v.scale=x.v.scale, y.kappa=x.kappa,
      y.h.scale=x.h.scale, unif.b=2,
           
      drop.distr=function(x) x * 80, drop.name='whittle',
      drop.scale=1, drop.kappa=2, drops=1,
           
      selected.dist = 2/3, front.factor=2, method = NULL,
      endpoint.tolerance=0, measure=function(x) x^2,
           
      simu.method=if (type \%in\% c("dependent", "all")) "TBM3",
      drop.simu.method = NULL, register=c(1,2),
      old.paths=NULL, PrintLevel=RFparameters()$Print, wait=FALSE,
      raw=FALSE, compress=TRUE, max.points = 5500000)
}
\arguments{
  \item{type}{character, one of \code{'identical'}, \code{'unif'},
    \code{'independent'}, \code{'dependent'}, \code{'all'}.
     Type of model. See details.
   }
  \item{length.profile}{positive integer. Length of the profile.}
  \item{depth}{positive integer. Depth of the profile.}
  \item{width.slice}{positive integer. Thickness of the slice. Usually
    a small number.}
  \item{delta.x}{enlargement of the simulation block to the left and
    right of the profile to avoid edge effects.}
  \item{delta.y}{enlargement of the thickness of slice in both
    directions to avoid edge effects.}
  \item{lambda.path}{intensity of the Poisson point process at the
    surface; if \code{grid=TRUE} and
    \code{len.x} or \code{len.y} are missing they are
    approximated by the number the number of points of the process.
    If \code{grid=TRUE} and \code{len.x} and \code{len.y} are not
    missing, the value of \code{lamda.path} is ignored.
  }
  \item{len.x}{integer.  Number of grid points in x-direction.  Only used
  if \code{grid=TRUE}.  Then it overwrites \code{n.path}}
  \item{len.y}{integer.  Number of grid points in y-direction.  Only used
    if \code{grid=TRUE}.  Then it overwrites \code{n.path}}
  \item{grid}{logical.  If \code{FALSE} the locations are uniformly and
    independently distributed.}
  \item{x.name}{character.  Name of the covariance function used to
    simulate the path increments.  Only used if \code{type} is different
    from \code{'unif'}, see %\command{\link[RandomFields]{PrintModelList}} or
    \command{\link[RandomFields]{CovarianceFct}}.
  }
  \item{x.var}{numeric.  Variance of the covariance function.}
  \item{x.v.scale}{numeric.  Scale in vertical direction.}
  \item{x.kappa}{Vector of additional parameters for the covariance
    function, if there are any.
  }
  \item{x.h.scale}{numeric.  Scale for the two horizontal directions.
    Only used if \code{type='dependent'} or \code{type='all'} -
    otherwise, the path are independent.
  }
  \item{y.name}{analogously to \code{x.name} for the y-direction}
  \item{y.var}{analogously to \code{x.var} for the y-direction}
  \item{y.v.scale}{analogously to \code{x.v.scale} for the y-direction}
  \item{y.kappa}{analogously to \code{x.kappa} for the y-direction}
  \item{y.h.scale}{analogously to \code{x.h.scale} for the y-direction}
  \item{unif.b}{numeric greater than; upper endpoint of the uniform
    distributon; relevant only if \code{type='unif'}; see Details.
  }
  \item{drop.distr}{The stained path length reached by a single drop is
    calculated as follows.  First a Gaussian random field is simulated
    with the parameters given below. Then the Gaussian random field is
    marginally transformed to uniform distribution by the
    Gaussian distribution function. Then the inverse of the required
    distribution function is applied.  \code{drop.distr} gives this
    inverse.  It should allow for vectors as input variables.
  }
  \item{drop.name}{Name of the covariance function for the horizontal
    dependence of the Gaussian random field that is used to model the
    stained path length.  Only used if \code{type='all'} - otherwise
    a Gaussian random field with pure nugget effect is used. See
    % \command{\link[RandomFields]{PrintModelList}} or
    \command{\link[RandomFields]{CovarianceFct}} for possible values.
  }
  \item{drop.scale}{numeric.  Scale in (horizontal) direction.}
  \item{drop.kappa}{Vector of additional parameters for the covariance
    function, if there are any.}
  \item{drops}{Number of drops running through each path. Then the stained
    part equals the longest stained part for all drops.
  }
  \item{selected.dist}{number in \eqn{(0,1)} or a vector of integers.
    \code{selected.dist} used for estimating the risk index, see
    \command{\link{risk.index}}.
  
    \code{selected.dist} gives the
    distances for which the form parameter of the pareto distribution
    is estimated.
    If \code{selected.dist} is a number in \eqn{(0,1)} the distances are
    \code{1:(1 + round(max(data[,1]-1) * selected.dist))}.
    Otherwise the integers are interpreted as indices for \code{data}.
  }
  \item{front.factor}{
    \code{front.factor} is
    used for estimating the risk index, see
    \command{\link{risk.index}}.

     The upper bound
    for the upper endpoint equals the \code{front.factor} times the
    largest distance for which at least stained pixel is observed.

    The value should not be changed if the algorithm is not understood.
  }
  \item{method}{character, 'fix.m' or 'optim.m'.
    \code{method} is
    used for estimating the risk index, see
    \command{\link{risk.index}}.
  
    If 'fix.m', the estimated number of path equals the maximum number of
    ob stained pixels below a given threshold. If 'optim.m' the number
    of paths is estimated by optimisation.

    If \code{method} is a vector containing both values, the risk index
    is calculated for both methods.

    If \code{method} is \code{NULL} nothing is estimated.
  }
  \item{endpoint.tolerance}{numeric.
    \code{endpoint.tolerance} is
    used for estimating the risk index, see \command{\link{risk.index}}

     If the shape
    parameter is negative then the distribution has a finite upper
    endpoint.  Hence, mathematically, the lowest upper end point of the
    Pareto distribution is given as the largest distance for which at
    least one stained pixel is observed. 
    For stability reasons and because the observed data might be a
    scale mixure of Pareto distribution it is advantageous to allow for
    some tolerance of the minimal upper end point.
    
    If \code{endpoint.tolerance} is positive then the lower threshold
    for the upper end point is the largest distance for which the
    number of observed stained pixels is larger than
    \code{endpoint.tolerance}.

    If \code{endpoint.tolerance} is negative then the lower threshold
    equals largest distance for which at least one stained pixel is
    observed minus the modulus of \code{endpoint.tolerance}.
  }
  \item{measure}{function that measures the distance between points; the
    function gets and returns a vector; used for estimating the risk
    index, see \command{\link{risk.index}}}
  \item{simu.method}{\code{NULL} or string;
    Method used for simulating
    the Gaussian random fields for the paths,
    see \link[RandomFields]{RFMethods}
   % or type \command{\link[RandomFields]{PrintMethodList}()}
    to get all options.
    If \code{model} is given as list then \code{method} may not be
    set if \code{model[[i]]$Method}, \eqn{i=1,3,..} is given, and vice
    versa

    \dQuote{TBM3} allows for some special effect:
    grid points and arbitrary points can taken from
    the same realisation by two subsequent calls of the function.
  }
  \item{drop.simu.method}{\code{NULL} or string; Method used for simulating
    the Gaussian random fields for the drops, see \code{simu.method}}
  \item{register}{vector of 2 components. The first compontent gives the
    register used in \command{\link[RandomFields]{GaussRF}} to simulate the
    paths, the 
    second component gives the register for the drops.}
  \item{old.paths}{\code{NULL} or \code{list}. If \code{list} it must have the
    same format as the return list of \code{flowpattern}.
    Then the flow pattern is recalculated where only the \code{drop.distr} 
    is changed.}
  \item{PrintLevel}{
    \itemize{
      \item <=0 : no messages are printed\cr
      \item >=1 : the higher the number the more information is given 
    }
  }
  \item{wait}{logical. If \code{PrintLevel>1} than the simulated profile
    is shown. If \code{wait=TRUE} the system waits for return after
    plotting.
  }
  \item{raw}{logical. If \code{grid=TRUE} then the simulated
    random fields for the x- and y-coordinates of the paths are
    also returned. Ignored otherwise.
  }
  \item{compress}{logical. If \code{grid=FALSE},
    then paths where no part is within the considered slice
    are deleted at an early stage.
    Ignored otherwise.
  }
  \item{max.points}{logical. If \code{grid=FALSE}, \code{type \%in\%
      c('dependent', 'all')} and \code{compress=TRUE},
    then the simulation of the
    the coordinates of the paths is done in several steps and
    ``compression'' is performed between the steps. 
  }
}

\value{
  it returns the following list
  \item{c.r}{empty list if \code{length(method)==0}.  Otherwise a list
    with \code{length(method)} components.  Each of them is the output
    of \command{\link{risk.index}}.
  }
  \item{intermediate}{list containing \code{xx}, a matrix of x-coordinates
    of the vertices of the paths, \code{yy}, a matrix of y-coordinates
    of the vertices of the paths, and \code{LEN}, the field of all the
    stained path lengths, after marginal transformation to uniform
    \eqn{U([0,1])} distribution.
  }
  \item{dist}{1:depth}
  \item{freq}{the number of stained pixels (after projection of the
    stained parts to the paths to a grid with grid length 1) as a
    function of \code{dist}}
  \item{i.x}{x-coordinates of the stained pixels}
  \item{i.d}{z-coordinates of the stained pixels}
  \item{input}{list of all input variables}
}
\details{
  \code{type} allows for the following values
  \describe{
    \item{\code{'identical'}}{All paths are identical, the
      increments of the path curve in x-direction are given by a
      Gaussian random field with parameters
      \code{x.name} \code{x.var},
      \code{x.v.scale} and \code{x.kappa}.
      The increments in y-direction by an independent Gaussian random
      field
      with parameters \code{y.name} \code{y.var},
      \code{y.v.scale} and \code{y.kappa}.
      The drops are i.i.d. according to \code{drop.distr}.
    }
    \item{\code{'unif'}}{For each path the angles of the segments
      to the horizontal plane are constant, such that the length of the
      paths are uniformly distributed
      \eqn{\sim}{~}\eqn{U([1,\code{unif.b}])}\eqn{*} [units of the
      vertical axis].
      The paths are variied only in the x-direction, but not in the
      y-direction.

      The drops are i.i.d.
    }
    \item{\code{'independent'}}{The paths and the drops
      are independently distributed.}
    \item{\code{'dependent'}}{
      The paths increments in x-direction are spatially dependent given
      by a Gaussian random field with parameters \code{x.name} \code{x.var},
      \code{x.h.scale} and \code{x.kappa}. Indeed, a three dimensional
      random field of geometrical anisotropy is generated, which
      includes the vertical direction with parameter \code{x.v.scale}.
      Analogously for the increments in y-direction. The x-y-coordinates
      of the sampling points of the random field are the starting points
      of the paths at the surface.

      The drops are i.i.d. according to \code{drop.distr}.}
    \item{\code{'all'}}{The paths are simulated as for
      \code{type='dependent'}. Here the drop distributions
      are spatially dependent, given by a two-dimensional
      Gaussian random field with
      parameters \code{drop.name}, \code{drop.scale} and
      \code{drop.kappa}.
    }
  }
  See the references for some background.
}
\references{
  \itemize{
    \item
    Schlather, M. and Huwe, B. (2005a)
    A stochastic model for 3-dimensional flow patterns in dye tracer
    experiments. \emph{J. Hydrol.} \bold{310}, 17-27.
    
    \item
    Schlather, M. and Huwe, B. (2005b)
    A risk index for characterising flow pattern in soils
    using dye tracer distributions.
    \emph{J. Contam. Hydrol.} \bold{79}, to appear.
  }
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/~schlather}}
\seealso{\command{\link[RandomFields]{RandomFields}}, \command{\link{plotFlow2d}}, \command{\link{plotFlow3d}}}
\examples{

############################################################
##                     simple example                    ##
############################################################
%fertig
%    library(SoPhy)
%library(RandomFields)
%    source("~/R/SOPHY/SoPhy/R/3dplot.R")

%%% muss wieder rueckgeaendert werden !!!
RFparameters(TBM3.lines=100, TBM3.linesimustep=0.3,
             TBM3.linesimufactor=0, Print=1)
fp3 <- flowpattern(type="all", lambda=0.1, length.profile=50,
                  depth=80, width=50, Print=20)
fp2 <- plotFlow3d(fp3, inf=c(1300, -500, 1/2000),
                  sun=c(100, -1500, 0))
plotFlow2d(fp2)
}
\keyword{spatial}

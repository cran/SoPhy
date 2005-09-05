\name{simulateHorizons}
\alias{simulateHorizons}
\title{
  Simulation of the stochastic part of the definition of the horizons 
}
\description{
  \command{simulateHorizons} simulates the Gaussian random fields that
  define the Miller-similar medium, the stones and the roots according
  to the definition of the horizons.
}
\usage{
simulateHorizons(h, first = 1, what = 'all',
                  PrintLevel = RFparameters()$Print,
                  message = function(s) {print(s)}, stone.trials=50)
}

%- maybe also `usage' foar other objects documented here.
\arguments{
  \item{h}{a list of the same format as the output of \command{\link{xswms2d}}}
  \item{first}{\code{first} gives the first horizon for which a
    Gaussian random field should be generated.  It is assumed that the
    preceding horizons have already been simulated.

    \code{first}
    greater than 1 is suitable in the following situation.
    Assume the random fields for all horizons of the profile
    have been generated.  Then the random field parameter of the \eqn{n}th
    horizon is changed.  For updating the simulated random fields it
    suffices now to call \command{simulateHorizons} by \code{first=n}}
  \item{what}{\code{what} is a second parameter afer \code{first} that
    is used to determine the starting point of the updating procedure:
    \describe{
      \item{'all'}{the simulations of the Gaussian random fields
	(according to \code{first}), the simulation of the stones and
	the simulation of the roots are performed}	       
      \item{'stone'}{only the simulations for the stones and the roots
	are performed}
      \item{'root'}{only the simulations for the roots are performed}
      \item{'randomfield'}{only the simulations for the Gaussian random
	fields are performed}
    }
  }
  \item{PrintLevel}{If non-positive nothing is printed.  The higher the
    number the more information is given on during the calculations.}
  \item{message}{function that has a string as argument.  \code{message}
    is called for printing error messages and for the current stage of
    the calculations.  \code{message} is called independently of the
    value of \code{PrintLevel}.  \code{message} is intended to be used
    for messages on graphical devices.
  }
  \item{stone.trials}{number of trials to find a realisation for the spatial
    distributions of the stones before \code{simulate.stones} gives up.}
}

%\details{ }

\value{
  List of the same format as the output of \command{\link{xswms2d}}.
  The following elements are changed in general

  \item{random.seed}{The random seed for the simulation of the Gaussian
    random fields.  It is changed when \code{what='all'} or
    \code{random.seed} has been \code{NULL}}

  \item{stone.random.seed}{The random seed for the simulation of the
    stones.  It is changed when \code{what\%in\% c('all','stone')} or
    \code{stone.random.seed} has been \code{NULL}}

  \item{root.random.seed}{The random seed for the simulation of the
    roots.  It is changed when \code{what=\%in\% c('all','stone','root')} or
    \code{root.random.seed} has been \code{NULL}}

  \item{RF}{the simulations of the Gaussian random fields}

  \item{Stone.RF}{the simulations of the profile including the values for the 
    stones, see also \command{\link{create.stones}}}

  \item{Root.RF}{the simulations of the profile including the root
    values if \code{$rf.Kf=TRUE}, see also \command{\link{create.roots}}}

  \item{rf.complete}{logical.  \code{rf.complete} is \code{FALSE} if the
    simulation has not been completed (because of an error).
  }
}

\note{
  As a side effect the global variable \code{.Random.seed} is reset.
}

\references{
  Schlather, M. (2001) Simulation and analysis of random fields
  \emph{R News} \bf{1} (2), 18-20.
}
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}


\seealso{\code{\link{SoPhy}}, \command{\link{xswms2d}}}

\examples{
%  library(SoPhy, lib="~/TMP"); library(RandomFields)
% source("/home/schlather/R/SOPHY/SoPhy/R/Sophy.R")
h <- xswms2d(xlim=c(1, 50), ylim=c(1, 50), step=1, new=NULL)
h$n <- 2            ## define a second
h$H2 <- h$H1
h$H2$type <- "H"    ## genuine horizon, not polygon 

## coordinates of the boundary segments between the horizons 
h$H2$points <- list(x=seq(1, 100, 5))             ## x coordinates
h$H2$points$y <- 40 + 20 * cos(h$H2$points$x / 5) ## y coordinates

## Gaussian random field of the second horizon is pure nugget effect:
h$H2$model <- list(model=list(list(model='nugget', var=0.5, scale=1)),
                   mean=0.5)
h$H2$stone$lambda <- 0.005  ## intensity of the stones in the 2nd horizon
h <- calculate.horizons(h)
for (i in 1:2) { 
  plotRF(simulateHorizons(h)) ##
  readline("Press return")
}

}
\keyword{spatial}

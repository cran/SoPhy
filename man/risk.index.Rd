\name{risk.index}
\alias{risk.index}
\title{Risk index for soil leaching}
\description{
  Estimation of the index of potential risk to groundwater
}
\usage{
risk.index(data, selected.dist=0.95, 
           selected.rate=cbind(c(0.5, 0.8), c(0.4, 0.9), c(0.3, 1.0)),
           weights=1, measure=function(x) x^2,
           method=c('fix.m', 'optim.m', 'ml'),
           min.neg.xi = -10, max.neg.xi = -0.1, max.pos.xi = 10,
           endpoint.tolerance = 0, front.factor = 2,
           # min.no.paths=max(data[, 2]),
           max.no.paths=10 * max(data[, 2]),
           PrintLevel=RFparameters()$Print, max.rate=TRUE)
}

\arguments{
  \item{data}{matrix of two columns. First column gives the distances
    (depths in the profile measured from the surface) and the second
    column the number of observed blue pixels.}
  \item{selected.dist}{scale or vector with values in \eqn{(0,1)}
    or a vector of integers. 
    Distances for which the form parameter of the pareto distribution
    is estimated; see Details. 
    If \code{selected.dist} is a number in \eqn{(0,1)} the distances are
    \code{1:(1 + round(max(data[,1]-1) * selected.dist))}. If it is a
    vector with values in (0,1), then the vector must have an even
    number of elements and pairs of elements are interpreted as intervals.
    
    Otherwise the integers are interpreted as indices for \code{data}.
  }
  \item{selected.rate}{vector of matrix of nrows.
    additionally to the indices given by
    \code{selected.dist}, the form parameter is estimated also for
    those distances where the corresponding number of observed relative stained
    pixels (w.r.t. to the maximum number of observed pixels)
    is within the interval given by the first column of
    \code{selected.rate}. 
    The risk index is calculated as the median of the estimated form
    parameters.
    
    In case no values are in the given interval of the first column
    the second column is considered etc, i.e. the first row should
    contain decreasing values and the second row increasing values.
  }
  \item{weights}{the estimation algorithm is based on a weighted least
    square algorithm; \code{weights} is usually either \eqn{1} or a vector
    of length \code{nrow(data)}.}
  \item{measure}{instead of the default least squares another distance
    function can be given.}
  \item{method}{the number of observed paths is a free parameter when
    fitting the Pareto distribution. It can either be set as the maximum
    number of stained pixels for the currently considered distances or
    depths (\code{'fix.m'}) or fitted within the optimisation algorithm
    (\code{'optim.m'}). Usually, it is not worth using the slower
    \code{'optim.m'} option.
    See also the Details.
  }
  \item{max.neg.xi}{optimisation parameter : largest negative
    value that is allowed as
    shape parameter of the Pareto distribution, i.e. a negative
    value close to 0.
  }
  \item{min.neg.xi}{optimisation parameter : smallest negative
    value that is allowed as
    shape parameter of the Pareto distribution}
  \item{max.pos.xi}{optimisation parameter : largest allowed 
    shape parameter of the Pareto distribution}
  \item{endpoint.tolerance}{optimisation parameter.  If the shape
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
  \item{front.factor}{optimisation parameter . The upper bound
    for the upper endpoint equals the \code{front.factor} times the
    largest distance for which at least one stained pixel is observed.
    
    The value should best not be changed.
 }
 \item{max.no.paths}{the number of paths is estimated as nuisance
   parameter when estimating the risk index; 
   \code{max.no.paths} give the upper bound for the nuisance
   parameter in the optimisation.}
 \item{PrintLevel}{
    The higher the value of \code{PrintLevel} the more tracing
    information is given. Up to value 1, no information is given.
%     \itemize{
%     \item \eqn{<}=0 : no messages are printed\cr
%     \item 1 : failures of \code{optim} are announced
%     \item 2 : occurences of partial failures
%     \item 3 : the percentage of performed work is shown.
%     \item \eqn{>}=4 : further tracing information
%   }
   Note that if \code{PrintLevel>=2} a running counter is shown that
   includes the printing of backspaces
   (\eqn{\mbox{\textasciicircum}}{^}H).
   The backspaces may have undesirable interactions with some few other
   R functions, e.g. \command{\link[utils]{Sweave}}.
   See package \link[RandomFields]{RandomFields} for the default option
   \command{\link[RandomFields]{RFparameters}}\code{()$Print}.
 }
 \item{max.rate}{logical. If \code{TRUE} then 
   the lines for which \eqn{m(D) / m(0)} is in \code{selected.rate}
   are used to calculate the final risk index.
   Here \eqn{m(D)} gives the maximum of \eqn{p(d)}, \eqn{d=D, D+1,
     \ldots}{d=D, D+1, ...} where  \eqn{p(d)} the number of stained
   pixels in depth \eqn{d}.
   If \code{FALSE} then the criterion \eqn{m(D) / m(0)} is replaced
   by \eqn{p(D) / p(0)}.
   }
}
\details{
  Denote by \eqn{f(d)} the number of blue pixels registered at depth
  \eqn{d} (or distances from the soil surface). 
  Then, the risk index is by definition a shape parameter of
  \eqn{f(d)} for large distances \eqn{d}.  Since the term \emph{large} cannot be 
  defined precisely, the shape parameter is calculated for the
  function values \eqn{f(d)} for distances \eqn{d \ge d_i}{d>=d_i}
  and several fixed starting distances \eqn{d_i}. 
  The distances \eqn{d_i} are given by \code{selected.rate}.
  (The approach is similar to that for analyzing extremal events.)

  The selection criterion \eqn{m(D) / m(0)} is \bold{always} based
  on \code{ method='fix.m'}, whatever \code{method}
  is chosen to estimate \eqn{\xi}{xi}.
}
\value{
  list of the following components
  \item{par}{matrix of estimated parameters; first row: risk index;
    second row: scale parameter; third row: estimated maximum number of
    paths \eqn{m(D)} except \eqn{m(0)} that is given by \code{max.freq}
    and is always set to the maximum number of pixels.
    forth row: \eqn{D} (\code{sel.dist}).
  }
  \item{data}{the input data except for some reordering}
  \item{weights}{the input weights expcept for some reordering}
  \item{selected.dist}{the selected distances in form of indices
    (in clear text, in case
    they were given in form of a real value in \eqn{(0,1)}).}
  \item{selected.rate}{range of the selected number of stained pixels}
  \item{sel.rate}{index set for the \code{data} where the observed
    number of stained pixels are within \code{selected.rate}}
  \item{sel.dist}{the index set containing \code{selected.dist} and
    \code{sel.rate}}
  \item{max.freq}{maximum number of observed stained pixels}
  \item{values}{the minimal least squares values}
  \item{method}{the input parameter \code{method}}
  \item{measure}{the input parameter \code{measure}}
  \item{raw.risk.index}{risk index calculated as median of the
    estimated form parameters for \code{selected.rate}
  }
  \item{risk.index}{
    the median is calculated only for values greater than \eqn{0.999
      \code{min.neg.xi}} and less than \eqn{0.999
      \code{max.pos.xi}}
  }
}

\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/~schlather}}

\references{
  Embrechts, P., Klueppelberg, C. and Mikosch, T. (1997)
  \emph{Modelling Extremal Events.} Berlin: Springer.

  Schlather, M. and Huwe, B. (2005)
  A risk index for characterising flow pattern in soils
  using dye tracer distributions.
  \emph{J. Contam. Hydrol.} \bold{79}, to appear.
 }

\seealso{
  \code{\link{SoPhy}},
  \command{\link{analyse.profile}},
  \command{\link{xswms2d}}
}

\examples{
%  library(SoPhy, lib="~/TMP"); source("~/R/SOPHY/SoPhy/R/extremal.coeff.R")
sample.depth <- 1 : 100
d <- rexp(1000, 1/25)
freq <- numeric(length(sample.depth))
for (i in 1:length(sample.depth)) freq[i] <- sum(d>=sample.depth[i])
cr <- risk.index(cbind(sample.depth, freq),
                 selected.rate=c(0.95, 0.9), 
                 endpoint.tolerance=20,  method="fix.m")
cr$risk.index ## the true value is 0
}
\keyword{spatial}

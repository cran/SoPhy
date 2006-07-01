\name{tortuosity}
\alias{tortuosity}
\title{Variability of the tortuosity}
\description{
  The function estimates the ratio between the mean tortuosity
  and the smallest tortuosity from profiles of Brilliant Blue tracer
  experiments
}
\usage{
  tortuosity(depth, freq,range.depth, range.xi, len.xi=100,
             lower.bound.m = 1.005, tolerance.m=0.01, PrintLevel=0) 
}

\arguments{
  \item{depth}{depth in pixels}
  \item{freq}{number of stained pixels found in depth \code{depth}}
  \item{range.depth}{range of depths for which \eqn{1-H^*} is fitted } 
  \item{range.xi}{range of values of \eqn{\xi}{xi} for which \eqn{1-H^*}
    is fitted} 
  \item{len.xi}{\code{range.xi} and \code{len.xi} define the grid
    \code{seq(min(range.xi), max(range.xi), len=len.xi)} for which
    \eqn{1-H^*} is fitted.
  }
  \item{lower.bound.m}{lower bound for \eqn{m} used in the numerical
    optimisation}
  \item{tolerance.m}{
    non-negative value. If  \code{tolerance.m=0} only
    the optimal \eqn{\xi}{xi} is used to estimate \eqn{m}.
    If \code{tolerance.m>0} then all estimated \eqn{m(D, \xi)}{m(D, xi)}
    are considered for which the span of \eqn{m} is less than or equal to
    (1 + \code{tolerance.m}) times the minimal span value.
  }
  \item{PrintLevel}{integer. If the value is greater than 1 then tracing
  information is given}
}
%\details{
%}
\value{
  \command{tortuosity} returns a list of the following elements:

  \item{xi}{vector of length \code{len.xi}. It contains the values
    of \eqn{\xi}{xi} used in the algorithm}  
  \item{raw.m}{matrix with \code{len.xi} rows. The number of colums
    equals the number of depth values that are between the given
    bounds \code{range.depth}. It contains
    the fitted values of \eqn{m}}
  \item{span.xi}{vector of length \code{len.xi} contains span of all
    \code{raw.m} values for the respective value of \eqn{\xi}{xi}}
  \item{span.m}{vector of length \code{len.xi}.
    It contains the difference between the largest and the smallest
    value of \code{raw.m} for given \eqn{\xi}{xi}
  }
  \item{opt}{list of the optimal values: \code{xi}, \code{m}, \code{s},
    \code{mspan} (the \code{span.xi} found at \code{xi}) and \code{D}
    (threshold depth)
  }
  \item{input}{list of two vectors: \code{freq} and \code{depth}}
  \item{fitted}{list of two vectors that give the best fitted
    curve: \code{depth} and \code{p}; the latter cannot be given
    if \eqn{\code{tolerance.m}>0}}
}

\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/~schlather}}
\references{
  Schlather, M. and Huwe, B. (2005)
  A characteristic for the variability of tortuosity.
  \emph{Submitted.}
}

\examples{
%  library(SoPhy, lib="~/TMP")
data(F04)
path <- paste(system.file(package='SoPhy'),"tracer", sep="/")
F04$name <- paste(path, "F04.G.tif", sep="/")
m <- analyse.profile(F04, estimate.all=FALSE, method="fix.m",
                     selected.rate=c(0, 0), selected.dist=NULL,
                     interactive=FALSE, Print=0)
stained <- m$fct(m$picture, m$param)
freq <- m$absfr
truedepth <- m$r.i$data[, 1] + ncol(stained) - m$loc[[1]]$y[2]
tol = 0

mt <- tortuosity(truedepth, freq,  range.depth=c(85, 97),
                 range.xi=c(-2.5, -0.6), Print=2, tol=tol)
str(mt$opt)

# fig. 3A in Schlather and Huwe (2005)
matplot(mt$xi, mt$raw.m, xlab=expression(xi), ylab="m",
        type="l", col=1, lty=1)     

# fig. 3B in Schlather and Huwe (2005)
plot(mt$xi, mt$span.m, xlab=expression(xi), ylab="span m",
     type="l", col=1, lty=1)
for (cex in c(9, 9.5, 10))
  points(mt$opt$x, mt$opt$mspan, cex=cex, col="darkgrey")

# fig. 2 in Schlather and Huwe (2005)
plot(mt$input$depth, mt$input$freq, xaxs="i", yaxs="i",
     pch=16, ylab="stained pixels", xlab="depth [pixels]",
     col="blue", cex=2)
lines(mt$fitted$depth, mt$fitted$p, lwd=3, col=2)
}
\keyword{spatial}

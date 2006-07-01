\name{Pareto}
\alias{Pareto}
\alias{pareto}
\alias{dpareto}
\alias{ppareto}
\alias{qpareto}
\alias{rpareto}
\title{The Pareto Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the Pareto distribution with form parameter \code{xi}
  and scale parameter \code{s}
}
\usage{
ppareto(q, xi, s=1, lower.tail=TRUE)

dpareto(x, xi, s=1)

qpareto(p, xi, s=1, lower.tail=TRUE)

rpareto(n, xi, s=1)
}
\arguments{
  \item{xi}{form parameter of the Pareto distribution}
  \item{s}{scale parameter of the Pareto distribution}
  \item{lower.tail}{logical; if \code{TRUE} (default),
    probabilities are \eqn{P[X <= x]}, otherwise, \code{P[X > x]}}
  \item{q}{vector of quantiles}
  \item{x}{vector of quantiles}
  \item{p}{vector of probabilities}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
          taken to be the number required.
}
}
\value{
  \code{dpareto} gives the density, \code{ppareto} gives the distribution
     function, \code{qpareto} gives the quantile function, and \code{rpareto}
     generates random deviates.
}
\references{
  Embrechts, P., Klueppelberg, C. and Mikosch, T. (1997)
  \emph{Modelling Extremal Events}. Springer: Berlin.
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/~schlather}}
\keyword{distribution}

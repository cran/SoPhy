\name{plotRGB}
\alias{plotRGB}
\title{Plotting RGB figures}
\description{
  The function plots an RGB array
}
\usage{
plotRGB(picture, x=1:dp[1], y=if (reverse) dp[2]:1 else 1:dp[2],
         reverse=TRUE, cex.axis=1, ...)
}
\arguments{
   \item{picture}{three dimensional array where the third dimension must
have length 3 or 4; the forth component is ignored if there is any}
   \item{x}{x-coordinates; the length must match the length of the first
dimension of \code{picture}}
   \item{y}{y-coordinates; the length must match the length of the first
dimension of \code{picture}; see details}
   \item{reverse}{logical; see details}
   \item{cex.axis}{graphical parameter, see \command{\link{par}}}
   \item{...}{additional graphical parameters for \command{\link[base]{image}}}
}
\value{
 NULL
}
\details{
 \code{reverse=TRUE} is for the plot of soil profiles, where the
  numbering of the y-axis start at the top and is counted positively
  downwards.  

  \code{reverse=FALSE}
  The usual presentation of the vertical axis.
 
  The functions checks whether \code{y} is in decreasing or increasing
  order, and plots the y axis accordingly. 
}
\seealso{\command{\link{read.picture}}, \command{\link{write.picture}}}

\examples{
 %library(SoPhy)
fig <- read.picture(paste(system.file(package='SoPhy'),
                    'tracer', 'K06', sep="/"))
plotRGB(fig)

}
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}
\keyword{hplot}

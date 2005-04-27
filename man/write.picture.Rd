\name{write.picture}
\alias{write.picture}
\title{RGB to TIF}
\description{
  The function writes an RGB-array to a TIF file. 
}
\usage{
write.picture(picture, tif)
}
\arguments{
  \item{picture}{three dimensional array where the third dimension must
    have length 3 or 4}
  \item{tif}{character.  Name of the tif file.}
}
\value{
  None.
}
\seealso{
  \code{\link{plotRGB}}, \code{\link{read.picture}}
}
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}
\keyword{IO}
\keyword{file}


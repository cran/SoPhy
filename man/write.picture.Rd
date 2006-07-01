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
  \command{\link{plotRGB}}, \command{\link{read.picture}}
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/~schlather}}
\keyword{IO}
\keyword{file}


\name{read.picture}
\alias{read.picture}
\title{Reading TIF pictures and conversion to RGB matrix}
\description{
  The function reads in TIF formated pictures and returns an RGB-A
  matrix, where A is the associated alpha matting information, if there
  is any (and identically 255 otherwise). If the format is different
  from TIF, it is first convert-ed to TIF.
}
\usage{
read.picture(picture, extensions=c("tif", "tiff",
                      if (.Platform$OS.type=="unix") c("gif", "jpeg")),
             PrintLevel=RFparameters()$Print, tmp.dir=".")
}
\arguments{
  \item{picture}{name of the graphic file, see details}
  \item{extensions}{standard extensions for the graphic file}
  \item{PrintLevel}{
    \itemize{
      \item \eqn{\le}{<=}0 : no messages are printed\cr
     \item \eqn{>0} : the higher the number the more information is given.
    }
  }
  \item{tmp.dir}{NULL or character.
    The variable is used if the picture file does not end with \sQuote{tif} or
    \sQuote{tiff}. Then the picture file is \emph{convert}ed into a tif
    file before being read in. The storage path is given by
    \code{tmp.dir}. An existing tif file in \code{tmp.dir} is not
    overwritten, but this file used instead and a warning is given.
  }
}
\value{
  three-dimensional array, where the first two dimensions are the
  dimensions of the graphic, and the third dimension has length 4, for
  red, green, blue and alpha.
}
\details{
  
  The procedure of finding the graphic is the following.
  \itemize{
    \item  if the file is given with one of the standard extension, it is
    assumed that this file is wanted

    \item else it is assumed that only the base name was given. Then,
    as above, the standard extensions are added sequentially.

    \item else the algorithm fails
  }
  If the file is neither TIF nor TIFF it is first convert-ed.
}
\seealso{\command{\link{plotRGB}, \command{\link{write.picture}}}
}
\examples{
fig <- read.picture(paste(system.file(package='SoPhy'),
                    'tracer', 'K06', sep="/"), Pr=3)
plotRGB(fig)
}
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}
\keyword{IO}
\keyword{file}




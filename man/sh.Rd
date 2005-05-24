\name{sh}
\alias{sh.jh}
\alias{sh.jch}
\title{Sample code}
\description{
  The functions provide the code for the figures and simulations
  in the publications of Schlather and Huwe
}
\usage{
sh.jch(input=NULL, dev=2, pspath="./", txt.result.dir = "txt/",
       simu.path="simu/", final=TRUE, PrintLevel=0, bw=TRUE,
       readlines=TRUE)
sh.jh(input=NULL, dev=2, pspath="./", final=TRUE, PrintLevel=0,
       readlines=TRUE, low.resolution=TRUE)
}
\arguments{
  \item{input}{vector of values for automatic choice of the submenues}
  \item{dev}{device, see \command{\link[RandomFields]{Dev}}}
  \item{pspath}{path used for all postscript files in \command{sh.jh}
    and for all postsript files in \command{sh.jch} except those based
    on the extended simulation study
  }
  \item{txt.result.dir}{path used for non-graphic results
  }
  \item{simu.path}{path used for the postscript files based
    on the extended simulation study
  }
  \item{final}{logical. If \code{TRUE} all the parameters
    are identical to those used in the papers; then at least 1GB of
    memory on a Linux machine is needed. If \code{FALSE}
    some few parameters are changed so that the simulations
    run faster and less memory is needed.
  }
  \item{PrintLevel}{If \eqn{\code{PrintLevel}>1} some information on
    the course of the algorithm is printed
  }
  \item{bw}{logical. if \code{TRUE} then all graphics will be printed
    in black and white, only.
  }
  \item{readlines}{logical. Only used if \code{is.numeric(dev)}. If
    \code{readlines=FALSE} then the system
    waits a second before the next plot is calculated or shown.
    Otherwise the system waits for return.
  }
  \item{low.resolution}{logical. Figures of low resolution are created.
    This takes time, but reduces size by about factor 10. This feature
    only works on unix/linux systems.
  }
}
\details{
  Additionally to the figures in Schlather and Huwe (2004a,b),
  sh.jh allows for the simulation of profiles given in
  Crestana and Posadas (1998) and
  Schwartz et al. (1999).
}
\value{
  \code{NULL}
  
  Depending on the value of \code{dev}, see
  \command{\link[RandomFields]{Dev}}, postscript or pdf files are created.
}
\seealso{
  \code{\link{SoPhy}}, \command{\link{plotFlow2d}}, \command{\link{plotFlow3d}}
}
\references{
  Schlather, M. and Huwe, B. (2005a)
  A stochastic model for 3-dimensional flow patterns in dye tracer
  experiments. \emph{J. Hydrol.} \bold{In Press}, .
    
  Schlather, M. and Huwe, B. (2005b)
  A risk index for characterising flow pattern in soils
  using dye tracer distributions
  \emph{Submitted to J. Contam. Hydrol.} % \bold{}, .

  Crestana, S. and  Posadas, A.N.D. (1998)
  k 2-d and 3-d fingering in unsaturated soils investigated by fractal
  analysis, invasion percolation modeling and non-destructive image processing.
  In Baveye, P., Parlange, J.Y., and Stewart, B.A., editors,
  \emph{Fractals in Soil Science}. Boca Raton: CRC Press.

Schwartz, R.C., McInnes, K.J., Juo, A.S.R., and Cervantes, C.E. (1999)
The vertical distribution of a dye tracer in a layered soil.
\emph{Soil Sci.}, \bold{164}, 561--573.
}
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}
\examples{

% library(SoPhy, lib="~/TMP"); 
% source("/home/schlather/R/SOPHY/SoPhy/R/sh.R")

  final <- FALSE # final <- TRUE
  printlevel <- 1 + !interactive()
  sh.jh(input=if (!interactive()) c(1, 2, 5, 0), final=final,
       Print=printlevel)
  sh.jch(input=if (!interactive()) c(1:14, 0), final=final,
         Print=printlevel)

}
\keyword{spatial}


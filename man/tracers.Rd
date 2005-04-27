\name{Tracer experiments}
\alias{BrilliantBlue}
\alias{brilliantblue}
\alias{tracer}
\alias{K06}
\alias{F06}
\alias{F04}
\title{Brilliant Blue tracer profile}
\usage{
   path <- paste(system.file(package='SoPhy'), 'tracer', sep="/")
   F04  <- read.picture(paste(path, "F04", sep="/"))
   F04g <- read.picture(paste(path, "F04.G", sep="/"))
   F06  <- read.picture(paste(path, "F06", sep="/"))
   F06g <- read.picture(paste(path, "F06.G", sep="/"))
   K06  <- read.picture(paste(path, "K06", sep="/"))
   K06g <- read.picture(paste(path, "K06.G", sep="/"))
}
\description{
  Digitised pictures of profiles that have been
  excavated after Brilliant Blue tracer
  experiments.
  The picture were taken by the soil physics group at ETH Zurich, see
  Source and References below.
}
\format{
  3-dimensional array, where the last dimension codes red, green, blue,
  and alpha.
}
\source{
  The data were collected by the soil physics group at ETH Zurich,
  \url{http://www.ito.umnw.ethz.ch/SoilPhys/Fliessmuster/}
}
\references{
  Data
  \itemize{
    \item  Flury, M., Leuenberger, J., Studer, B. and Jury, W.A.
    and  Fluehler, H. (1992-94)
    \emph{Prafeerentielle Fliesswege in Ackerboeden - ein screening-Test unter
      Feldbedingungen.}
    \url{http://www.ito.umnw.ethz.ch/SoilPhys/Fliessmuster/projekt_flury.html}

    \item von Albertini, N., Laeser, H.P., Leuenberger, J. and
    H. Fluehler (1994-97)
    \emph{Strukturerholung eines verdichteten Ackerbodens unter einer
      mehrjaehrigen Kunstwiese.}
    \url{http://www.ito.umnw.ethz.ch/SoilPhys/Fliessmuster/projekt_vonalbertini.html}
  }

  Papers
  \itemize{
    \item von Albertini, N.,  Leuenberger, J., Laeser, H.P. and
    H. Fluehler (1995)
    Regeneration der Bodenstruktur eines verdichteten Ackerbodens unter
    Kunstwiese 
    \emph{Bodenkundliche Gesellschaft der Schweiz, Dokument} \bold{7}, 10-16
     
    \item Flury, M. and  Fluehler, H. (1994)
    Brilliant Blue FCF as a Dye Tracer for Solute Transport Studies - A
    Toxicological Overview 
    \emph{J. Environm. Quality} \bold{23}, 1108-1112    

    \item Flury, M. and Fluehler, H. (1995)
    Tracer characteristics of Brilliant Blue FCF. \emph{Soil
      Sci. Soc. Am. J.} \bold{59}, 22-27.

    \item Flury, M., Fluehler, H., Jury, W.A. and  Leuenberger,
    J. (1994)
    Susceptibility of Soils to Preferential Flow of Water: A Field
    Study.
    \emph{Water Resour. Res.} \bold{30}, 1945-1954

    \item German-Heins, J. and  Flury, M. (2000)
    Sorption of Brilliant Blue FCF in soils as affected by pH and
    ionic strength. \emph{Geoderma} \bold{97}, 87-101.

    \item Schlather, M. and Huwe, B. (2004)
    A risk index for characterising flow pattern in soils
    using dye tracer distributions
    \emph{Submitted to J. Contam. Hydrol.} % \bold{}, .
  }  
}
\examples{
% library(SoPhy)
%    source("/home/schlather/article/R/SOPHY/SoPhy/R/rgb.R")
##########################################################################
##  Note: these data can currently be treated only under unix systems!  ##
##########################################################################
#ifdef unix
path <- paste(system.file(package='SoPhy'), 'tracer', sep="/")
for (nm in c("F04", "F06", "K06")) {
   cat("\n", nm, ":\nOriginal picture. ", sep="")
   p <- read.picture(paste(path, "/", nm,  sep=""))
   plotRGB(p)
   readline("Press return.")
   p <- read.picture(paste(path, "/", nm, ".G", sep=""))
   plotRGB(p)
   readline("Picture after first preparation. Press return.")
}
#endif
}
\keyword{datasets}

\name{SoPhy}
\alias{SoPhy}
\alias{SOPHY}
\alias{SoilPhysics}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Soil physics tools}
\description{
  This package provides a collection of tools that might be useful in 
  soil physics. It is a first version of work in progress.
 }

\details{
  Currently, the following functionalities are provided:

  \itemize{
    \item \bold{Physical modelling}
    \itemize{
      \item \code{\link{swms2d}} : R interface to SWMS\_2d, an algorithm by
      Simunek et. al based on finite elements to simulate water flux and
      solute transport in soil
      \item \code{\link{xswms2d}} : interactive plot to define
      soil structure and to simulate water flux and solute transport
      \item various functions for the modified simulation of the water flux or
      the solute transport in a profile; see
      \code{\link{calculate.horizons}}
      for an example
    }
    \item \bold{Stochastic modelling}
    \itemize{
      \item \code{\link{flowpattern}}
      stochastic simulation of dye tracer profiles
      %  \item \code{\link{invasion}} : simulation of percolation
      %   invasion models
    }
    \item \bold{Estimation of the risk index}
    \itemize{
      \item \code{\link{analyse.profile}} : analysis of brilliant blue
      experiments w.r.t. the contamination risk index
      \item \code{\link{BrilliantBlue}} : examples of profile pictures
      \item \code{\link{Pareto}} : the Pareto distribution
      \item \code{\link{risk.index}} : esimation of the risk index 
    }
    \item \bold{Graphics}
    \itemize{
      \item \code{\link{my.legend}} : legend giving a colour spectrum
      \item \code{\link{quader}} : tools for 3d perspective plot
      \item \code{\link{read.picture}}, \code{\link{plotRGB}} :
      Conversion of tiff pictures to RGBA array and plotting RGB figures
    }
    \item \bold{Code used in publications by Schlather and Huwe}
    \itemize{
      \item \code{\link{sh.jh}} : Schlather, M. and Huwe, B. (2004a)
      \item \code{\link{sh.jch}} : Schlather, M. and Huwe, B. (2004b)
    }
  } 
}


\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}

  The Fortran code has been published as \sQuote{swms\_2d} 
  by Jirka Simunek, T. Vogel and Martinus Th. van Genuchten, 
  \url{http://www.ussl.ars.usda.gov/MODELS/MODELS.HTM};
  the file ORTHOFEM.f was written by Edward A. Sudicky  
  and Carl A. Mendoza,
  based on code by Frank W. Letniowshi and contributions by Rene Therrien,
  and modified by Jirka Simunek.
}

\references{
  Schlather, M. (2005)
  SOPHY - a package for the simulation and the analysis of water flux in soil
  \emph{In preparation}% \bold{}, .

  Schlather, M. and Huwe, B. (2004a)
  A stochastic model for 3-dimensional flow patterns in dye tracer
  experiments. \emph{Submitted to J. Hydrol.}% \bold{}, .
    
  Schlather, M. and Huwe, B. (2004b)
  A risk index for characterising flow pattern in soils
  using dye tracer distributions
  \emph{Submitted to J. Contam. Hydrol.} % \bold{}, .
   
  Schlather, M. and Huwe, B. (2005)
  Modelling inhomogeneous soils: Part I, Theory \emph{In preparation}% \bold{}, .

  Schlather, M. and Huwe, B. (2005)
  SOPHY - an interactive programme for water flow modelling
  \emph{In preparation}% \bold{}, .

  Schlather, M., Zeilinger, J. and Huwe, B. (2005)
  Modelling inhomogeneous soils: Part II, Applications
  \emph{In preparation}% \bold{}, .
 
  Simunek., J., Vogel, T., and van Genuchten, M.Th. (1994)
  \emph{The SWMS\_2D code for simulating water flow and solute transport
    in two-dimensional variably saturated media, Version 1.21.}
  Research Report No. 132, 197 p.,
  U.S. Salinity Laboratory, USDA, ARS, Riverside, California.
}

\section{Acknowledgement}{
  The author thanks Rien van Genuchten, Carl A. Mendoza, Rene Therrien,
  and Edward Sudicky for kindly changing the former copyright to
  the GNU copyleft licence, hence making their code directly available
  for this package.

  The work has been supported by the German Federal Ministry of Research
  and Technology (BMFT) grant PT BEO 51-0339476C during 2000-03.
}

\keyword{spatial}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line



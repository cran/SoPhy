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
      \item \command{\link{swms2d}} : R interface to SWMS\_2d, an algorithm by
      Simunek et. al based on finite elements to simulate water flux and
      solute transport in soil
      \item \command{\link{xswms2d}} : interactive plot to define
      soil structure and to simulate water flux and solute transport
      \item various functions for the modified simulation of the water flux or
      the solute transport in a profile; see
      \command{\link{calculate.horizons}}
      for an example
    }
    \item \bold{Stochastic modelling}
    \itemize{
      \item \command{\link{flowpattern}}
      stochastic simulation of dye tracer profiles
      %  \item \command{\link{invasion}} : simulation of percolation
      %   invasion models
    }
    \item \bold{Estimation of the risk index}
    \itemize{
      \item \command{\link{analyse.profile}} : analysis of Brilliant Blue
      tracer experiments w.r.t. the contamination risk index
      \item \command{\link{BrilliantBlue}} : examples of profile pictures
      \item \command{\link{Pareto}} : the Pareto distribution
      \item \command{\link{risk.index}} : esimation of the risk index 
    }
    \item \bold{Estimation of the variability of the tortuosity}
    \itemize{
      \item \command{\link{tortuosity}} : estimation of tortuosity
      from  Brilliant Blue tracer
      experiments
    }
    \item \bold{Graphics}
    \itemize{
      \item \command{\link{my.legend}} : legend giving a colour spectrum
      \item \command{\link{quader}} : tools for 3d perspective plot
      \item \command{\link{read.picture}}, \command{\link{plotRGB}} :
      Conversion of tiff pictures to RGBA array and plotting RGB figures
    }
    \item \bold{Code used in publications by Schlather and Huwe}
    \itemize{
      \item \command{\link{sh.jh}} : Schlather, M. and Huwe, B. (2005a)
      \item \command{\link{sh.jch}} : Schlather, M. and Huwe, B. (2005b)
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
  Schlather, M. (2006)
  SOPHY - a package for the simulation and the analysis of water flux in soil
  \emph{In preparation}% \bold{}, .

  Schlather, M. and Huwe, B. (2005a)
  A stochastic model for 3-dimensional flow patterns in dye tracer
  experiments. \emph{J. Hydrol.} \bold{310}, 17-27.
    
  Schlather, M. and Huwe, B. (2005b)
  A risk index for characterising flow pattern in soils
  using dye tracer distributions
  \emph{J. Contam. Hydrol.} \bold{79}, to appear.
   
  Schlather, M. and Huwe, B. (2006)
  Modelling inhomogeneous soils: Part I, Theory \emph{In preparation}% \bold{}, .

  Schlather, M. and Huwe, B. (2006)
  SOPHY - an interactive programme for water flow modelling
  \emph{In preparation}% \bold{}, .

  Schlather, M., Zeilinger, J. and Huwe, B. (2006)
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

  The work has been financially
  supported by the German Federal Ministry of Research
  and Technology (BMFT) grant PT BEO 51-0339476C during 2000-03.
}

\section{Installation}{

  LINUX: 
  The tiff-library is avalailable from
  \url{http://www.remotesensing.org/libtiff}.
  SuSe-users install libtiff and libtiff-devel.
  
  WINDOWS: 
  The binaries are available from
  \url{gnuwin32.sourceforge.net/packages/tiff.htm}.
  Before starting the R session, add in the 
  system variable \$PATH the path to the
  libtiff3.dll file; the standard path to be added 
  is C:\\Programme\\GnuWin32\\bin.
  Further, make sure that zlib1.dll and jpeg62.dll
  are present.

%  After installation copy the files in GnuWin32/lib and GnuWin32/include
%  to MinGW
}
     
}

\keyword{spatial}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line



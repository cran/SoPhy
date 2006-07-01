\name{swms2d}
\alias{swms2d}
\alias{SWMS2D}
\title{SWMS2D function call}
\description{
  Modelling of water flow and solute transport by SWMS2D
}
\usage{
swms2d(d, max.iteration = 1e+05, iter.print = max.iteration, ShortF = TRUE, 
       message = NULL, breakpoint = 1e+10, intermediate.result = NULL)
}
\arguments{
  \item{d}{list; see details}
  \item{max.iteration}{maximum number of incremental time steps}
  \item{iter.print}{number of incremental time steps after which a short
    message of status is given}
  \item{ShortF}{logical. If \code{TRUE} information is provided only at
    selected time points given by \code{d$TPrint}}
  \item{message}{\code{NULL} or a function that takes a scalar time
    argument and returns a logical value.
    The function is called in case \code{max.iteration} is exceeded.
    The function may allow the user to decide on the basis of the current
    time whether the simulation should be continued.}
  \item{breakpoint}{number of incremental time steps after which
    \code{intermediate.result} is called}
  \item{intermediate.result}{\code{NULL} or a function with two
    parameters: the first one gives the time, the second is a \eqn{(6 \times
      \code{NumNP}) x 1}{(6 x \code{NumNP} x 1)}-array where the first
    dimension gives
    the (1) current pressure h, (2) discharge rate Q,
    (3) water content theta, (4) x-component vx of Darcian flux,
    (5) z-component vz of Darcian flux, 
    and (6) the solute concentration;
    \code{NumNP} gives the number of nodes; \code{intermediate.result}
    can be used to create movies of water flow at calculation time.}
%  \item{NumNPD}{maximum number of nodes; default is unlimited}
%  \item{NumElD}{maximum number of finite elements; default is unlimited}
%  \item{NumBPD}{maximum number of boundary points; default is unlimited}
%  \item{NMatD}{maximum number of materials; default is unlimited}
  }

\details{
  \code{d} is a profile definition of a type as the return list of
  \code{link{xswms2d}}, or it is 
  a list that should contain the following elements
  \itemize{
    \item{\code{Kat}}{type of flow system; 0:horizontal, 1:axisymmetric;
      2:vertical.}
    
    \item{\code{MaxIt}}{maximum number of iteration during any time step.}
    
    \item{\code{TolTh}}{maximum desired change of water content within one step.}
    
    \item{\code{TolH}}{maximum desired change of pressure head within one step.}
    
    \item{\code{lWat}}{logical. If \code{TRUE} transient water flow else
      steady state.}
    
    \item{\code{FreeD}}{logical. \code{TRUE} if free drainage at the bottom.}
    
    \item{\code{NLay}}{Number of subregions for water balances.}
    
    \item{\code{hTab}}{\code{c(hTab1, hTabN)}, interval of pressure heads
      within which a table of hydraulic properties is generated.}
    
    \item{\code{Par}}{vector of 9 elements or matrix of 9 columns, each row
      representing a different kind of material.}
    
    \item{\code{dt}}{initial time increment.}
    
    \item{\code{dtMinMax}}{\code{c(dtmin, dtmax)}, minimum and maximum permitted
      time increment.}
    
    \item{\code{DMul}}{\code{c(dMul, dMul2)},
      \eqn{\code{dMul}\ge 1}{\code{dMul >= 1}},
      \eqn{\code{dMul2}\le 1}{\code{dMul2 <= 1}};
      if number of required iterations is less than 4 or greater than 6
      then the next time
      step is multiplied by \code{dMul} and \code{dMul2}, respectively.
    }
    
    \item{\code{TPrint}}{vector of increasing print-times.}
    
    \item{\code{NP}}{seepage information: matrix or list.
      matrix of one or two rows and arbitrary columns or
      each row defines a seepage face, the columns
      give the node numbers (or are \code{NA}) or
      a list of vectors;
    }
    \item{\code{DrCorr}}{reduction factor for drainage; see the SWMS2D manual.}
    
    \item{\code{ND}}{vector of up to two elements given the global node number
      of the drain.}
    
    \item{\code{EfDim}}{\code{(length{ND} x 2)-matrix} drain information: first
      column gives the effective diameter; second column gives the
      dimension of the square in the finite element mesh.}
    
    \item{\code{KElDr}}{drainage information: matrix or list;
      matrix of \code{length(ND)} rows and arbitrary colums;
      each row defines the surrounding finite
      elements of the drain; values are either the global finite
      element numbers or \code{NA} or a list of vectors.}
    
    \item{\code{Epsi}}{temporal weighing coefficient: 0: explicit scheme,
      0.5:Crank-Nicholson, 1:fully implicit.}
    
    \item{\code{lUpW}}{logical. Upstream weighing formalation used if \code{TRUE}
      else Galkerin formulation.}
    
    \item{\code{lArtD}}{logical. \code{TRUE} if artificial dispersion is to be
      added to fulfill the stability criterion \code{PeCr}.}
    
    \item{\code{PeCr}}{Stability criterion, see the SWMS2D manual; zero if
      \code{lUpW}.}
    
    \item{\code{ChPar}}{vector of 9 element if \code{Par} is a vector; other
      wise (\code{nrow(Par)} x 9)-matrix; chemical material properties:
      bulk density, ionic diffusion coeffcient in free water,
      longitudinal dispersivity, transverse dispersivity, Freundlich
      isotherm coefficient, first-order rate constant for dissolved
      phase, first-order rate constant for solid phase, zero-order rate
      constant for dissolved phase, zero-order rate constant for solid
      phase.}
    
    \item{\code{KodCB}}{vector of \code{sum(nCodeM[,2]!=0)} elements, defining
      the chemical boundary conditions, see the SWMS2D manual.}
    
    \item{\code{cBound}}{vector of six elements: jth element gives the 
      concentration for boundary nodes \code{i} with chemical boundary
      condition \code{abs(KodCB[i])==j}.}
    
    \item{\code{tPuls}}{time duration of the concentration pulse.}
    
    \item{\code{nCodeM}}{matrix of 11 or 12 columns and arbitrary number of rows;
      \code{nCodeM} gives the nodal information;
      the 12 columns contain: (1) nodal number (might be omitted), (2) code
      giving the boundary condition (see the SWMS2D manual), (3) x-coordinate,
      (4) z-coordinate, (5) initial pressure head, (6) initial concentration,
      (7) prescribed recharge or discharge rate, (8) material number,
      (9) water uptake distribution value, (10) pressure head scaling
      factor, (11) conductivity scaling factor, (12) water content scaling
      factor.
    }
    
    \item{\code{KXR}}{matrix of 8 columns and arbitrary rows; gives the finite
      element information; first 4 colums give the nodes in
      counter-clockwise order (last node is repeated if element is a
      triangle), 5th to 7th column give anisotropy parameters of the
      conductivity tensor (angle, first and second principal component),
      8th column: subregion number.
    }
    
    \item{\code{Width}}{scalar or vector of \code{sum(nCodeM[,2]!=0)} elements;
      width of the boundary associated with the boundary nodes
      (in the order given by nCodeM).}
    
    \item{\code{rLen}}{width of soil surface associated with transpiration}
    
    \item{\code{Node}}{vector of observation nodes for which information is
      collected at each time level.}
    
    \item{\code{SinkF}}{logical; water extraction from the root zone if
      \code{TRUE}.} 
    \item{\code{qGWLF}}{logical; If \code{TRUE} then the discharge-groundwater
      level relationship is used, see the SWMS2D manual.}
    
    \item{\code{GWL0L}}{reference position of groundwater table (usually the
      z-coordinate of the soil surface).}
    
    \item{\code{Aqh}}{parameter in the discharge-groundwater
      level relationship.}
    
    \item{\code{Bqh}}{second parameter in the discharge-groundwater
      level relationship.}
    
    \item{\code{tInit}}{starting time of the simulation.}
    
    \item{\code{hCritS}}{maximum allowed pressure head at soil surface.}
    
    \item{\code{atmosphere}}{matrix of 10 columns and arbitrary rows;
      each row presents the atmospherical conditions at a certain
      instance;
      the columns are (1) the instance at which the time period ends
      (2) precipitation, (3) solute concentration of rainfall water, (4)
      potential evaporation rate, (5) potential transpiration rate,
      (6) absolute value of minimum allowed pressure head at the soil
      surface, (7) drainage flux across the bottom boundary
      (for \code{abs(nCodeM[,2])==3}), (8) ground water level,
      (9) concentration of drainage flux
      (for \code{(abs(nCodeM[,2])==3) & (KodCB<0)})
      (10)  concentration of drainage flux
      (for \code{(abs(nCodeM[,2])==3) & (KodCB>0)}).
    }
    
    \item{\code{root}}{vector of 6 elements or (\code{nrow(nCodeM)} x 6)-matrix.
      \code{c(P0, P2H, P2L, P3, r2H, r2L)} information about
      water uptake by roots, see the SWMS2D manual.}
    
    \item{\code{POptm}}{pressure heads, below which roots start to extract water
      at maximum possible rate;
      if \code{is.vector(root)} then \code{POptm} is a 
      scalar if \code{Par} is a vector, and a vector of
      \code{nrow{Par}} elements otherwise;
      if \code{!is.vector(root)} then \code{POptm} has
      \code{nrow(nCodeM)} elements.
    }
  }
}
\value{
  If an error occured a string containing the error message is returned.
  Otherwise the result is a list of the following elements
  \item{hQThFlC}{array(dim=c(6, \code{NumNP}, n)) where \eqn{n\le
      \code{MPL} + 2}{n <= \code{MPL} + 2} (simulation output);
    first dimension gives (i) water pressure H, (2) discharge rate Q,
    (3) water content theta, (4) x-component vx of Darcian flux,
    (5) z-component vz of Darcian flux, 
    and (6) concentration;
    the third dimension gives the time points given by TPrint except the
    first that gives the initial conditions.
  }
  \item{TlSolObs}{matrix with \eqn{46 + 3 * \code{NObs}} columns
    (simulation output);
    the columns are (1) time point, (2) rAtm, (3) rRoot,
    (4-10) vK[Atm, Root, 3, 1, Seep, 5, 6],
    (11-17) hK[Atm, Root, 3, 1, Seep, 5, 6],
    (18) CumQAP, (19) CumQRP, (20-26) CumQ[Atm, Root, 3, 1, Seep, 5, 6],
    (27) dt, (28) Iter, (29) ItCum, (30) Peclet, (31) Courant,
    (32) CumCh0, (33) CumCh1, (34) ChumR, (35-40) ChemS[1:6],
    (41-46) SMean[1:6], hNew[1:\code{Nobs}] if \code{Nobs} > 0,
    ThNew[1:\code{Nobs}] if \code{Nobs} > 0,
    ConcNew[1:\code{Nobs}] if \code{Nobs} > 0;
    see the SWMS2D manual for details.
  }
  \item{atmOut}{(\code{MaxAl} x 9)-matrix (simulation output);
    (1) AtmTime, (2) CumQAP, (3) CumQRP, (4) CumQA, (5) CumQR, (6) CumQ3,
    (7) hAtm, (8) hRoot, (9) hKode3; see the SWMS2D manual for details.
  }
  \item{balance}{matrix of \code{balance.ncol} columns and at most
    \code{boundary.n} rows (simulation output);
    first 5 columns are (1) time, (2) WatBalT, (3) WatBalT, (4) CncBalT, (5)
    CncBalR; for each of the subjects (i) Area, (ii) Volume, (iii)
    InFlow, (iv) hMean, and (v) ConcVol and (vi) cMean if solute transport,
    \code{NLay +1} subsequent columns are given for the total value and
    the value of each subregion; see the SWMS2D manual for details.
  }
  \item{boundary}{array(dim=c(\code{NumBP}, 11, n)) where \eqn{n\le
      boundary.n}{n <= boundary.n} (simulation output);
    the second dimension gives (1) i, (2) n, (3) x, (4) z, (5) Code,
    (6) Q, (7) v, (8) h, (9) th, (10) Conc, (11) time;
    see the SWMS2D manual for details.
  }
  \item{flux}{\code{d$flux};  different from \code{NULL} if the function
    is called with the output of \code{create.waterflow}
    or with a profile definition, see the return value of
    \command{\link{xswms2d}}.} 
  \item{water.x}{\code{d$water.x};  different from \code{NULL}
    if the function is called with the output of \code{create.waterflow}
    or with a profile definition, see the return value of
    \command{\link{xswms2d}}.} 
  \item{water.y}{\code{d$water.y};  different from \code{NULL} if
    the function is called with the output of \code{create.waterflow}
    or with a profile definition, see the return value of
    \command{\link{xswms2d}}.} 
}
\references{
  \itemize{
    \item
    Schlather, M. and Huwe, B. (2004)
   The use of the language interface of R: two examples for modelling water
   flux and solute transport. \emph{Computers \& Geosciences} \bold{30}, 197-201.

   \item
    Simunek., J., Vogel, T., and van Genuchten, M.Th. (1994)
    \emph{The SWMS2D code for simulating water flow and solute transport
      in two-dimensional variably saturated media, Version 1.21.}
    Research Report No. 132, 197 p.,
    U.S. Salinity Laboratory, USDA, ARS, Riverside, California.
    
  }
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/institute}

  The Fortran code was published as \sQuote{swms\_2d}
  by Jirka Simunek, T. Vogel and Martinus Th. van Genuchten, 
  \url{http://www.ussl.ars.usda.gov/MODELS/MODELS.HTM};
  the file ORTHOFEM.f was written by Edward A. Sudicky  
  and Carl A. Mendoza,
  based on code by Frank W. Letniowshi and contributions by Rene Therrien,
  and modified by Jirka Simunek.
}

\section{Acknowledgement}{
  The author is grateful to Rien van Genuchten, Carl A. Mendoza, Rene Therrien,
  and Edward Sudicky for kindly changing the copyright of the modified
  SWMS2D code, as used in this package, to the GNU copyleft
  licence.

  Note that, however, HYDRUS2D, an advanced development of SWMS2D,
  and its code is proprietary.

  The work has been financially
  supported by the German Federal Ministry of Research
  and Technology (BMFT) grant PT BEO 51-0339476C during 2000-2002.
}

\seealso{\command{\link{xswms2d}},
         \command{\link{read.swms2d.table}},
         \code{\link{SoPhy}}
       }

\examples{
##################################################################
##  SWMS2D, Example 1, stochastically modified material properties
##################################################################
# E = exp(mu + 1/2 * sd^2)
# Var = exp(2*mu + sd^2)*(exp(sd^2) - 1)
%# print(c("log.sd=",sqrt(exp(2*mu + sd^2)*(exp(sd^2) - 1))))

path <- paste(system.file(package='SoPhy'), 'swms2d', sep="/")
x <- read.swms2d.table(path)
x$TPrint <- 20 * 60
sd <- 0.3 ## then sd of log.gauss is 0.306878
mu <- log(x$nCodeM[, 10]) - 0.5 * sd^2
n <- if (interactive()) 100 else 10 # better 10000
front <- integer(n)

for (j in 1:n) {
  cat(j,"\n")
  x$nCodeM[, 10] <- rlnorm(length(x$nCodeM$z), m = mu, sd = sd)
  x$nCodeM[, 11] <- x$nCodeM[, 10]^(-2)
  z <- swms2d(x)$hQ[3, , ]
  front[j] <- (min(which(z[-1:-2, 2] / z[-1:-2, 1] < 1.005)) \%/\% 2) + 1
}

idx <- 24 
cat('The probability that the front advances at most', 
    max(x$nCodeM$z) - x$nCodeM$z[2 * idx],'cm is about',
    format(mean(front <= idx), dig=2),'\n')


par(cex=1, mar=c(4.2,4.2,1.2,0.4))
distance <- max(x$nCodeM$z) - x$nCodeM$z[2 * front]
h <- sort(unique(distance))
dh <- min(diff(h) / 2)
hist(distance, freq=FALSE,
     breaks=seq(h[1]-dh, h[length(h)]+dh, 2 * dh),
     main="Histogramm for the depth of the water front",
     xlab="depth", cex.axis=1.5, cex.main=1.2, cex.lab=1.5)

}
\keyword{spatial}


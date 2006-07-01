\name{read.swms2d.table}
\alias{read.swms2d.table}
\title{Standard SWMS2D input files}
\description{
 This function reads in the standard SWMS2D input files
}
\usage{
read.swms2d.table(path, selct.in = "SELECTOR.IN", 
                  grid.in = "GRID.IN", atm.in = "ATMOSPH.IN")
}

\arguments{
  \item{path}{directory of the input files}
  \item{selct.in}{file name for general input data, see the manuscript of SWMS2D}
  \item{grid.in}{file name for information about the finite elements}
  \item{atm.in}{file name for atmospheric data}
}
\value{
  the function returns a list of the following components;
  see the SWMS2D manual for details
  \item{Units}{\code{c(LUnits, TUnits, MUnits, BUnits)},
    units of length, time, mass.}
  \item{Kat}{type of flow system; 0:horizontal, 1:axisymmetric; 2:vertical.}
  \item{MaxIt}{maximum number of iteration during any time step.}
  \item{TolTh}{maximum desired change of water content.}
  \item{TolH}{maximum desired change of pressure head.}
  \item{lWat}{logical. If \code{TRUE} transient water flow else steady state.}
  \item{FreeD}{logicl. \code{TRUE} if free drainage at the bottom.}
  \item{NLay}{Number of subregions for water balances.}
  \item{hTab}{\code{c(hTab1, hTabN)}, interval of pressure heads
    within which a table of hydraulic properties is generated.}
  \item{Par}{vector of 9 elements or matrix of 9 columns, each row
    representing a different kind of material.}
  \item{dt}{initial time increment.}
  \item{dtMinMax}{\code{c(dtmin, dtmax)}, minimum and maximum permitted
    time increment.}
  \item{DMul}{\code{c(dMul, dMul2)},
    \eqn{\code{dMul}\ge 1}{\code{dMul >= 1}},
    \eqn{\code{dMul2}\le 1}{\code{dMul2 <= 1}};
    if the number of required iterations is less than 4 or greater than 6
    then the next time
    step is multiplied by \code{dMul} and \code{dMul2}, respectively.
  }
  \item{TPrint}{vector of increasing print-times.}
  \item{NP}{matrix of one or two rows and arbitrary columns;
    seepage information: each row defines a seepage face, the columns
    give the node numbers (or are \code{NA}).}
  \item{DrCorr}{reduction factor for drainage; see the SWMS2d manual}
  \item{ND}{vector of up to two elements given the global node number
    of the drain.}
  \item{EfDim}{\code{(length{ND} x 2)-matrix} drain information: first
    column gives the effective diameter; second column gives the
    dimension of the square in the finite element mesh.}
  \item{KElDr}{matrix of \code{length(ND)} rows and arbitrary colums;
    drainage information: each row defines the surrounding finite elements
    of the drain; values are either the global finite element numbers
    or \code{NA}.}
  \item{Epsi}{temporal weighing coefficient: 0: explicit scheme,
    0.5:Crank-Nicholson, 1:fully implicit.}
  \item{lUpW}{logical. Upstream weighing formalation used if \code{TRUE}
  else Galkerin formulation.}
  \item{lArtD}{logical. \code{TRUE} if artificial dispersion is to be
    added to fulfill the stability criterion \code{PeCr}.}
  \item{PeCr}{Stability criterion, see the SWMS2d manual; zero if
    \code{lUpW=TRUE}.} 
  \item{ChPar}{vector of 9 element if \code{Par} is a vector; other
    wise (\code{nrow(Par)} x 9)-matrix; chemical material properties:
      bulk density, ionic diffusion coeffcient in free water,
      longitudinal dispersivity, transverse dispersivity, Freundlich
      isotherm coefficient, first-order rate constant for dissolved
      phase, first-order rate constant for solid phase, zero-order rate
      constant for dissolved phase, zero-order rate constant for solid
      phase.} 
  \item{KodCB}{vector of \code{sum(nCodeM[,2]!=0)} elements, defining
    the chemical boundary conditions, see the SWMS2d manual.}
  \item{cBound}{vector of six elements: jth element gives the 
    concentration for boundary nodes \code{i} with chemical boundary
    condition \code{abs(KodCB[i])==j}.}
  \item{tPuls}{time duration of the concentration pulse.}
  \item{nCodeM}{matrix of 12 columns and arbitrary number of rows;
    \code{nCodeM} gives the nodal information;
    the 12 columns contain: (1) nodal number, (2) code
    giving the boundary condition (see the SWMS2D manual), (3) x-coordinate,
    (4) z-coordinate, (5) initial pressure head, (6) initial concentration,
    (7) prescribed recharge or discharge rate, (8) material number,
    (9) water uptake distribution value, (10) pressure head scaling
    factor, (11) conductivity scaling factor, (12) water content scaling
    factor.
  }
  \item{KXR}{matrix of 8 columns and arbitrary rows; gives the finite
    element information; (1)-(4) give the nodes in
    counter-clockwise order (last node is repeated if element is a
    triangle); (5)-(7) give the anisotropy parameters of the
    conductivity tensor (angle, first and second principal component),
    (8) subregion number.
  }
  \item{Width}{vector of \code{sum(nCodeM[,2]!=0)} elements;
    width of the boundary associated with the boundary nodes
    (in the order given by nCodeM).}
  \item{rLen}{width of soil surface associated with transpiration.}
  \item{Node}{vector of observation nodes for which information is
    collected at each time level.}
  \item{SinkF}{logical; water extraction from the root zone if \code{TRUE}.}
  \item{qGWLF}{logical; If \code{TRUE} then the discharge-groundwater
    level relationship is used, see the SWMS2d manual.}
  \item{GWL0L}{reference position of groundwater table (usually the
    z-coordinate of the soil surface).}
  \item{Aqh}{parameter in the discharge-groundwater
    level relationship.}
  \item{Bqh}{second parameter in the discharge-groundwater
    level relationship.}
  \item{tInit}{starting time of the simulation.}
  \item{hCritS}{maximum allowed pressure head at soil surface.}
  \item{atmosphere}{matrix of 10 columns and arbitrary rows;
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
    (for \code{(abs(nCodeM[,2])==3) & (KodCB>0)}).}
  \item{root}{\code{c(P0, P2H, P2L, P3, r2H, r2L)} information about
    water uptake by roots, see the SWMS2D manual.}
  \item{POptm}{scalar if \code{Par} is a vector, a vector of
    \code{nrow{Par}} elements otherwise;
    pressure heads, below which roots start to extract water
    at maximum possible rate.}
  \item{path}{the input parameter.}
  \item{selct.in}{the input parameter.}
  \item{grid.in}{the input parameter.}
  \item{atm.in}{the input parameter.}
}
\references{
  \itemize{
    \item
    Schlather, M. and Huwe, B. (2004)
    The use of the language interface of R: two examples for modelling water
    flux and solute transport.
    \emph{Computers \& Geosciences} \bold{30}, 197-201.
    
    \item
    Simunek, J., Vogel, T., and van Genuchten, M.Th. (1994)
    \emph{The SWMS2D code for simulating water flow and solute transport
      in two-dimensional variably saturated media, Version 1.21.}
    Research Report No. 132, 197 p.,
    U.S. Salinity Laboratory, USDA, ARS, Riverside, California.
  }
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/institute}
}


\seealso{\command{\link{swms2d}}%, \code{\link{SoPhy}}
  }

\examples{
##################################################################
##  SWMS2D, Example 1
##################################################################
%##col <- rev(topo.colors(100))
%##col <- grey(seq(1,0.15,-0.01))

%library(lib="~/TMP", SoPhy); source("../R/swms2d.R");

path <- paste(system.file(package='SoPhy'), 'swms2d', sep="/")
x <- read.swms2d.table(path)
x$TPrint <- seq(10, 5400, 10)
par(cex=1, mar=c(4.2,4,0.2,0.2))

z <- swms2d(x)$hQ[3, , ]
i <- ((length(x$nCodeM$z) / 2):1) * 2 - 1
image(x=c(0, x$TPrint) / 60, y=x$nCodeM$z[i] - max(x$nCodeM$z),
      z=t(z[i, ]), xlab='time [min]', ylab='z [cm]',
      col=grey(seq(1, 0.15, -0.01)), cex.lab=1.5, cex.axis=1.5)

my.legend(0, -max(x$nC$z[i]), zlim=range(z), y.i=0.02,
          col=grey(seq(1,0.15,-0.01)), cex=1.5)
}

\keyword{spatial}

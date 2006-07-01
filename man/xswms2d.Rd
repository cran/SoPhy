% R CMD Rdconv -d -t txt SoPhy/man/xswms2d.Rd

% R CMD Rdconv -t txt xswms2d.Rd
\name{xswms2d}
\alias{sophy}
\alias{xswms2d}
\title{Modelling of water flux and solute transport}
\description{
  Interactive surface that is used to define a soil profile and
  to give the physical and chemical properties for each horizon.
  It simulates the water flux and solute transport using swms2d by
  J. Simunek,  T. Vogel, and M.Th. van Genuchten.
}
\usage{
sophy(...)

xswms2d(h, xlim, ylim, step, picture=NULL,
  %## material and water flux parameters
  water=list(
    %## printing
    TPrint=500, red=4, print=7, mesh=TRUE, lim.rf = 1,
    lim.swms2d = 1, 
    %## simulation control
    TolTh=0.0001, TolH=0.01, lWat=TRUE, dt=1, dtMinMax = c(0.01, 60),
    max.iteration=1000,  max.iter.prec=50000,
    iter.print=100, breakpoint=100,
    %## boundary conditions
    top.bound=2, top.value=0, bottom.bound=1, bottom.value=0
  ),
  materials=list(
     thr=.02, ths=0.35, tha=0.02, thm=0.35, Alfa=0.041, n=1.964,
     Ks=0.000722, Kk=0.000695, thk=0.2875,
     first=1, second=1, angle=0, 
     Hslope=0, Hseg=-100, POptm=-25, sharpness=1, 
     Bulk.d=1500, Diffus=0, longiDisper=1, transvDisper=0.5,
     Adsorp=0.0004, SinkL1=-0.01, SinkS1=-0.01, SinkL0=0, SinkS0=0
  ),
  Hinit = function(Hseg, depth) Hseg,
  model=list(model="exp", param=c(0, 0.25, 0, diff(ylim) * 0.1)),
  anisotropy=NULL,
  miller.link=function(rf) exp(sign(rf) * sqrt(abs(rf))),
  millerH = function(lambda) lambda, 
  millerK = function(lambda) lambda^-2, 
  millerT = function(lambda) 1,

  %## chemical data
  chemical=list(
     lChem=FALSE, Epsi=2, lUpW=1, lArtD=FALSE, PeCr=10, tPulse=500,
     top.bound=2, top.value=1, bottom.bound=1, bottom.value=0,
     root.uptake=0, intern.source=0
  ),
  
  %## atmospherical data
  atmosphere=list(
     AtmInf=TRUE, tInit=0, Aqh=-.1687, Bqh=-.02674,
     hCritS=1.e30, GWL0L=ylim[2], rLen=diff(xlim)
  ),
  atm.periods=1, 
  atm.data=c(
     end=100000000, prec=0, cPrec=0, rSoil=0, rRoot=0, 
     hCritA= 1000000, rGWL=0, GWL=0, crt=0, cht=0
  ),

  %## stones
  stone=list(
     %## basics
     value=NA, lambda=0,
     %## location
     no.upper=FALSE, no.lower=TRUE, no.overlap=FALSE,
     %## size 
     main.distr=rnorm,  main.mean=diff(ylim)/20, main.s=diff(ylim)/400,
     sec.distr=rnorm, sec.mean=diff(ylim)/50, sec.s=diff(ylim)/400,
     phi.distr=rnorm, phi.mean=1, phi.s=0
  ),

  %## roots
  plant.types=1,
  Kf=function(plant, distances, depth, rf, 
              param=list(field.value=1, field.factor=0))
     rep(param$field.value, length(rf)) + param$field.factor * rf,
  beta=function(plant, distances, depth, rf, 
	      param=list(beta.value=1, beta.factor=0))
    rep(param$beta.value, length(rf)) + param$beta.factor * rf,
  root=list(
     plants.lambda=0,
     plants.mindist=diff(xlim) / 20,
     mean=sqrt(diff(ylim) * diff(xlim)),
     sd=sqrt(diff(ylim) * diff(xlim)) / 5,
     %## root growth
     knot.probab=0.1,
     knot.mindist = 5 * step,
     shoots.3=0.4, shoots.4=0,
     stop.probab=function(knot, dist, m, s) 1 - exp(-m + s * dist),
     stop.m=0, stop.s=0,
     rf.link=function(v, m, s){
               v <- s * v
               m * (exp(sign(v) * sqrt(abs(v))))^(-2)
             },
     rf.m=0, rf.s=0,
     %## advanced root growth   
     no.own.root=TRUE, age.bonus=1.2, depth.bonus=5.2,
     side.bonus=5.1, diagonal=TRUE, dir.ch=3.71, dir.ch.s=0.5,
     %## water uptake / change of conductivity
     rf.Kf=FALSE, P0=-10, P2H=-200, P2L=-800, P3=-8000, r2H=0.5,
     r2L=0.1, root.condition=3, root.uptake=-150
     ),
  %##
  root.zone = NULL,
  %## images
  col.rf = NULL, col.simu = NULL,
  %## control parameters for swms2d
  MaxIt=20, hTab=c(0.001,200), DMul=c(1.1, 0.33),
  %## menue, parameter input
  col.rect="red", col.bg="yellow", col.sep="gray", col.left="red",
  col.mid="darkgreen", col.right="darkgreen", col.line="red",
  col.txt="black", col.submenue="darkgreen", col.subord="steelblue",
  col.forbid="gray88", col.bg.forbid="lightyellow", col.flash="blue",
  ## drawing configure CFLAGS="-O2 -Wall -pedantic"


  col.hor=c("#000000", "#996600", "#660000",  "#CC9933", "#666600",
            "#CCCC99", "#CCCCCC", "#990000", "#FFCC00", "#FFFFCC"),
  col.draw="green", col.mess="red",
  %## finite element mesh
  col.mesh="blue", col.mesh.pt="green",
  col.exception = c("brown", "lightgray"),
  cex.eval = 0.8,
  areas=is.null(h$picture), PrintLevel=RFparameters()$Print,
  new=TRUE, X11.width=9.5, X11.height=X11.width * 0.789,
  frequent.reset = TRUE, update=FALSE, waterflow=FALSE, 
  zlim = NULL, Zlim = NULL,
  %## postscript output
  print.par=list(ps="sophy", height=3, titl=TRUE, legend = TRUE)   
)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
 
  \item{h}{a list of the output format, see the Details.
    If \code{h} is given the parameters \code{xlim}, ..., \code{col.simu}
    are ignored.}

  \item{xlim}{vector of two elements.  The horizontal extension of the
    profile, see also 
    \code{picture}.  If \code{diff(xlim)} is not an integer multiple of
    step then \code{xlim[2]} is decreased suitably.}

  \item{ylim}{vector of two elements.  The vertical extension of the
    profile, see also 
    \code{picture}.  If \code{diff(ylim)} is not an integer multiple of
    step then \code{ylim[2]} is decreased suitably.}

  \item{step}{numeric.
    the grid spacing for both directions.  The finite element
    method is essentially based on a quadratic mesh.}

  \item{picture}{array or file name for a digitised profile in tif or jpeg
    format; if \code{picture}
    is an array it must have three dimensions, where the third
    dimension has 3 or 4 components (RGB or RGBA format).

    If \code{picture} is given then either \code{xlim} or
    \code{ylim} might be missing and is calculated from the other lim
    parameter and the dimension of the picture.}
  %
  \item{water}{
    \{swms2d, water\} list of the following components
    \describe{
      \item{\code{TPrint} \{TPrint\}}{scale; 
	modelling time after which the result is printed.
	
	Note that for directly calls of \command{\link{create.waterflow}}
	and subsequent call of \command{\link{swms2d}}, \code{TPrint}
	might also be a vector, see \command{\link{read.swms2d.table}}
	for futher information.
      }
      
      \item{\code{red} \{size reduction\}}{
	 positive integer.
	to increase the simulation spead the
	size of the finite element grid can coarsened}
      
      \item{\code{print} \{printed variable\}}{
	integer 1,...,6. 
	One of the following
	spatial data sets can be plotted:
	the water potential \eqn{H},
	discharge/recharge rates \eqn{Q} for internal sink/sources,
	the water contents \eqn{theta},
	the x-components of the Darcian flux vector vx,
	the z-components of the Darcian flux vector vz,
	or the soulte concentration Conc.  These variables
	are numbered
	consecutively form 1 to 6.
      }
      
      \item{\code{mesh} \{show FE mesh\}}{
	 logical.  If \code{TRUE} the currently
	used finite element mesh in a swms2d simulation
	is shown in the upper left drawing area.}
       
      \item{\code{lim.rf} \{upper quantile for random field plot\}}{
	numeric in \eqn{[0,1]}.  The value of \code{zlim[2]} in
	\code{image} equals the \code{lim.rf} quantile of
	the random field values when the transformed random field is plotted.}

      \item{\code{lim.swms2d} \{quantile for swms2s plot\}}{
	numeric in \eqn{[0,1]}. See \command{\link{plotWater}} for details.}
      
      \item{\code{TolTh} \{TolTh\}}{
	maximum desired change of water content
	within one step of the FE method.}
      
      \item{\code{TolH} \{TolH\}}{
	 maximum desired change of pressure head within
	one step of the FE method.}
      
      \item{\code{lWat} \{LWat\}}{
	 logical. If \code{TRUE} transient water flow else
	steady state.}
      
      \item{\code{dt} \{dt - initial time step\}}{
	 initial time increment of the FE method.}

      \item{\code{dtMinMax} \{minimal time step\} and \{maximal time step\}}{
	\code{c(dtmin, dtmax)}, minimum and maximum permitted
	time increment of the FE method.}
      
      \item{\code{max.iteration} \{max. iterations\}}{
	 maximum number of
	incremental time steps for coarsed finite element meshes.
	% coarsening by factor one is also possible
	Once \code{max.iteration} is reached, the current time and
	the aimed time is shown.  The user is asked whether it should be
	continued.}
      
      \item{\code{max.iter.prec} \{max. iter. (precise)\}}{
	same as
	\code{max.iteration} except that this value is used if the main
	menue bottom \dQuote{precise waterflow} is pressed.}
      
      \item{\code{iter.print} \{swms message, loop period\}}{number of
	incremental time steps after which a short message of status is
	given.}
      
      \item{\code{breakpoint} \{swms, image, loop period\}}{
	 number of incremental
	time steps after which the current value of the spatial variable
	is shown in a two dimensional plot.  If
	\code{breakpoint} is chosen appropriately the user gets the
	impression of a movie.  Note that pictures are presented
	after a certain amount of iteration loops,
	so the length of the simulated time intervals between the
	pictures differ.}
      
      \item{\code{top.bound} \{top\}}{
	 surface boundary condition: 1 (impermeable),
	2 (Dirichlet), 3 (Neumann), 4 (variable H), 5 (variable Q), or
	6 (atmospheric).  See the SWMS2D manual for details.}
      
      \item{\code{top.value} \{top value\}}{
	 the value of the \code{top.bound} condition, if
	appropriate.}
      
      \item{\code{bottom.bound} \{bottom\}}{
	 boundary condition at the bottom:
	1 (free drainage), 2 (deep drainage), 3 (impermeable), 4 (Dirichlet), or
	5 (Neumann).  See the SWMS2D manual for details.}
      
      \item{\code{bottom.value} \{bottom value\}}{
	 the value of the \code{bottom.bound} condition, if
	appropriate.}
    }
  }
  %
  \item{materials}{\{material (phys) / material (chem)\};
    \code{thr},...,\code{sharpness} are found in the menue for the
    physical 
    properties, \code{Bulk.d},...,\code{SinkS0} in the menue for
    the chemical
    ones.  List has the following components
    \describe{
      \item{\code{thr} \{\eqn{\theta_r}{theta_r}\}}{
	residual water contents
	\eqn{\theta_r}{theta_r}}
      
      \item{\code{ths} \{\eqn{\theta_s}{theta_s}\}}{
	staturated water contents
	\eqn{\theta_s}{theta_s}}
      
      \item{\code{tha} \{\eqn{\theta_a}{theta_a}\}}{
	fictious parameter
	\eqn{\theta_a\le\theta_r}{theta_a<=theta_r}, introduced
	into the van Genuchten model for more flexibility (to allow for
	a non-zero air-entry value), see page 10 of the SWMS2d manual
	for details.}
      
      \item{\code{thm} \{\eqn{\theta_m}{theta_m}\}}{
	\eqn{\theta_m\ge\theta_s}{theta_m>=theta_s}, introduced 
	into the van Genuchten model for more flexibility (to allow for
	a non-zero air-entry value)}
      
      \item{\code{Alfa} \{\eqn{\alpha}{alpha}\}}{
	 factor \eqn{\alpha}{alpha} in
	 the van Genuchten model}
      
       \item{\code{n} \{\eqn{n}{n}\}}{
	 exponent \eqn{n} in the van Genuchten model}
       
       \item{\code{Ks} \{\eqn{K_s}{K_s}\}}{
	 saturated hydraulic conductivity \eqn{K_s}{K_s}}
       
       \item{\code{Kk} \{\eqn{K_k}\}}{
	 non-saturated hydraulic conductivity
	 \eqn{K_k} measured for some \eqn{\theta_k}{theta_k}}
       
       \item{\code{thk} \{\eqn{\theta_k}{theta_k}\}}{see \eqn{K_k}{K_k}}

       \item{\code{first} \{1st principal comp\}}{
	 anisotropy parameter for the
	 conductivity tensor}

%}}}\section{XXX}{{{
      
      \item{\code{second} \{2nd principal comp\}}{
	 anisotropy parameter for the
	conductivity tensor}
      
      \item{\code{angle} \{angle (degrees)\}}{
	  anisotropy parameter for the
	conductivity tensor}
      
      \item{\code{Hslope} \{intial H, slope\}}{
	 The initial matrix potential is
	calculated from the
	\code{miller.link}ed and then \code{millerH}-transformed
	Gaussian random field by linear transformation. \code{Hslope}
	and and a function of \code{Hseg}, see \code{Hinit},
	give the slope and the segment of that linear
	transformation.}
      
       \item{\code{Hseg} \{initial H, segment\}}{
	 see \code{Hslope} and \code{Hinit}}
      
      \item{\code{POptm} \{POptm (root)\}}{
	 highest (negative) matrix potential
	(so close to saturation) for which the water uptake is still
	optimal; see \code{root$P2H} for further information}
      
      \item{\code{sharpness} \{sharpness\}}{
	The values for a random field (of the hydrolic conductivity)
	of a simulated horizon are obtained by a
	linear combination of a conditional simulation (towards random
	fields of the preceding horizons) and an
	independent simulation.
	\code{sharpness} is the factor for the independent simuation
	and
	\eqn{\sqrt{1-\code{sharpness}^2}}{(1-\code{sharpness}^2)^{1/2}}
	the factor for the conditional simulation.
	The mean of the random field to be simulated
	has been substracted beforehand from the boundary
	conditions.
	
	\code{sharpness} is not used for the initial horizon, thus not
	given as a menue point in this case.}

      %
      \item{\code{Bulk.d} \{bulk density\}}{
	 bulk density}
      
      \item{\code{Diffus} \{diffusion\}}{
	 ionic diffusion coeffcient in free water}
    
      \item{\code{longiDisper} \{longitud. Dispers.\}}{
	 longitudinal dispersivity}
      
      \item{\code{transvDisper} \{transvers. Dispers.\}}{
	  transverse dispersivity}
      
      \item{\code{Adsorp} \{Freundlich isoth coeff.\}}{
	Freundlich isotherm coefficient}
      
      \item{\code{SinkL1} \{1st-order const./dissolved\}}{
	 first-order rate
	constant \eqn{\mu_w}{mu_w} for dissolved phase, see page 15 of
	the SWMS2d manuscript}
      
      \item{\code{SinkS1} \{1st-order const./solid\}}{
	 first-order rate constant
	\eqn{\mu_s}{mu_s} for solid phase}
      
      \item{\code{SinkL0} \{0-order const./dissolved\}}{
	 zero-order rate
	constant \eqn{\gamma_w}{gamma_w} for dissolved phase}
      
      \item{\code{SinkS0} \{0-order const./solid\}}{
	 zero-order rate constant
	\eqn{\gamma_s}{gamma_s} for solid phase}
    }
  }
  %
  \item{Hinit}{function of two variables; the first varible is
    a vector of \code{Hseg}, the second the depths \eqn{d}.
    the initial \eqn{h} values are calculated as
    \deqn{\code{Hinit}(\code{Hseg}, d) + \code{Hslope} * \code{millerH}(Z)}  
    where \eqn{Z} is the simulated random field.
    The initial \eqn{h} value may afterwards be modified for certain
    segments due to the given boundary conditions and to a given given
    Dirichlet condition for the roots.
  }
  \item{model}{
    \{structure\} the covariance model for the Gaussian random
    fields used for the definition of the Miller similar media.
    See \command{\link[RandomFields]{CovarianceFct}} for details on the
    definiton of the covariance model.
    
    If the mean equals \code{NA} and if it is the last horizon,
    the area is interpreted as air.}
  
  \item{anisotropy}{
    logical or \code{NULL}. If logical it overrides the
    anisotropy condition given by model.  If \code{model} was
    anisotropic and \code{anisotropic=FALSE}, the first anisotropic
    scale parameter is used as scale for the isotropic model.
    If \code{model} is isotropic, only the representation is changed,
    and the user has the opportunity to change the values to an
    truely anisotropic model.}
  
  \item{miller.link}{
    function that transforms the Gaussian
    random field to a 
    field of non-negative values.  The argument and the value are
    matrices (of the same size).  If \code{NULL} the \code{miller.link}
    is the identity.
    }

  \item{millerH}{function that transforms the \code{miller.link}ed
    field into a field of water potential \eqn{H}.
    This field is used for areas of constant potential, and for the
    potential at starting time.  For the latter the \code{millerH}
    transformed field is further transformed linearly by
    a function of
    \code{$matrials$Hslope}, see \code{Hinit}, and \code{$matrials$Hseg}.
    The argument and the value are
    matrices (of the same size).}
  
  \item{millerK}{function that transforms the \code{miller.link}ed
    field into a field of scaling factors that are
    associated with the saturated hydraulic
    conductivity.  The argument and the value are
    matrices (of the same size).}

  \item{millerT}{function that transforms the \code{miller.link}ed
    field into a field of scaling factors that are
    associated with the water content. 
    The argument and the value are
    matrices (of the same size).}
  
  %
  \item{chemical}{
    \{swms2d (chem)\} list of the following components
    \describe{
      \item{\code{lChem} \{Solute transport\}}{ 
	 logical that indicates whether solute
	transport should be modelled}
      
      \item{\code{Epsi} \{scheme\}}{integer 1,2,3. 
	 Temporal weighing coefficient:
	1 (explicit: weight=0), 2 (Crank-Nicholson: weight=0.5) or
	3 (implicit: weight=1)}
      
      \item{\code{lUpW} \{formulation\}}{
	 1 (upstream weighing formulation)
	or 2 (Galerkin formulation)}
      
      \item{\code{lArtD} \{added artific. disp\}}{
	 logical.  If \code{TRUE}
	artificial dispersion is added to satisfy the stability
	criterion \code{PeCr}}
      
      \item{\code{PeCr} \{PeCr\}}{
	 Stability criterion (Peclet number * Courant
	number); usually 2, but can be increased to 5 or 10;
	see page 44 of the SWMS2D manuscript.}
      
      \item{\code{tPulse} \{pulse duration\}}{
	 time duration of the
	concentration pulse}

      %
      \item{\code{top.bound} \{top\}}{
	 boundary condition for soil surface:
	1 (impermeable), 2 (Dirichlet), 3 (Neumann), 4 (variable
	Concentration), 5 (variable Q), 6 (atmospheric)}
      
      \item{\code{top.value} \{top value\}}{
	 value of the boundary condition if
	appropriate}
      
      \item{\code{bottom.bound} \{bottom\}}{
	 boundary condition for the bottom:
	1 (free drainage), 2 (deep drainage), 3 (impermeable)}
      
      \item{\code{bottom.value} \{bottom value\}}{
	 value of the boundary condition
	if appropriate}
      
      \item{\code{root.uptake} \{root uptake\}}{
	 - not implemented yet}
      
      \item{\code{intern.source} \{internal source\}}{
	  - not implemented yet}
      %
    }
  }
  
  %
  \item{atmosphere}{\{atmosphere, control\}
    list of the following components 
    \describe{
      \item{\code{AtmInf} \{use atmospheric data\}}{
	 logical.  Determines whether or not atmospherical data should
	 be used.}
      
      \item{\code{tInit} \{start time of simu\}}{
	 starting time of the
	simulation.  This information is necessary since the
	\code{atm.data} (see below)
	give the end of an atmospheric period only.
	So for the first period the starting time is missing.}
      
      \item{\code{Aqh} \{\eqn{A_{qh}}\}}{
	 only used if \code{water$bottom.bound}
	equals 2 (deep drainage).  Then the modulus of the discharge rate
	equals 
	\code{step}*\eqn{A_{qh}}\eqn{\exp(B_{qh} |h -
	  \mbox{GWL0L}|)}{exp(B_qh |h - GWL0L|)} where \eqn{h} is the
	local 
	pressure head.} 
      
      \item{\code{Bqh} \{\eqn{B_{qh}}\}}{
	 see \code{Aqh}}
      
      \item{\code{hCritS} \{max pressure at surface\}}{
	 maximum allowed pressure
	head at the soil surface}
      
      \item{\code{GWL0L} \{ref. pos of groundwater\}}{
	 Reference position of the
	groundwater table (usually the z-coordinate of the soil surface);
	only used if
	\code{water$bottom.bound} equals 2 (deep drainage)}

      \item{\code{rLen} \{width of soil (transpiration)\}}{
	width of soil surface associated with transpiration.  Only used
	in case of atmospheric root uptake.
	}
    }
  }
  \item{atm.periods}{
    (maximum) number of atmospherical periods. Only
    used if \code{atm.data} is a vector.  Otherwise \code{atm.periods}
    is set to the number of the rows of \code{atm.data}.}
  \item{atm.data}{\{atmosphere, data\}
    vector of 10 components or matrix of 10 columns.
    The 10 components are
    \describe{
      \item{\code{end} \{end point\}}{
	 end point of the atmospherical
	period; the value of \code{end} for the last atmospherical
	period must be greater than the value of \code{TPrint}}
      
      \item{\code{prec} \{precipitation\}}{
	precipitation %[L/T]
      }
      
      \item{\code{cPrec} \{solute, precipitation\}}{
	 solute concentration
	 of rainfall water %[\eqn{M/L^3}]
       }
      
      \item{\code{rSoil} \{potential evaporation\}}{
	potential evaporation rate %[L/T]
      }
      
      \item{\code{rRoot} \{potential transpiration\}}{
	potential transpiration rate %[L/T]
      }
      
      \item{\code{hCritA} \{abs. min. pressure at surface\}}{
	absolute
	value of the minimum allowed pressure head at the soil surface}
      
      \item{\code{rGWL} \{drainage flux (drain or var. Q)\}}{
	 drainage
	flux or other time-dependent prescribed flux boundary
	condition. For nodes with atmospheric flux condition.}
      
      \item{\code{GWL} \{ground water level\}}{
	 groundwater level or other
	time-dependent prescribed head boundary condition. For nodes
	with atmospheric head condition.  Pressure is then
	\code{GWL}+\code{GWL0L}.}
      
      \item{\code{crt} \{conc. flux (drainage)\}}{
	 time-dependent concentration flux}
      
      \item{\code{cht} \{conc. \sQuote{pressure} (drain/var. H)\}}{
	 time-dependent
	concentration of the first-type boundary condition.  See the
	SWMS2D manual for details.}
    }
  }
  %

  \item{stone}{\{stones\} list of the following components
    \describe{
      \item{\code{value} \{value\}}{
	 numeric or \code{NA}. \code{NA} is used to
	model stones, i.e. water cannot penetrate.
	If not \code{NA},
	the (Gaussian) random field gets, at the place of the
	\sQuote{stone}, 
	the value \code{value}.  This might be of interest if small
	lenses are modelled by means of \sQuote{stones}.}
      
      \item{\code{lambda} \{intensity\}}{
	 intensity of the stones}
      
      \item{\code{no.upper} \{!overlap into upper horizons\}}{
	 if \code{TRUE}
	overlap into any upper horizon is forbidden}
      
      \item{\code{no.lower} \{!overlap into lower horizons\}}{
	 if \code{TRUE}
	overlap into any lower horizon is forbidden}
      
      \item{\code{no.overlap} \{!overlap of stones\}}{
	 if \code{TRUE} the
	 overlap of stones is forbidden. This point
	 is clear from a practical
	 point of view.  However the currently implemented algorithm
	 (SSI/RSA)
	 for non-overlapping stones is slow and it is furthermore not
	 guarantied that it will not fail (due to \sQuote{bad} random
	 numbers). 
	If overlapping is allowed a fast and simple algorithm is
	used (Boolean model).
	See for both, SSI and Boolean model, the references in Stoyan et
	al., for example.}
      
      \item{\code{main.distr}}{
	distribution for the length of the main axis of
	the ellipse}
      
      \item{\code{main.mean} \{main axis, mean\}}{
	 the parameter for the mean of
	the distribution}
      
      \item{\code{main.s} \{main axis, sd\}}{
	 the parameter for the standard
	deviation of the distribution}
      
      \item{\code{sec.distr}}{
	distribution for the length of the second axis of
	the ellipse}
      
      \item{\code{sec.mean} \{second axis, mean\}}{
	 the parameter for the mean of
	the distribution}
      
      \item{\code{sec.s} \{second axis, sd\}}{
	 the parameter for the standard
	deviation of the distribution}
      
      \item{\code{phi.distr}}{
	distribution for the angle between the horizontal
	line and the main axis of the ellipse}
      
      \item{\code{phi.mean} \{second axis, mean\}}{
	 the parameter for the mean of
	the distribution}
      
      \item{\code{phi.s} \{second axis, sd\}}{
	 the parameter for the standard
	deviation of the distribution}
    }
  }
  
  \item{plant.types}{positive integer.  Number of different types of
    root systemes that will be generated}

  \item{Kf}{\code{function(plant, distances, depths, rf, 
      param=list(field.value=1, field.factor=0))}.
    According to the parameters the (Gaussian) random field get
    a new value if there is a root pixel.
    Here, \code{plant} is the type
    of plant (scalar integer value 1..\code{plant.types}),
    \code{distances} a vector of distances from the
    considered locations
    to the surface along the root, \code{depth} is a vector of depths
    of the considered locations, and \code{rf} is the value of the
    Gaussian 
    random field (including the stones).
    \code{param} is a list of named elements (possibly empty)
    whose values can be
    modified interactively; all parameters are real valued and the
    usual magnitude of the range is 1; the names may not match any
    name in the list for \code{root}, since the parameters are included
    in the \code{root} list, which is passed to \code{Kf} as a whole.
    \cr
    
    The function is only
    used if \code{root$rf.Kf} is \code{TRUE} for the specific
    plant type.  This function is still in an experimental stage.
  }


  \item{beta}{\code{function(plant, distances, depths, rf, 
                       param= list(beta.value=1, beta.factor=0))}.
    \code{beta} gives the raw potential root 
    water uptake according to a single \code{plant} type, the
    \code{distances} 
    to the beginning of the root, the \code{depths}, and \code{rf} (see
    also \code{Kf}).  The raw potential
    values are subsequently normed such that
    they sum up to one.  This gives
    the potential root water uptake.  \code{param} is a list of named
    elements (possibly empty) whose values can be can be 
    modified interactively; all parameters are real valued and the
    usual magnitude of the range is 1; the names may not match any
    name in the list for \code{root}
    or any parameter name within the parameter
    list in \code{Kf},  since the parameters are included
    into the \code{root} list, which is passed to \code{beta} as a whole.}
  
  \item{root}{\{root growth\} for \code{plants.lambda} to \code{dir.ch.s}
    and\cr \{root, water uptake\} otherwise. \cr Any plant type is set to the
    values of \code{root} at starting time. 
    \code{root} is a list of the following components
    \describe{
      \item{\code{plants.lambda} \{mean \# of plants\}}{
	 the number of
	plants is Poisson distributed with mean \code{plants.lambda}}
      
      \item{\code{plants.mindist} \{min. plant. dist.\}}{
	 plants of the
	same species are at least \code{plants.mindist} away.  Actually,
	the algorithm tries to find random positions so that this
	condition is satisfied, but will give up if not successful after
	a few attemps. Then the value of \code{plants.mindist} is
	ignored.}
      
      \item{\code{mean} \{aimed mean total length\}}{
	 mean total length of
	all the roots of a single plant.  For a long time run, the
	average will be about the \code{mean}, but the algorithm is not
	too sophisticated to guarantee this.}
      
      \item{\code{sd} \{aimed sd of total length\}}{
	standard deviation for the total length of all the roots of a
	single plant.  There is a check by the algorithm to be
	close to \code{sd}}
     
      \item{\code{knot.probab} \{knot probability\}}{
	 probability of any root
	pixel to be a knot, except if the distance from the position to
	the previous knot is less than \code{knot.mindist}.  Then the
	probability for a knot is 0.}
      
      \item{\code{knot.mindist} \{min. knot distance\}}{
	 see \code{knot.probab}}
    
      \item{\code{shoots.3} \{3 shoots probab.\}}{
	 positive number, see
	\code{shoots.4}}
      
      \item{\code{shoots.4} \{4 shoots probab.\}}{
	 positive number.
	\describe{
	  \item{\code{shoots.4}\eqn{\ge}{>=}1}{
	    knots have always 4 shoots.}
	  
	  \item{\code{shoots.4}\eqn{<}1
	    and \code{shoots.4}+\code{shoots.3}\eqn{\ge}{>=}1}{
	    knots have with probability \code{shoots.4} 4 shots and 3
	    shoots otherwise.}
	  
	  \item{\code{shoots.4}\eqn{<}1
	    and \code{shoots.4}+\code{shoots.3}\eqn{\ge}{<}1}{
	    knots have 4, 3, 2 shoots with probability \code{shoots.4},
	    \code{shoots.3} and 1-\code{shoots.4}-\code{shoots.3},
	    respectively}
	}
      }
      
      \item{\code{stop.probab}}{
	\code{function(knot, dist, stop.m, stop.s)}
	that returns a values in \eqn{[0,1]}, which is the probability
	that a given root pixel does not continue to grow. 
	The function should be able to take matrices for the first 2
	components.  \code{knot} is the number of knots preceding and
	including the current pixel,
	\code{dist} is the distance to the surface along the root,
	\code{stop.m} and \code{stop.s} are arbitrary
	additional parameters.}
      
      \item{\code{stop.m} \{stop probability, m\}}{
	see \code{stop.probab}}
      
      \item{\code{stop.s} \{stop probability, s\}}{
	see \code{stop.probab}}
      
      \item{\code{rf.link}}{
	\code{function(v, rf.m, rf.s)}. 
	Summand in calculating the priority for all the neighbouring
	pixels of all active root pixels (end points of the roots).
	The pixels with the highest priority will be the next new root
	parts. 
	\code{rf.link} transforms the value of the simulated
	Gaussian random field where \code{rf.m} and \code{rf.s}
	are additional parameters.}
      \item{\code{rf.m} \{random field link, m\}}{
	see \code{rf.link}}
      
      \item{\code{rf.s} \{random field link, s\}}{
	see \code{rf.link}}
      
      \item{\code{no.own.root} \{own root penetration\}}{
	if \code{TRUE} a
	pixel may contain at most one root part of the same plant.}
      
      \item{\code{age.bonus} \{age bonus\}}{
	factor by which the priority number of a pixel is
	multiplied after each step.  The \code{age.bonus}
	is usually greater than or equal to 1.}
      
      \item{\code{depth.bonus} \{depth bonus\}}{
	 summand added to the priority number
	 if the subsequent pixel is below the current root pixel}
      
      \item{\code{side.bonus} \{side bonus\}}{
	 summand added to the
	priority if the subsequent pixel has the same depth than
	the current root pixel;  see \code{depth.bonus} for deeper
	points.}

     \item{\code{diagonal} \{diagonal directions\}}{
        if \code{TRUE} then roots may also grow into
       diagonal directions.
       }
      
      \item{\code{dir.ch} \{no direction change\}}{
	The number of pixels in x direction and the number of pixels
	in y direction  by which the new direction of the
	root growth deviates from the previous is added and gives a
	value \eqn{d} between 0 and 4.
      
	E.g., if the sequence of root pixels is (1,3), (2,4)
	then the previous direction is (1,1). If the pixel (2,3)
	is considered the new direction is (0, -1), so the modulus of
	deviation is 1 and 2 in x and y direction, respectively. So,
	\eqn{d=1+2=3}.
	
	Let \eqn{v=d \code{dir.ch}}{v = d * \code{dir.ch}}. 
	Then the summand for the priority of a pixel is given by normal
	random variable with mean v and \eqn{sd=|v| \code{dir.ch.s}}.}
     
      \item{\code{dir.ch.s} \{no dir. change, rel. sd\}}{
	 see \code{dir.ch}}
      
      \item{\code{rf.Kf} \{root, water uptake\}\{root changes field
	value\}}{ 
	logical.  If \code{TRUE} the value of the (Gaussian) random
	field at a root segment
	is changed according to the function \code{Kf}.}

      \item{\code{field.value} \{root, water uptake\}\{field.value\}}{
	The list entries \code{field.value} and \code{field.factor} are
	present only if the function \code{Kf} has the default
	definition.  In general, the named list elements of the
	parameter \code{param} of \code{Kf} are given.  See also
	\code{Kf}.}
      
      \item{\code{field.factor} \{root, water uptake\}\{field.factor\}}{
	see \code{field.value}}
       
      \item{\code{beta.value} \{root, water uptake\}\{beta.value\}}{
      	The list entries \code{beta.value} and \code{beta.factor} are
	present only if the function \code{beta} has the default
	definition.  In general, the named list elements of the
	parameter \code{param} of \code{beta} are given.  See also
	\code{beta}.} 
      
      \item{\code{beta.factor} \{root, water uptake\}\{beta.factor\}}{
	see \code{beta.value}}
       
      \item{\code{P0} \{root, water uptake\}\{P0\}}{
	(negative) pressure
	above which too less 
	oxygene to allow water uptake by plants}
      
      \item{\code{P2H} \{root, water uptake\}\{P2H\}}{
	\code{P2H}, \code{P2L},
	\code{r2H} and \code{r2L} determine, in dependence of the
	potential transpiration, the lower bound for optimal water
	uptake (piecewise linear curves);
	the upper bound of optimal water uptake
	is given by \code{materials$POptm}.
	See the SWMS2D manual for details.}
      
      \item{\code{P2L} \{root, water uptake\}\{P2L\}}{
	 see \code{P2H}}
      
      \item{\code{P3} \{root, water uptake\}\{P3\}}{
	 (negative) pressure below which
	no water uptake is possible}
      
      \item{\code{r2H} \{root, water uptake\}\{r2H\}}{
	see \code{P2H}}
      
      \item{\code{r2L} \{root, water uptake\}\{r2L\}}{
	 see \code{P2H}}
      
       \item{\code{root.condition} \{root, water uptake\}\{root water
	 uptake\}}{
	 boundary condition for the root.  Orginially, only
	 \sQuote{atmospheric} was allowed.
	}
      
	\item{\code{root.uptake} \{root, water uptake\}\{water uptake
	  value\}}{
	  value at the root in case of \sQuote{dirichlet} or
	  \sQuote{neumann} boundary condition
	}
    }
  }
  %
  \item{root.zone}{\code{NULL} or function of one variable.
    If not \code{NULL} then \code{rootzone} defines for any given
    \eqn{x} coordinate the interval for the height of the rootzone.
    A plant of species 1 is defined so that  
    all segments within the so defined root zone contains a root of this
    plant.
  }
  \item{col.rf}{sequence of colours that are used to show the
    random field; if \code{NULL} then some internally
    defined color sequence is used}
  \item{col.simu}{
    sequence of colours that are used to present any
    simulation result by swms2d; if \code{NULL} then some internally
    defined color sequence is used
  }
  \item{MaxIt}{
    maximum number of iteration during any time step.}

  \item{hTab}{
    \code{c(hTab1, hTabN)}, interval of pressure heads
    within which a table of hydraulic properties is generated.}

  \item{DMul}{
    \code{c(dMul, dMul2)},
    \eqn{\code{dMul}\ge 1}{\code{dMul >= 1}},
    \eqn{\code{dMul2}\le 1}{\code{dMul2 <= 1}};
    if the number of required iterations is less than 4 or greater than 6
    then the next time
    step is multiplied by \code{dMul} and \code{dMul2}, respectively.}
  %
  \item{col.rect}{
    colour of the button for free input; see
    \command{\link[RandomFields]{eval.parameters}}}

  \item{col.bg}{
    colour of a interactive bar; see
    \command{\link[RandomFields]{eval.parameters}}}

  \item{col.sep}{
    colour of the separating line; see
    \command{\link[RandomFields]{eval.parameters}}}

  \item{col.left}{
    colour of preceding element; see
    \command{\link[RandomFields]{eval.parameters}}}

  \item{col.mid}{
    colour for the message; see
    \command{\link[RandomFields]{eval.parameters}}}  

  \item{col.right}{
    colour of subsequent element; see
    \command{\link[RandomFields]{eval.parameters}}}

  \item{col.line}{
    colour of the marking line in interactive bars of
    absolute choice; see \command{\link[RandomFields]{eval.parameters}}}

  \item{col.txt}{
    colour of headers; see
    \command{\link[RandomFields]{eval.parameters}}}

  \item{col.submenue}{
    colour for the text in the main menue}
  
  \item{col.subord}{
    colour for the text in the  title of a secondary menue}

  \item{col.forbid}{colour for
    text that indicate forbidden main menues, e.g. \sQuote{polygon} after
    having defined already the maximum number of allowed horizons}

  \item{col.bg.forbid}{background colour for forbidden menue points}

  \item{col.flash}{
    colour for hints or the titles of active windows}
 
  \item{col.hor}{vector of 10 components.
    Colours that are used to draw the horizons}

  \item{col.draw}{
    colour for drawing or showing selected horizons}

  \item{col.mess}{
    colour for warning or error messages}

  \item{col.mesh}{
    colour for the edges of the finite element mesh}

  \item{col.mesh.pt}{
    colour for the vertices of the finite element mesh}
  %
  \item{col.exception}{vector of 2 components; colour plotted values
    that are below or above the given range of \code{zlim},
    see \code{lim.swms2s} and \code{zlim}.}
  \item{cex.eval}{value for the parameters \code{cex} and \code{cex.i}
    of \command{\link[RandomFields]{eval.parameters}} and
    for the parameter \code{cex.names} in
    \command{\link[RandomFields]{ShowModels}}.
  }
  \item{areas}{
    logical. If \code{TRUE} the horizons are coloured.  If
    \code{FALSE} only the border lines between the horizons are given.}

  \item{PrintLevel}{
    If \eqn{\le}{<=}0 nothing is printed. The higher the
    value the more information is given.}

  \item{new}{
    logical or \code{NULL}.  If \code{TRUE} the interactive plot is
    opened in a new window.  If \code{NULL} then the interactive window
    is not opened and the standard definitions are returned.}

  \item{X11.width}{
    numeric.  Width of the window.  Only used if
    \code{new=TRUE}.} 

  \item{X11.height}{numeric.  Height of the window.  Only used if
    \code{new=TRUE}.} 

  \item{frequent.reset}{logical.
    Background: the current implementation of \command{xswms2d}
    uses \command{\link[graphics]{split.screen}} to build the interactive plot;
    \command{\link[graphics]{split.screen}} memorises all changings in the plot and
    replays them all if the X11 window has to be rebuilt.
    If \code{frequent.reset=TRUE} then the memory effect is reduced
    at the cost of some flickering of the X11 window in the main menu.
  }
  
  \item{update}{
    \{updating\} logical.  If \code{TRUE} the simulations are
    updated after any changing of the parameters.

    Independently of the value of \code{update} an update of the
    simulation is performed if the user clicks on any of the given
    titles or subtitles of the submenu.
  }

  \item{waterflow}{
    \{water flow\} logical.  If \code{TRUE} swms2d is run}
  \item{zlim}{vector of 2 components or \code{NULL} which gives the
    limits for the plotted values in the \command{\link[graphics]{image}}s. If
    \code{NULL} then \code{zlim} is calculated internally, cf.
    \command{\link{plotWater}}.
  }
  \item{Zlim}{vector of 2 components or \code{NULL} which gives the
    limits for the plotted values in \command{\link[RandomFields]{ShowModels}}.
  }
  %
  \item{print.par}{
    \{postscript\} list of the following components
    \describe{
      \item{\code{ps} \{ps base name\}}{
	 character. Base name of the
	postscript file.}
      
      \item{\code{height} \{picture height\}}{
	 numeric. Height of the
	postscript figure in inches.}
    
      \item{\code{title} \{plot title\} }{
	logical. If \code{TRUE} a
	title is plotted.}

      \item{\code{legend} \{legend\} }{
	logical. If \code{TRUE} a legend is added to the plot.
      }
    }
  }
  %
  \item{...}{\command{sophy} is alias for \command{xswms2d} and
    takes the same arguments as \command{xswms2d}.}
}
\details{
  \bold{The interactive plot} is organised as follows:
  \describe{
    \item{left top}{drawing area: here, a soil profile is shown
      (if available, see parameter \code{picture})
      and the user draws the boundaries of the horizons by
      a sequence of left mouse button clicks.}
    \item{left bottom}{
        Random field of hydraulic conductivity.  That is, the 
        Gaussian random field including stones and roots, to which
	\code{miller.link} and then \code{millerK} has been applied. 
        (In the standard setting, the modulus of the 
        Gaussian random field is anti-proportional to the root of the
	hydraulic conductivity \eqn{K} and proportional to 
	the initial
        water potential \eqn{H}.)}
    \item{centre top}{Main menue; this area is also used for defining the
	parameters of the Gaussian random fields}
    \item{centre bottom}{plotting area for the water pressure \eqn{H},
        the charge rate \eqn{Q}, water contents \eqn{theta}, the flux vector
	components or the solute concentration as
        calculated by swms2d for coarsened grids.}
    \item{right}{area used by submenues}
    \item{right bottom}{plotting area for the water pressure etc for the
      orginal, fine grid, see main menue point \sQuote{precise
	waterflow}
      below}
  }

  \bold{The menue points} of the interactive plot are
  \describe{
    \item{postscript}{plotting
      the water potential \eqn{H},
      discharge/recharge rates \eqn{Q} for internal sink/sources,
      the water contents \eqn{theta},
      the x-components of the Darcian flux vector vx,
      the z-components of the Darcian flux vector vz,
      or
      the concentration;
      in grey it is noted whether the
      \sQuote{precise water flow} or the result of the approximate
      \sQuote{water flow}
      is given (according to the main menu point
      \sQuote{precise waterflow}). 
      Further the random field might be plotted. 
      For the additional parameters, see the input variable
      \code{print.par}.}
    
    \item{polygon}{A sequence of left mouse button clicks
      in the upper left
      drawing area defines the vertices of the polygon.  Up to 512
      vertices are allowed.  The right mouse button terminates the
      drawing.
      A polygon must consist of at least three points.  The points may
      be placed such that the connecting lines cross.}
    
    \item{horizon}{A sequence of left mouse button clicks in the upper
      left drawing area defines the border of the horizon.  Up to 512
      vertices are allowed. The right mouse button
      terminates the drawing. 
      A sequence of chosen points is only accepted to define a horizon
      if the first point has a smaller x-coordinate than all the other
      points 
      and the last point has a greater x-coordinate than all the other
      points;
      except for this restriction and the fact that at least 2 points
      must be given, any sequence of points is allowed, points may even
      leave the plotting area. 
      If the first or the last point is within the area then the horizon
      boundary is linearly extrapolated by a line through the first and
      second or the two last points, respectively.
      The new horizon gets as default parameters for the stones the
      input stone parameters, and for the material the material
      parameters of the preceding horizon.}
    
    \item{structure}{Definition of the random fields that underly the
      Miller similar medium. If more than one horizons or at least one
      polygon is defined the user selects first the part the user likes
      to define by a click into the appropriate part of the drawing
      area.
      Afterwards a new menue is opened in the upper right part of the
      window.  Within this new menue we have, see 
      \command{\link[RandomFields]{ShowModels}} for details,
      
      \describe{
	\item{bottom left}{a simulation of the Gaussian random field
	  with the specified parameters.}
	
	\item{top left}{the covariance function or the variogram for
	  the random field and the the transformed random field according
	  to \code{miller.link}.  If the \code{model} is anisotrope or
	  \code{anisotrop=TRUE} is given explicitely, the models are shown
	  for the two main axis.}
	
	\item{right}{Interactive menue. 
	  From top to bottom:
	  \itemize{
	    \item \sQuote{simulate} is given 
	    if \sQuote{updating} of the
	    main menue is \code{FALSE}
	    \item  name of the model
	    \item the formula.  Note the formula is not given
  	    for all the coded variogram models.
	    \item the specific parameters if there are any
	    \item variance
	    \item nugget effect
	    \item the scale parameter or the anisotropy paramters
	    \item mean.
	    \item If \sQuote{practical range}=\code{TRUE} then the
	    (isotropic) covariance  model is
	    rescaled such that at distance 1 the value of the model is
	    approximately 0.05.
	    \item \sQuote{variogram} switches between the
	    presentation of the variogram and the corresponding covariance
	    function.
	    \item \sQuote{variogram angle} gives the angle for the
	    direction of the main variogram. If \sQuote{variogram
	      angle}\code{=NA} the angle equals the angle of the anisotropy
	    specification.
	    \item If \sQuote{show transformed field}=\code{TRUE} then
	    the same transformation as for the hydraulic conductivity
	    parameter \eqn{\alpha_k}{alpha_k} is applied, i.e.,
	    \code{millerK(miller.link( ))}
	  }
	}
      }
     
      This menue is left by the right mouse button. A menue for
      choosing another variogram model is reached. Again,
      the right mouse button 
      leaves the menue.
      
      The first horizon gets as default the model and the
      parameter values passed to \code{xswms2d} by \code{model} and
      \code{anisotropy}.
      For the following horizons the values of the previously considered
      horizon are taken over, when considered first.}
    
    \item{stones}{definition of the stones}
    
    \item{root growth}{definition of the root growth; if
      \code{plant.types} is 
      greater than 1, the user first chooses among the different types,
      then enters
      the submenue. \sQuote{plant type} is used to name a species only,
      and
      is used nowhere else. The decription of all the other variables
      (\sQuote{mean \# of plants},...,\sQuote{no dir. change, rel. sd})
      is given above; see
      the input variable \code{root}.}
    
    \item{material[phys]}{Definition of the physical material constants,
      see variable \code{materials}.}

    \item{material[chem]}{Definition of the chemical material constants,
      see variable \code{materials}.}
    
    \item{root water uptake}{
      Definition of the water uptake by the roots.  See the variable
      \code{root}.}
    
    \item{atmosphere, data}{
      See the variable \code{atm.data} for a description.
      If \code{nrow(atm.data)} or \code{atm.periods} is greater than 1
      then the user chooses the atmospherical period first.}
    
    \item{swms2d [water]}{see the variable \code{water} for a description}

    \item{swms2d [chem]}{see the variable \code{chemical} for a
      description}
    
    \item{atmosphere control}{see the variable \code{atmosphere}
      for a description}
    
    \item{new simulation}{If parameters are changed, simulations are
      redone reusing previous simulations
      as much as possible.  That is, if a
      parameter of SWMS2D is changed the simulation of the Gaussian
      random field is reused.  Further the redone simulations are based
      on the old random seeds.
      
      Here, a simulation is redone from scratch, based on a new random
      seed.}
    
    \item{precise waterflow}{By \sQuote{swms2d water},
      submenue \sQuote{size reduction} (variable \code{water$red}) the
      finite 
      element mesh can be coarsened for faster simulations.  If the
      coarsening is genuine, i.e. \code{water$red}\eqn{>}1, then this
      menue point is active. If pressed, a simulation is performed on
      the original finite element mesh.}

    \item{water flow}{yes/no.  If \sQuote{no}, only the Gaussian random field,
      the stones and the roots are simulated.}

    \item{updating}{yes/no.  If yes the simulation is redone after any
      changing of any parameter, except for \sQuote{structure} definition,
      where the new simulation is performed when the submenue
      \sQuote{structure} is left.}

    \item{end}{terminates the interactive surface}
  }
 
  \bold{General information}
  \itemize{
    \item
    Currently, the number of horizons is limited to \code{h$max.horizon}=10.

    \item
    The horizons are build up from the bottom, i.e. first, the
    boundary for the C horizon should be drawn then those for the B
    horizons, etc.

    \item
    If the mean of the last genuine(!)
    horizon (polygons excluded) is defined \code{NA} within the menue
    \sQuote{structure}, then this part of the profile is interpreted as
    air and 
    the previously defined horizons are genuine ones.
    Any other horizon but the last one should not have mean
    \code{NA}.

    \item
    Polygons are used to define large lenses or to circumvent
    the predefined ordering of horizons.

    \item
    The ordering (as the sequence in which the horizons and
    polygons are defined, not their position in the profile)
    of the horizons and polygons is used in the conditional simulation
    of the Gaussian random fields, where conditioning happens with
    respect to preceding horizons.

    \item
    First the Gaussian random fields are simulated, starting with
    the first horizon.  After simulation of the complete profile, the
    stones are simulated, then the roots.
    If \sQuote{water flow}=\sQuote{yes}
    then the SWMS2D simulation is performed.  See
    \command{\link{simulateHorizons}}
    and the references therein for further details on the stochastic
    simulation.  See the SWMS2D manuscript for details on SWMS2D.

    \item
    If a parameter was changed, usually only a part of the
    simulations is redone, namely that level of simulation that
    corresponds to the parameter (random field for a specific horizon,
    stones, roots, SWMS2D) and all subsequent levels.
    For example, if a stone parameter has been changed, the simulation of
    the stones and the roots is done, and SWMS2D is performed. 
  }
}
\note{
 % In the following \{...\} indicates the name of the
 %   corresponding button (or the sequence of buttons)
 %   in the interactive plot. Generally the
 %   buttons for the elements of lists can be found in submenus.
 %   The button to enter the submenu is given at the beginning of the list.
 The phrases in brackets in the argument section of this documents
  give the respective menue points of the interactive plot.}
\value{
  \code{sophy} is a synonym for \code{xswms2d}
  A list is returned, where
  the first 10 components contain the definitions of the
  horizons; these 10 components
  are called by their position in the list.  Any other list component
  does not have a predefined location and may be called only by name.
  
  Note that if \command{sophy} or \command{xswms2d} is called with
  parameter \code{new=NULL} the interactive window is not opened, but
  the initial list (of \code{h}) is returned immediately.

  \item{1:10}{
    Each of the first 10 components is a list of the following
    elements: 
    \describe{
      \item{\code{type}}{
        character. Either \sQuote{Start}, \sQuote{H}, or \sQuote{P}.
	Exactly the 
        first element has the value \sQuote{Start}, all the others may be 
        horizons (\sQuote{H}) or polygons (\sQuote{P}).}
      
      \item{\code{cut.x}, \code{cut.y}}{
        \code{cut.x} and \code{cut.y} define the horziontal and
        vertical position of a rectangle, respective, which encloses the
        whole horizon or polygon.  \code{cut.x} and \code{cut.y} should be
        as tight as possible, to accelerate the simulation of the random
        fields for the horizons.}
            
      \item{\code{stone}}{
        list of stone parameters, see the input variable
        \code{stone}.  The input parameters are the default parameters for
        all stone specifications}
      
      \item{\code{materials}}{
        list of material constants, see the input variable
        \code{materials}.  The default parameters for the first horizon
        are the input parameters;  the default parameters for any other
        horizon (or polygon)
        are the parameters of the preceding horizon (or polygon)}
      
      \item{\code{model}}{
        list or \code{NULL}. Covariance model in list format, see
        \command{\link[RandomFields]{CovarianceFct}}.
        The \code{model} is \code{NULL} if \sQuote{structure} of the
	main menue 
        has not been called yet for the respective horizon. The
	\sQuote{Start} 
        horizon gets the value of the input parameter at starting time.
        If \sQuote{undo} is called when only the \sQuote{Start} horizon
	exists then 
        \code{model} is set to \code{NULL} even for the \sQuote{Start}
	horizon.} 
      
      \item{\code{border}}{
        matrix of two columns or \code{NULL}.
        boundary points between two horizons (or polygons); used for the
        conditional simulation of the Gaussian random fields.
        Further, if for the selection of a horizon the user clicks on such
        a boundary points, the user is warned. 
        \code{border} is \code{NULL} if
	only the \sQuote{Start} horizon exists.}
      
      \item{\code{idx}}{
        matrix of logical values with \code{diff(cut.x)+1} rows
        and \code{diff(cut.y)+1} columns.  \code{idx} indicates which
        pixels of the \code{cut.x }\eqn{\times}{ x }\code{cut.y} area belongs
        to the present horizon.}
    }
  }
  
  \item{grid.x}{\code{seq(xlim[1], xlim[2], step)}}

  \item{grid.y}{\code{seq(ylim[1], xlim[2], step)}}

  \item{idx.rf}{integer. Matrix of size
    \code{length(grid.x) }\eqn{\times}{x }\code{length(grid.y)}
    indicates the affiliation to a horizon or polygon (number is
    decremented by 1, i.e., \eqn{0,...,\code{h$n}-1}), see also
    the output variable \code{border} within the specific horizons.
    
    It allows also for the indication of the presence of a stone; then,
    if a stone is present, the horizon number is
    incremented by  \code{h$max.horizon}.
    The modulus by \code{h$max.horizon} gives the by 1 decremented
    number of the horizon.}
  
  \item{n}{current number of horizons}

  \item{step}{input parameter \code{step}}

  \item{max.horizons}{maximum number of horizons; currently 10.
                      The value of this variable may never be changed!}

  \item{beta}{input parameter \code{beta}}

  \item{root}{list, where \code{length(root)} equals the
    input parameter \code{plant.types}.
    Each component of the list is a list with the same components as
    the input parameter \code{root}. Additionally, it has the
    components
     \itemize{
       \item
       \code{plant.type}, a user defined name of the species.

       \item
       the parameters within the \code{param} list of the function
       \code{Kf}

       \item
        the parameters within the \code{param} list of the function
       \code{beta}
     }
   }
  
  \item{atmophere}{input parameter \code{atmosphere}}

  \item{atm.data}{matrix of 10 columns, corresponding to the input
    parameter \code{atm.data} and \code{atm.periods}.}

  \item{water}{input parameter \code{water}}

  \item{chem}{input parameter \code{chemical}}

  \item{Kf}{input parameter \code{Kf}}

  \item{RF}{matrix of size
    \code{length(grid.x) }\eqn{\times}{x }\code{length(grid.y)} or
    \code{NULL}. 
    \code{RF} contains the Gaussian random field resulting from all
    definitions for the horizons.
    \code{RF} is \code{NULL} if
    \code{simulate.horizon} has not been called yet.}

  \item{Stones.RF}{matrix of size
    \code{length(grid.x) }\eqn{\times}{x }\code{length(grid.y)} or
    \code{NULL}. 
    The value of \code{Stones.RF} is
    the Gaussian random field modified by the
    stones (concerning hydraulic conductivity).
    Stones have usually value \code{NA}, whereas air has value
    \code{NaN}.  
    
    \code{Stones.RF} is \code{NULL} if
    \code{simulate.horizon} has not been called yet or terminated
    with an error. Further \code{Stones.RF} is set to \code{NULL} if
    \code{materials$sharpness} has been changed, or \sQuote{structure},
    \sQuote{undo}, \sQuote{horizon}, or \sQuote{polygon} has
    been called.}

  \item{Root.RF}{matrix of size
    \code{length(grid.x)}\eqn{\times}{x}\code{length(grid.y)} or
    \code{NULL}. 
    The value of \code{Root.RF} is the value of \code{Stones.RF}
    modified by the roots (concerning hydraulic conductivity).
    That is, if \code{root$rf.Kf} is \code{TRUE}
    the values at the respecitve root segments are given by
    \code{h$KF}.
    
    \code{Root.RF} is \code{NULL} if
    \code{simulate.horizon} has not been called yet or terminated
    with an error. Further \code{Root.RF} is set to \code{NULL} if
    \code{materials$sharpness} has been changed, or \sQuote{structure},
    \sQuote{undo}, \sQuote{horizon}, or \sQuote{polygon} has
    been called.}

  \item{m.link}{link function used in the submenue \sQuote{structure},
    currently the function is identical to \code{miller.link}}

  \item{miller.link}{input parameter \code{miller.link}}

  \item{col.rf}{input parameter \code{col.rf} if not \code{NULL} and
    the internally defined sequence of colors otherwise}
  
  \item{col.simu}{input parameter \code{col.simu} if not \code{NULL} and
    the internally defined sequence of colors otherwise
  }

  \item{random.seed}{list or \code{NULL}.
    Each element contains the random seed for the
    simulation of the Gaussian random field in the respective horizon.
    \code{random.seed} is \code{NULL} if
    \command{\link{simulateHorizons}}
    has not been called yet.
  }
  \item{stone.random.seed}{The random seed for the simulation of the
    stones} 
  \item{root.random.seed}{The random seed for the simulaton of the roots}

  \item{rf.complete}{logical. \code{TRUE} if the stochastic simulation
    has been performed completely.}

  \item{plants}{list or \code{NULL}.
    Each component contains the information on the roots of
    a single plant.  Each component is a matrix of 8 columns and rows
    according to the number of root segments.
    The 8 columns are
    \enumerate{
      \item
      horizontal coordinate of a root segment.
      The coordinate is given in pixel units. 

      \item
      vertical coordinate of a root segment. 
      For the first segment it equals
      1 plus the vertical position of the first
      pixel that is not air at the chosen horizontal position.

      \item
      row index for the parent root pixel

      \item
      row index where the subsequent root pixel is found; \code{NA} if
      it is a terminating root segment. If the pixel is a knot with k
      children then the value gives the position of the first child.
      The other children are given in consecutive rows.

      \item
      number of children in case it is a knot, 1 otherwise

      \item
      the number of preceding knots until and including this position

      \item
      aimed random total length of the roots if it is the first
      pixel, and the current distance to the surface along the roots,
      otherwise.

      \item distance to the preceding knot
    }
    \code{plants} is \code{NULL} if \command{\link{simulateHorizons}}
    has not been called yet, or any parameter in \code{root}
    has been changed.
  }
  \item{plants.idx}{vector or \code{NULL}.  \code{plant.idx[i]} gives the
    plant type (\eqn{1,\ldots,\code{plant.types}}{1,...,\code{plant.types}})
    of the \code{i}th simulated plant.
    \code{plants.idx} is \code{NULL} if
    \command{\link{create.roots}} or \command{\link{simulateHorizons}}
    has not been called yet.}
  
  \item{water.x}{the coordinates of \code{grid.x} after thinning by
    factor \code{water$red}}
  
  \item{water.y}{the coordinates of \code{grid.y} after thinning by
    factor \code{water$red}}
  
  \item{flux}{logical vector of
    \code{length(water.x) }\eqn{\times}{x }\code{length(water.y)} elements. 
    \code{flux} indicates which pixels of the
    \code{length(water.x)}\eqn{\times}{x}\code{length(water.y)}
    grid are used in the SWMS2D simulation.  (Pixel are left out
    if they indicate stones (\code{NA}) or air (\code{NaN}).)}
  
  \item{swms2d}{three dimensional array or \code{NULL}. 
    The first dimension has 6 components which are
    \enumerate{
      \item the water potential \eqn{H}
      \item discharge/recharge rates \eqn{Q} for internal sink/sources
     \item the water contents \eqn{theta}
      \item the x-components of the Darcian flux vector vx
      \item the z-components of the Darcian flux vector vz
      \item the concentration
    }
    The second component corresponds to the sequence of pixels
    used in the SWMS2D simulation and has length \code{sum(flux)}.
    
    \code{swms2s} is \code{NULL} if
    the water flux has not been simulated yet
    because no simulation has been performed or one of the stochastic
    simulations (Gaussian random field, stones, roots) failed,
    or \sQuote{water flow} has been
    set to \sQuote{no}. Further \code{swms2d} is set to \code{NULL}
    if an error occurs when trying to plot
    the SWMS2D simulation results (by \command{\link{plotWater}}); this
    happens if the result has no finite values.
    Finally, \code{swms2d} is set to \code{NULL} if a relevant
    parameters was changed, such as those in
    \code{atmosphere}, \code{atm.data}, \code{materials}, \code{chem}
    or \code{water}.
  }
  
  Any of the lists that contain input parameters may have the additional
  component \code{.history}, an internal variable used by the
  interactive menue algorithm that gives the last entries by the user;
  see \command{\link[RandomFields]{evalpar}}. 
}
\section{Warning}{The user should consider the following parameters
  of \code{h[[i]]} as being read-only:
  \itemize{
    \item \code{cut.x}
    \item \code{cut.y}
    \item \code{border}
    \item \code{idx}
  }
  Further, \code{h$max.horizon} may never be changed.

}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/institute}}
\references{
  \itemize{
   \item
    Simunek., J., Vogel, T., and van Genuchten, M.Th. (1994)
    \emph{The SWMS2D code for simulating water flow and solute transport
      in two-dimensional variably saturated media, Version 1.21.}
    Research Report No. 132, 197 p.,
    U.S. Salinity Laboratory, USDA, ARS, Riverside, California.

    \item
    Stoyan, D., Kendall, W.S., Mecke, J. (1995)
    \emph{Stochastic Geometry and its Applications} Chichester: Wiley, 2nd
    edition.

    \item
    Schlather, M., Huwe, B. (2005)
    SOPHY - an interactive programme for water flow modelling.
    \emph{In preparation}
  }
}
\seealso{\command{\link{calculate.horizons}},
  \command{\link{create.roots}},
  \command{\link{create.stones}},
  \command{\link{create.waterflow}}, 
  \command{\link{draw.horizons}},
  \command{\link{modify.horizons}},
  \command{\link{plotRF}},
  \command{\link{plotWater}},
  \command{\link{simulateHorizons}},
  \code{\link{SoPhy}}, 
  \command{\link{swms2d}}
}

\examples{
%  library(SoPhy, lib="~/TMP"); 
% source("/home/schlather/R/SOPHY/SoPhy/R/simu.R")
% source("/home/schlather/R/SOPHY/SoPhy/R/Sophy.R")
% source("/home/schlather/R/SOPHY/SoPhy/R/swms2d.R")
% source("/home/schlather/R/NEW.RF/RandomFields/R/ShowModels.R");


% for unknown reasons the check runs for a too long time
\dontrun{
### without underlying profile
h0 <- xswms2d(xlim=c(1, 100), ylim=c(1, 100), step=1)
# or, equivalently:
# h1 <- sophy(xlim=c(1, 100), ylim=c(1, 100), step=1)

% for unknown reasons the check runs for a too long time
### underlying profile
### the profile was taken by M. Flury, J. Leuenberger,
### B. Studer, W.A. Jury and H. Fl\"uhler, see also URL
### \url{http://www.ito.umnw.ethz.ch/SoilPhys/Fliessmuster/projekt_flury.html}
pic <- paste(system.file(package='SoPhy'), 'tracer', 'K06', sep="/")
h <- xswms2d(xlim=c(0,160), step=1, aniso=TRUE,
             update=FALSE, waterflow=FALSE, pict=pic)

### repeated call
h <- xswms2d(h=h, update=TRUE, waterflow=TRUE)

% for unknown reasons the check runs for a too long time
### an example for non-atmospheric root uptake
h <- sophy(xl=c(1,10), yl=c(1,10), step=0.1, new=NULL) ## get standard
h$root[[1]]$root.condition <- 1   ## dirichlet condition for roots
h$root[[1]]$root.uptake <- -200   ## value for dirichlet condition
h$root[[1]]$plants.lambda <- 0.07 ## intensity of plants
h$water$TPrint <- 2     ## end point of simulation (film lasts about 2 min)
h$water$red <- 1        ## precise simulation
h$water$breakpoint <- 3 ## high image frequency (close to a film)
h$water$top.bound <- 2  ## dirichlet
h$water$top.value <- -8 ## value on the boundary
h <- xswms2d(h=h, update=TRUE, waterflow=TRUE)
}

% source("/home/schlather/R/RF/RandomFields/R/evalpar.R"); source("/home/schlather/R/RF/RandomFields/R/ShowModels.R"); source("/home/schlather/R/SOPHY/SoPhy/R/simu.R"); source("/home/schlather/R/SOPHY/SoPhy/R/Sophy.R");try(repeat dev.off()); h1 <- sophy(xlim=c(1, 100), ylim=c(1, 100), step=1)
    
% source("/home/schlather/R/RF/RandomFields/R/evalpar.R"); source("/home/schlather/R/RF/RandomFields/R/ShowModels.R"); source("/home/schlather/R/SOPHY/SoPhy/R/simu.R"); source("/home/schlather/R/SOPHY/SoPhy/R/Sophy.R");try(repeat dev.off());h <- xswms2d(xlim=c(0,160), step=1, aniso=TRUE, update=FALSE, waterflow=FALSE, pict=pic)

}

\keyword{spatial}

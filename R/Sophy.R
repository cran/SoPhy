### inconsistent: h kann nicht durchgezogen werden!



sophy <- function(...) xswms2d(...)

xswms2d <- 
  function(h, xlim, ylim, step, picture = NULL,

           ## material and water flux parameters
           water = list(
             ### printing
             TPrint = 500,
             red = 4,
             print = 7,
             mesh = TRUE,
             lim.rf = 1, ## upper limit for z in plotRF
             lim.swms2d = 1, ## upper limit for z in plotWater
             #mesh = FALSE,
	     ### simulation control
             TolTh = 0.0001,
             TolH  = 0.01,
             lWat = TRUE,
             dt = 1,
             dtMinMax = c(0.01, 60),
             max.iteration = 1000,
             max.iter.prec = 50000,
             iter.print = 100,
             breakpoint = 100,
	     ### boundary conditions
             top.bound = 2,
             top.value = 0,
             bottom.bound = 1,
             bottom.value = 0
             ),
           materials=list(thr=.02, ths=0.35, tha=0.02, thm=0.35, 
             Alfa=0.041, n=1.964, Ks=0.000722, Kk=0.000695, thk=0.2875,
             first=1, second=1, angle=0, ## geometry for conductivity
             Hslope=0, Hseg=-100, POptm=-25, sharpness=1,
             
             Bulk.d=1500, Diffus=0, longiDisper=1, transvDisper=0.5,
             Adsorp=0.0004, SinkL1=-0.01, SinkS1=-0.01, SinkL0=0, SinkS0=0
             ),
           Hinit=function(Hseg, depth) Hseg,
           model= list(model="exp", param=c(0, 0.25, 0, diff(ylim) * 0.1)),
           anisotropy = NULL,
           miller.link = function(rf) exp(sign(rf) * sqrt(abs(rf))),
           millerH = function(lambda) lambda, # Zavattaro et al, SSSA 1999;
           millerK = function(lambda) lambda^-2, # SWMS, p. 12
           millerT = function(lambda) 1,

           ## chemical data
           chemical=list(
             lChem = FALSE,             
             Epsi=2, ## epsi is transformed in create.water!!!
             lUpW=1, ## lUpW is transformed in create.water!!!
             lArtD=FALSE, PeCr=10, tPulse=500,
             top.bound = 2,
             top.value = 1,
             bottom.bound = 1,
             bottom.value = 0,
             root.uptake = 0,
             intern.source = 0
             ),
   
           ## atmospherical data
           atmosphere = list(
             AtmInf = TRUE,
             tInit = 0, # starting point of simu. tAtm gives the END of a period
             ## except for the very last entry, where the end point is ignored
             Aqh = -.1687, ## 0 if qGWLF=false
             Bqh = -.02674,## 0 if qGWLF=false
             hCritS = 1.e30,
             GWL0L = ylim[2],
             rLen = diff(xlim)
             ),
           atm.periods = 1, 
           atm.data = c(
             end   = 100000000, # end point of period, for which the following
             ##           values are used
             ## !!!!!!!!!!! value must be greater than TPrint
             prec  = 0, # precipitation [0.5]
             cPrec = 0, # precipitation, solute, 0 if lChem=false
             rSoil = 0, # potential evapo
             rRoot = 0, # transpiration [0.3]
             hCritA= 1000000, # abs value of min allowed pressure head at surface
             rGWL  = 0, # drainage flux (for Kode = -3), 0 if no Kode(n)=-3
             GWL   = 0, # ground water level
             crt   = 0, # concentration, Kode=-3, KodCB < 0
             cht   = 0 # concentration, Kode=+/-3, KodCB > 0
             ),

           stone=list(
             ### basics
             value = NA,
             lambda = 0,
             ### location
             no.upper=FALSE,
             no.lower=TRUE,
             no.overlap=FALSE,
             ### size 
             main.distr=rnorm, main.mean=diff(ylim)/20, main.s=diff(ylim)/400,
             sec.distr=rnorm, sec.mean=diff(ylim)/50, sec.s=diff(ylim)/400,
             phi.distr=rnorm, phi.mean=1, phi.s=0
             ),
     
           ### roots
           plant.types = 1,
           Kf=function(plant, distances, depth, rf, 
                       param=list(field.value=1, field.factor=0))
             rep(param$field.value, length(rf)) + param$field.factor * rf,
           beta=function(plant, distances, depth, rf, 
		         param=list(beta.value=1, beta.factor=0))
             rep(param$beta.value, length(rf)) + param$beta.factor * rf,
           root=list(
             plants.lambda = 0,
             plants.mindist = diff(xlim)/20,	
             mean = sqrt(diff(ylim) * diff(xlim)),
             sd = sqrt(diff(ylim) * diff(xlim)) / 5,
           ### root growth
             knot.probab = 0.1,
             knot.mindist = 5 * step,
             shoots.3  = 0.4,
             shoots.4  = 0,
             stop.probab = function(knot, dist, m, s) 1 - exp(-m + s * dist),
             stop.m = 0,
             stop.s = 0,
             # rf.link = function(v, m, s) m * h$millerK(h$miller.link(s * v)),
             rf.link = function(v, m, s) {
               v <- s * v
               m * (exp(sign(v) * sqrt(abs(v))))^(-2)
             },
             rf.m = 0,
             rf.s = 0,
	  ### advanced root growth   
             no.own.root = TRUE,
             age.bonus = 1.2,
             depth.bonus = 5.2,
             side.bonus = 5.1,
             diagonal = TRUE,
             dir.ch = 3.71,
             dir.ch.s = 0.5,
          ### water uptake / change of conductivity
             rf.Kf=FALSE, # alternative value for random field at a root?
             ## ordering is  P3 (up) h3 (constant) POptim (down) P0
             ## where h3 is a function of the potential transpiration rate
             ## piecewise linear, depending on P2L, P2H, r2H and r2L
 	     P0=-10, 
	     P2H=-200,
	     P2L=-800, 
	     P3=-8000,
             r2H=0.5, 
	     r2L=0.1,
             root.condition = 3, # atmospheric
             root.uptake = -150
             ),
           root.zone = NULL,
           col.rf = NULL,
           col.simu = NULL,
           
	   ## control parameters for swms2d 
	   MaxIt = 20, 
	   hTab = c(0.001,200),
	   DMul = c(1.1, 0.33),
           
           ## menue, parameter input
           col.rect="red",
           col.bg="yellow",
           col.sep="gray", 
           col.left="red", col.mid="darkgreen", col.right="darkgreen",
           col.line ="red", 
           col.txt="black",

           col.submenue="darkgreen",
           col.subord="steelblue", # box colors both able.
           col.forbid = "gray88",
           col.bg.forbid="lightyellow",# disabled menue point; 
           col.flash="blue",
           
          ## drawing 
           col.hor=c("#000000", "#996600", "#660000",  "#CC9933",
             "#666600",  
             "#CCCC99", "#CCCCCC", "#990000", "#FFCC00", "#FFFFCC"),
           col.draw="green",
           ## messages
           col.mess = "red",
           ## images
           ## alternatively: 
           ##  col.rf = paste("#0000", HEX[1 + rev(l.col) / 16],
           ##     HEX[1 + rev(l.col) %% 16], sep=""),
           ## col.rf = rainbow(200),
           ## finite element mesh
           col.mesh="blue", col.mesh.pt="green",
           col.exception = c("brown", "lightgray"),
           cex.eval = 0.8,
           areas = is.null(h$picture), ## drawing picture as image or
           ##                             as line scetch?
           PrintLevel=RFparameters()$Print,
           new=TRUE,
           X11.width=9.5, X11.height=X11.width * 0.789,
           frequent.reset = TRUE,
           update = FALSE,
           waterflow = FALSE,
           zlim=NULL,
           Zlim=NULL,

           ## postscript output
           print.par = list(ps="sophy", height=3, titl=TRUE, legend=TRUE)
           )
{
  simu.dev <- draw.dev <- showmodel.dev <- menue.dev <- param.dev <-
    waterflow.dev <- message.dev <- prec.water.dev <- percolation.dev <-
      few.param.dev <- draw.lim <- NULL

  HEX <- c(0:9, LETTERS[1:6])
  l.col <- as.integer(c(seq(255, 230, len=15), 230:101, seq(100, 0, len=5)))
  l.col2 <- as.integer(c(255:100, seq(100, 0, len=50)))
  if (is.null(col.simu))
    col.simu <-
      c(paste("#", HEX[1 + l.col2 / 16], HEX[1 + l.col2 %% 16],
              HEX[1 + l.col2 / 16], HEX[1 + l.col2 %% 16], "FF", sep=""),
        paste("#0000", HEX[1 + l.col / 16], HEX[1 + l.col %% 16], sep="")[-1])
  if (is.null(col.rf))
    col.rf <- paste("#", HEX[1 + rev(l.col) / 16],
                    HEX[1 + rev(l.col) %% 16],HEX[1 + rev(l.col) / 24],
                    HEX[1 + rev(l.col) %% 16],"00", sep="")
  
  old.options <- options()[c("warn", "locatorBell")]
  options(warn=1, locatorBell=FALSE)
  nlocator <- 512 ## maximum number of points by which a horizon is drawn

  if (!missing(h) &&!is.list(h)) stop("h must be a list!")
  
 
  atm.names <- c("end", "prec", "cPrec", "rSoil", "rRoot",
                   "hCritA", "rGWL","GWL", "crt", "cht")
  if (missing(h)) {
    h <- list()
    max.horizons <- RFparameters()$maxmodels
    h[[max.horizons + 5]] <- NA
    names(h) <- c(paste("H", 1:max.horizons,sep=""), 
                  "grid.x",
                  "grid.y",
                  "idx.rf",
                  "n",
                  "step")
    h$max.horizons <- max.horizons
    stopifnot(step>0)
    
    if (!missing(xlim)) { ## only allowed if picture != NULL
      stopifnot(length(xlim)==2)
      h$grid.x <- seq(xlim[1], xlim[2], step)
      xlim[2] <- h$grid.x[length(h$grid.x)]
    }
    if (!missing(ylim)) { ## only allowed if picture != NULL
      stopifnot(length(ylim)==2)
      h$grid.y <- seq(ylim[1], ylim[2], step)
      ylim[2] <- h$grid.y[length(h$grid.y)]
    }

#    if (.Platform$OS.type!="unix" && && !is.null(picture) &&
#         is.character(picture)) {
#      if (PrintLevel>0)
#        cat("\npictures can currently be loaded only under unix systems\n")
#      picture <- NULL
#    }
    if (!is.null(picture)) {
      if (is.character(picture))
        picture <- read.picture(picture, PrintLevel=PrintLevel)
      else stopifnot(is.array(picture))
      dp <- dim(picture)
      if (length(dp)==0) stop("file not found or invalid picture format")
      else stopifnot(length(dp)==3, dp[3] %in% c(3,4))
      ly <- length(h$grid.y)
      if ((lx <- length(h$grid.x))==0) {
        h$grid.x <- seq(0, diff(ylim) / dp[2] * dp[1] - step, step)
        xlim <- c(0, h$grid.x[length(h$grid.x)])
        lx <- length(h$grid.x)
      } else if (ly==0) {
        h$grid.y <- seq(0, diff(xlim) / dp[1] * dp[2] - step, step)
        ylim <- c(0, h$grid.y[length(h$grid.y)])
        ly <- length(h$grid.y)
      }     
      h$picture <-
        picture[1 + dp[1] * 0:(lx - 1) / lx, 1 + dp[2] * 0:(ly - 1)/ ly, ]
      picture <- NULL
    }
    
    if (length(h$grid.x) < 4)
      stop("grid in x direction must have at least 4 points")
    if (length(h$grid.y) < 4)
      stop("grid in y direction must have at least 4 points")
    h$step <- step
    h$n <- as.integer(1)    
    h$idx.rf <- as.double(NA)
    h$rf.complete <- FALSE
    
    
    stopifnot(length(eval(formals(Kf)))==5, length(eval(formals(beta)))==5)
    if (length(cc <- as.character(formals(Kf)$param))==0 || cc=="" )
      stop("argument 'param' in Kf may not be missing, neither its value")
    if (length(cc <- as.character(formals(beta)$param))==0 || cc=="" )
      stop("argument 'param' in beta may not be missing, neither its value")
    Kf.param <- eval(formals(Kf)$param)
    beta.param <- eval(formals(beta)$param)
    stopifnot(is.list(Kf.param), is.list(beta.param))
    stopifnot(is.function(Kf), is.function(beta), is.function(Hinit),
              is.function(root.zone) || is.null(root.zone),
              is.function(root$stop.probab), is.function(root$rf.link)
              )
    h$Kf <- Kf
    environment(h$Kf) <- EmptyEnv()
    h$beta <- beta
    environment(h$beta) <- EmptyEnv()
    h$Hinit <- Hinit
    environment(h$Hinit) <- EmptyEnv()
    h$root.zone <- root.zone
    if (!is.null(root.zone)) environment(h$root.zone) <- EmptyEnv()

    h$root <- list()
    for (j in 1:plant.types) {
      h$root[[j]] <- c(root, Kf.param, beta.param)
      h$root[[j]]$plant.type <- j
      environment(h$root[[j]]$stop.probab) <- 
        environment(h$root[[j]]$rf.link) <- EmptyEnv()	     
    }
    
    h$atmosphere <- atmosphere
    if (is.matrix(atm.data)) {
      stopifnot(ncol(atm.data)==10)
      h$atm.data <- atm.data
      dimnames(h$atm.data) <- list(NULL, atm.names)
      atm.periods <- nrow(h$atm.data)
    } else {
      stopifnot(is.vector(atm.data), length(atm.data)==10)
      h$atm.data <- t(matrix(atm.data, nrow=10, ncol=atm.periods,
                             dimnames=list(atm.names, NULL)))
    }

    
    h$water <- water
    h$chem <- chemical
    h[[1]] <- list(type="Start",
                   cut.x=as.integer(c(1,length(h$grid.x))),
                   cut.y=as.integer(c(1,length(h$grid.y))), 
                   stone=stone,
                   materials=materials
                   )
    model <- PrepareModel(model=model)
    m.aniso<- model$anisotropy
    model <- convert.to.readable(model, allowed="list")
    if (!is.null(anisotropy)) {
      if ((m.aniso)!=anisotropy) {
        if (anisotropy) {
          model$model[[1]]$aniso <- diag(1/model$model[[1]]$scale,2)
          model$model[[1]]$scale <- NULL
        } else {
          model$model[[1]]$scale <- model$model[[1]]$aniso[1]
          model$model[[1]]$aniso <- NULL
        }
      }
    }
    h[[1]]$model <- model
    h$miller.link <- if (is.null(miller.link))  function(x) x else miller.link
    h$millerH <- if (is.null(millerH)) function(x) x else millerH
    h$millerK <- if (is.null(millerK)) function(x) x^(-2) else millerK
    h$millerT <- if (is.null(millerT)) function(x) 1 else millerT
    stopifnot(is.function(h$miller.link),
              is.function(h$millerH),
              is.function(h$millerK),
              is.function(h$millerT))
    environment(h$miller.link) <- # environment(h$m.link) <-
      environment(h$millerH) <- environment(h$millerK) <- 
      environment(h$millerT) <- EmptyEnv()
    h$col.rf <- col.rf
    h$col.simu <- col.simu
  } else {
    xlim <- range(h$grid.x)
    ylim <- range(h$grid.y)
    if (length(h$grid.x)<4 || length(h$grid.y)<4)
      stop("grid must consist of at least 4 points in each direction")
    step <- h$grid.x[2] - h$grid.x[1]
    # step <- h$grid.y[2] - h$grid.y[1]
    plant.types <- length(h$root)    ## other possibility:
    ##            max(length(h$root), plant.types)   + initialisation
    atm.periods <- nrow(h$atm.data) ## other possibility:
    ##            max(nrow(h$atm.data), plant.types)   + initialisation
    model <- h[[1]]$model
    if (!all(sapply(h$root, is.list)))
      stop("root misspecification -- not a list of lists")
    environment(h$miller.link) <- # environment(h$m.link) <-
      environment(h$millerH) <- environment(h$millerK) <- 
      environment(h$millerT) <- EmptyEnv()
  }
 
   for (i in 1:h$n) {
     if (!is.null(h[[i]]$model))
       h[[i]]$model <-
         convert.to.readable(PrepareModel(h[[i]]$model), allowed="list")
  }
  
  m.link <- function(x) h$millerK(h$miller.link(x))## used in Showmodels and
  ##                                          postscript

  .Call("GetHorizons", h, n=as.integer(c(1,h$n)), PACKAGE="SoPhy")
  # h$RF <- matrix(NA, nrow=length(h$grid.x), ncol=length(h$grid.y)) # 31.12.03
  if (is.null(new)) return(h)

  
  ENVIR <- environment()
  materials.entry <-
    list(	  
         ## what="water" is enough if no change in sharpness
         ## at long term put sharpness somewhere else
         list(name="Material Constants", val="simulate", col=col.subord),
         list(name=expression(theta[r]), var="materials$thr", 
              delta=TRUE, val=function(d, v) {quadratic(d=d, v=v, a=1) }),
         list(name=expression(theta[s]), var="materials$ths", 
              delta=TRUE, val=function(d, v) {quadratic(d=d, v=v, a=1) }),
         list(name=expression(theta[a]), var="materials$tha", 
              delta=TRUE, val=function(d, v) {quadratic(d=d, v=v, a=0.1) }),
         list(name=expression(theta[m]), var="materials$thm", 
              delta=TRUE, val=function(d, v) {quadratic(d=d, v=v, a=0.1) }),
         list(name=expression(alpha), var="materials$Alfa", 
              delta=TRUE, val=function(d, v) {quadratic(d=d, v=v, a=0.1) }),
         list(name="n", var="materials$n", 
              delta=TRUE, val=function(d, v) {quadratic(d=d, v=v, a=1) }),
         list(name=expression(K[s]), var="materials$Ks", 
              delta=TRUE, val=function(d, v) {quadratic(d=d, v=v, a=0.001) }),
         list(name=expression(K[k]), var="materials$Kk", 
              delta=TRUE, val=function(d, v) {quadratic(d=d, v=v, a=0.001) }),
         list(name=expression(theta[k]), var="materials$thk", 
              delta=TRUE, val=function(d, v) {quadratic(d=d, v=v, a=0.001) }),
         ##
         list(name="Conductivity: Geometry", val="simulate", col=col.subord),
         list(name="1st principal comp. (scale??)", var="materials$first",
              delta=TRUE, val=function(d, v) {quadratic(d=d, v=v, a=10,
                            mini=0.001) }),
         list(name="2nd principal comp.", var="materials$second", delta=TRUE,
                val=function(d, v){quadratic(d=d, v=v, a=10,
                  mini=0.001) }),
         list(name="angle (degrees)", var="materials$angle",
              delta=FALSE, val=function(d, v) pmin(180, pmax(0,d * 180))),
         ##
         list(name="Other", val="simulate", col=col.subord),
         list(name="initial H, slope", var="materials$Hslope", 
              delta=TRUE, val=function(d, v) {quadratic(d=d, v=v, a=100,
                            mini=-Inf, maxi=0) }),
         list(name="initial H, segment", var="materials$Hseg", 
              delta=TRUE, val=function(d, v) {quadratic(d=d, v=v, a=100,
                            mini=-Inf) }),
         list(name="POptm (root)", var="materials$POptm", 
              delta=TRUE, val=function(d, v) {quadratic(d=d, v=v, a=20,
                            maxi=0, mini=-Inf) }),
         list(name="sharpness", var="materials$sharpness", 
              delta=FALSE, val=function(d, v) pmin(1, pmax(0,d)),
              cond="type!='Start'", what="all")
         )
    chem.material.entry <-
    list(
         list(name="Chemical Constants", val="simulate", col=col.subord),
         list(name="bulk density", var="materials$Bulk.d", 
               delta=TRUE, val=function(d, v) {quadratic(d=d, v=v, a=1000) }),
         list(name="diffusion", var="materials$Diffus", 
              delta=TRUE, val=function(d, v) {quadratic(d=d, v=v, a=1) }),
         list(name="longitud. Dispers.", var="materials$longiDisper", 
              delta=TRUE, val=function(d, v) {quadratic(d=d, v=v, a=1) }),  
         list(name="transvers. Dispers.", var="materials$ transvDisper", 
              delta=TRUE, val=function(d, v) {quadratic(d=d, v=v, a=1) }),
         list(name="Freundlich isoth. coeff.", var="materials$Adsorp", 
              delta=TRUE, val=function(d, v) {quadratic(d=d, v=v, a=0.001) }),
         list(name="1st-order const./dissolved", var="materials$SinkL1", 
              delta=TRUE,
              val=function(d, v) {quadratic(d=d, v=v, a=0.01, mini=-Inf) }),  
         list(name="1st-order const./solid", var="materials$SinkS1", 
              delta=TRUE,
              val=function(d, v) {quadratic(d=d, v=v, a=0.01, mini=-Inf) }),
         list(name="0-order const./dissolved", var="materials$SinkL0", 
              delta=TRUE,
              val=function(d, v) {quadratic(d=d, v=v, a=0.01, mini=-Inf) }),    
         list(name="0-order const./solid", var="materials$SinkS0", 
              delta=TRUE,
              val=function(d, v) {quadratic(d=d, v=v, a=0.01, mini=-Inf) })
          )

  stones.entry <-
    list(
         list(name="Stones, Basics", val="simulate", col=col.subord),
         list(name="value", var="stone$value",
              delta=TRUE, val=function(d, v) {
                quadratic(d=d, v=v, a=100) }),
         list(name="intensity", var="stone$lambda",
              delta=TRUE, val=function(d, v) {
                quadratic(d=d, v=v, a=0.01 / h$step^2) }),
         ##         
         list(name="Stone Location", val="simulate", col=col.subord),
         list(name="overlap into upper horizons", var="stone$no.upper",
              val=FALSE),
         list(name="overlap into lower horizons", var="stone$no.lower",
              val=FALSE),
         list(name="overlap of stones", var="stone$no.overlap", val=FALSE),
         ##
         list(name="Stone Shape", val="simulate", col=col.subord),
         list(name="main axis, mean", var="stone$main.mean",
              delta=TRUE, val=function(d, v) {
                quadratic(d=d, v=v, a=diff(range(h$grid.x)) / 10,
                          mini=h$step/100) }),
         list(name="main axis, sd", var="stone$main.s", 
              delta=TRUE, val=function(d, v) {
                quadratic(d=d, v=v, a=diff(range(h$grid.x)) / 50)}),
         list(name="second axis, mean", var="stone$sec.mean",
              delta=TRUE, val=function(d, v) {
                quadratic(d=d, v=v, a=diff(range(h$grid.x)) / 10,
                          mini=h$step/100) }),
         list(name="second axis, sd", var="stone$sec.s",
              delta=TRUE, val=function(d, v) {
                quadratic(d=d, v=v, a=diff(range(h$grid.x))/50) }),     
         list(name="phi, mean", var="stone$phi.mean",
              delta=TRUE, val=function(d, v) {
                quadratic(d=d, v=v, a=pi/2, maxi=pi) }),
         list(name="phi, sd", var="stone$phi.s",
              delta=TRUE, val=function(d, v) {
                quadratic(d=d, v=v, a=10) })
       )
  atmosphere.entry <- 
    list(list(name="atmosph. control param.", val="simulate", col=col.subord),
         list(name="use atmospheric data", var="AtmInf", val=TRUE),
         list(name="start time of simu.", var="tInit",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1000)),  
         list(name=expression(A[qh]), var="Aqh",
              delta=TRUE, val=function(d,v)quadratic(d=d, v=v,a=0.5,mini=-Inf)),
         list(name=expression(B[qh]), var="Bqh",
              delta=TRUE, val=function(d,v)quadratic(d=d, v=v,a=0.5,mini=-Inf)),
         list(name="max pressure at surface", var="hCritS",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1e30)), 
         list(name="ref. pos of groundwater", var="GWL0L",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=ylim[2]/5)),
         list(name="width of soil (transpiration)", var="rLen",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=diff(xlim)/5))
        )
  atm.data.entry <- 
    list(
         list(name="end point", var="end",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1000)),
         list(name="surface information", val="simulate", col=col.subord),
         list(name="precipitation", var="prec",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=0.01)),
         list(name="solute, precipitation", var="cPrec",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1)),
          list(name="potential evaporation", var="rSoil",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1)),
         list(name="potential transpiration", var="rRoot",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=0.005)),
         list(name="abs min. pressure at surface", var="hCritA",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1000000)),
         ##
         list(name="groundwater & drainage", val="simulate", col=col.subord),
         list(name="drainage flux (drain. or var. Q)", var="rGWL",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1)),
         list(name="ground water level", var="GWL",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1)),
          list(name="conc. flux (drainage)", var="crt",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1)),        
         list(name="conc. `pressure' (drain/var. H)", var="cht",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1))
     )
  #
  uptake.entry <-
    list(list(name="Kf parameters", val="simulate", col=col.subord),
         list(name="root changes field value", var="rf.Kf", val=TRUE)
         )
  if (!is.null(nm <- names(formals(h$Kf)$param))) {
    stopifnot(length(nm)>1)
    for (i in 2:length(nm)) {      
      uptake.entry <-
        c(uptake.entry,
          list(list(name=nm[i], var=nm[i],
                    delta=TRUE, val=function(d, v) quadratic(d=d, v=v, a=1,
                                  min=-Inf))))
    }
  }
  if (!is.null(nm <- names(formals(h$beta)$param))) {
    stopifnot(length(nm)>1)
    uptake.entry <-
      c(uptake.entry,
        list(list(name="params of fctn beta", val="simulate", col=col.subord)))
    for (i in 2:length(nm)) {      
      uptake.entry <-
        c(uptake.entry,
          list(list(name=nm[i], var=nm[i],
                    delta=TRUE, val=function(d, v) quadratic(d=d, v=v, a=1,
                                  min=-Inf))))
    }
  }

  uptake.entry <- 
    c(list(
         list(name="plant type", var="plant.type", what="none")
         ),
      uptake.entry,
      list(
           list(name="Water uptake curve", val="simulate", col=col.subord),
           list(name="P0", var="P0",
                delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=10, 
                              mini=-Inf, maxi=0)),
           list(name="P2H", var="P2H",
                delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=400, 
                              mini=-Inf, maxi=0)),
           list(name="P2L", var="P2L",
                delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1000, 
                              mini=-Inf, maxi=0)),
           list(name="P3", var="P3",
                delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=10000, 
                              mini=-Inf, maxi=0)),
           list(name="r2H", var="r2H",
                delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1)),
           list(name="r2L", var="r2L",
                delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1)),
           list(name="root water uptake", var="root.condition",
                val=c("dirichlet", "neumann", "atmospheric", "none")),
           list(name="water uptake value", var="root.uptake", 
                delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=20,
                              mini=-Inf))
           )
      )
        
  root.entry <-
    list(
         list(name="plant type", var="plant.type", what="none"),
         list(name="Roots, Basics", val="simulate", col=col.subord),
         list(name="mean # of plants / unit length", var="plants.lambda",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v,
                            a=10 / diff(xlim))),
         list(name="aimed min. plant dist. [unit length]", var="plants.mindist",
              delta=TRUE, val=function(d,v)quadratic(d=d, v=v,a=diff(xlim)/15)),
         list(name="aimed mean total length [unit length]", var="mean",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v,
                            a=5 * sqrt(diff(ylim) * diff(xlim)))),      
	 list(name="aimed sd of total length [unit length]", var="sd",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v,
                            a=sqrt(diff(ylim) * diff(xlim)))),
         ##
         list(name="Root Growth", val="simulate", col=col.subord),
         list(name="knot probability", var="knot.probab",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1)),
         list(name="min. knot distance [unit length]", var="knot.mindist",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=10 * step)),
         list(name="3 shoots probab.", var="shoots.3", 
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1, maxi=1)),
         list(name="4 shoots probab.", var="shoots.4", 
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1, maxi=1)),
         list(name="stop probability, m", var="stop.m",
              delta=TRUE, val=function(d,v)quadratic(d=d, v=v, a=1, mini=-Inf)),
         list(name="stop probability, s", var="stop.s",
              delta=TRUE, val=function(d,v)quadratic(d=d, v=v, a=1, mini=-Inf)),
         list(name="random field link, m", var="rf.m",
              delta=TRUE, val=function(d,v)quadratic(d=d, v=v, a=1, mini=-Inf)),
         list(name="random field link, s", var="rf.s",
              delta=TRUE, val=function(d,v)quadratic(d=d, v=v, a=1, mini=-Inf)),
         ##
         list(name="Advanced Root Growth", val="simulate", col=col.subord),
         list(name="own root penetration", var="no.own.root", val=FALSE),
         list(name="age bonus", var="age.bonus",
              delta=TRUE, val=function(d,v)quadratic(d=d, v=v, a=1, mini=-Inf)),
         list(name="depth bonus", var="depth.bonus",
              delta=TRUE, val=function(d,v)quadratic(d=d, v=v, a=5, mini=-Inf)),
         list(name="side bonus", var="side.bonus",
              delta=TRUE, val=function(d,v)quadratic(d=d, v=v, a=5, mini=-Inf)),
         list(name="diagonal directions", var="diagonal", val=TRUE),
         list(name="no direction change", var="dir.ch",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=5)),
         list(name="no dir. change, rel. sd", var="dir.ch.s",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1))
       )
  print.list <- c("H", "Q", "theta", "vx", "vz", "Conc", "logH")
  water.entry <-
    list(
         list(name="Printing", val="simulate", col=col.subord),
         list(name="TPrint", var="TPrint",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1000)),
         list(name="size reduction", var="red", val=paste("  ",1:7)),
         list(name="printed variable", var="print",
              val=paste(" ", print.list),
              what="none", update=TRUE), ## what is used by `simulate'
         list(name="show FE mesh", var="mesh", val=TRUE, what="none"),
         list(name="upper quantile for random field plot", var="lim.rf",
              what="none", delta=FALSE, val=function(d, v) pmin(1, pmax(0, d)),
              update=TRUE),
         list(name="quantile for swms2s plot", var="lim.swms2d",
              what="none", delta=FALSE, val=function(d, v) pmin(1, pmax(0, d))),
         ##
         list(name="Swms2d Simulation Control", val="simulate", col=col.subord),
         list(name="TolTh", var="TolTh",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=0.01)),
         list(name="TolH", var="TolH",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1)),
         list(name="lWat", var="lWat", val=TRUE),
         list(name="dt -- initial time step", var="dt",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=0.1)),
         list(name="minimal time step", var="dtMinMax[1]",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=10)),
         list(name="maximal time step", var="dtMinMax[2]",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=10)),
         list(name="max. iterations", var="max.iteration",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1000),
              what="none"),
         list(name="max. iter. (precise)", var="max.iter.prec",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=10000),
              what="none"),   
         list(name="swms message, loop period", var="iter.print",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=10000),
              what="none"),
         list(name="swms, image, loop period", var="breakpoint",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=10000),
              what="none"),
         #
         list(name="Swms2d Boundary Conditions", val="simulate", col=col.subord),
         list(name="top", var="top.bound",
              val=c("impermeable", "H constant (Dirichlet)",
                "Q constant (Neumann)",
                "variable  H", "variable Q", "atmospheric")),
         list(name="top value", var="top.value", 
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=20,
                            mini=-Inf)),
         list(name="bottom", var="bottom.bound",
              val=c("free drainage", "deep drainage (atmosph)", "impermeable",
                "H constant (Dirichlet)", "Q constant (Neumann)")),
         list(name="bottom value", var="bottom.value",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=20,
                            mini=-Inf))
         )
  
  chemical.entry <-
    list(
         list(name="Chemical Constants", val="simulate", col=col.subord),
         list(name="Solute transport", var="lChem", val=TRUE),
         list(name="scheme", var="Epsi",
              val=c("explicit", "Crank-Nicholson", "implicit")),         
         list(name="formulation", var="lUpW", val=c("upstream","Galerkin")),
         list(name="added artific. disp.", var="lArtD", val=TRUE),
         list(name="PeCr", var="PeCr",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1, maxi=10)),
         list(name="pulse duration", var="tPulse",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=1, maxi=100)),
         ##
         list(name="Chem. Swms2d Boundary Cond.", var=NULL, col=col.subord,
              val="simulate"),
         list(name="top", var="top.bound",
              val=c("impermeable", "Conc. constant (Dirichlet)",
                "Q constant (Neumann)",
                "variable  Conc", "variable Q", "atmospheric")),
         list(name="top value", var="top.value", 
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=20,
                            mini=-Inf)),
         list(name="bottom", var="bottom.bound",
              val=c("free drainage", "deep drainage", "impermeable")),
         list(name="bottom value", var="bottom.value",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=20,
                            mini=-Inf)),
         list(name="root uptake", var="root.uptake",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=20,
                            mini=-Inf)),
         list(name="internal source", var="intern.source",
              delta=TRUE, val=function(d,v) quadratic(d=d, v=v, a=20,
                            mini=-Inf))
         )
  
  print.all.list <- 
    c(lapply(as.list(print.list), function(i)
             list(name=i, var=NULL, col=col.rect,
                  val=paste("pl.water(pr.par$ps, pr.par$height,pr.par$titl,",
                    "pr.par$legend, '", i, "')", sep=""))),
      list(list(name="all smws2d results", var=NULL, col=col.rect,
                val=paste("pl.water(pr.par$ps, pr.par$height, pr.par$titl,",
                  "pr.par$legend, c('", paste(print.list, collapse="','"),"'))",
                  sep=""))))

  
  choice <- function(s, txt="", cex=1, col="blue", h.f=1.4, sp="", adj=0) {
    ## function used to draw the menue for selection the atmosphericial period,
    ## or the plant type.
    ## it works. that is fine and do not ask why.
    screen(param.dev)
    par(mar=c(0.1, 0.1, 2, 0.1))
    plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), ann=FALSE, axes=FALSE,
         xaxs="i", yaxs="i")
    title(paste("Choose the", txt, "!"), col.main=col.flash)
    w <- max(strwidth(if (is.null(sp)) s else paste(sp,s,sp), "user", cex=cex))
    h <- max(strheight(s, "user", cex=cex)) * h.f
    cols <- max(1, as.integer(1/w))
    rows <- as.integer(1/h)
    if (length(s)> cols * rows) s <- s[1:(cols * rows)]
    y <- rep(rev(1 / rows * (1:rows) - 0.5 / rows), times=cols, length=length(s))
    x <- rep(1 / cols * (1:cols) - (1 - adj) / cols, each=rows, length=length(s))
    text(x, y, cex=cex, col=col, adj=c(adj, 0.5), lab=s)
    loc <- Locator(1)
    screen(param.dev)
    plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), ann=FALSE, axes=FALSE,
         xaxs="i", yaxs="i")
    if (length(loc)>0 && (loc$x>=0) && (loc$x<1) && (loc$y>0) && (loc$y<=1) &&
        ((n<-(1+ as.integer((1-loc$y) * rows) + rows * as.integer(loc$x * cols)))
          <=length(x))) return(n)
    else return(NA)
  }

  message <- function(s) {  
    screen(message.dev)
    assign("old.message", s, envir=ENVIR) 
    plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), ann=FALSE, axes=FALSE,
         xaxs="i", yaxs="i")
    text(0.5, 0.5, adj=c(0.5, 0.5), label=s, cex=1.2, col=col.mess)
    if (s!="") cat(s,"\n")
  }
 
  linear <- function(x1, y1, x2, y2, x) 
       ((y2-y1) * x + y1 * x2 - y2 * x1) / (x2 - x1)

  quadratic <- function(d, v, a, mini=0, maxi=Inf) {
    d <- pmin(1, pmax(0, d)) - 0.5
    d <- ((d>0) * 2 - 1) * d^2 * a * 4
    if (missing(v)) d else pmax(mini, pmin(maxi, v + d))
  }
  
  redraw.horizons <- function(areas=TRUE, select=0, active=FALSE, all=TRUE) {
    cur.dev <- screen()
    screen(draw.dev, new=FALSE)
    col <- c(col.hor, rep("white", h$max.horizons))
    col[select] <- col.draw
    lim <- draw.horizons(h=h, areas=areas, col.hor=col, border.col=NULL, #"white"
                  quadratic=TRUE, picture=h$picture, all=all)
    title("Drawing area", cex.main=2,
          col.main=if (active) col.flash else col.forbid)
    screen(cur.dev, new=FALSE)
    return(lim)
  }
 
  
  rf.menu <- c("randomfield", "stone", "root", "all", "auto") # do not
  ##                                                            change order
  simulate <- function(h, first=1, PrintLevel=1, what="all",
                       wf.dev=waterflow.dev,
                       wf.label="water flux (swms_2d) ...",
                       ps=NULL,
                       ...) {
    assign("precise", NULL, envir=ENVIR) ## if simulation fails then
    ## neither a precise nor a rough simuation has been performed

    
    if (any(is.na(pmatch(what, c(rf.menu, "water"), dup=TRUE))) &&
        (length(what)!=1 || what!="none")) {
      print(pmatch(what, c(rf.menu, "water"), dup=TRUE))
      print(is.na(pmatch(what, c(rf.menu, "water"), dup=TRUE)))
      stop(paste("simulation call with strange value(s) of what:",
                 paste(what, collapse=","))) ## programming error then!
    }
    rf.choice <- pmatch(what, rf.menu, dup=TRUE)
    rf.choice <- rf.choice[!is.na(rf.choice)]
    water.simu <- is.null(h$hQThFlC) || length(rf.choice)>0 || "water" %in% what

    # print(what);  print(rf.choice); print(is.null(h$RF))
    
    if ("water" %in% what && length(rf.choice)==0) {
      ## 31.12.03 -- see at the very end of the function definition for
      ## the alternative error message "random field not simulated yet"
      rf.choice <- 5 ## rf.menu : "auto"
    }

    ## simulate the random fields first, including stones and roots
    if (length(rf.choice)>0) {
      h <- simulateHorizons(h, first=first, PrintLevel=PrintLevel,
                             message=message, what=rf.menu[rf.choice],
                             stone.trials=50)
      if (is.null(h$Root.RF)) {
        ##message("simulation failed")#message already given by simulateHorizons
        return(h)
      }
      screen(simu.dev)
      if (is.character(mess <- plotRF(h, col.txt=col.txt, col.rf=col.rf,
                                      lim=h$water$lim.rf))) {
        message(mess)
        return(h)
      }
    }

    screen(wf.dev)
    if (waterflow && h$rf.complete) {
      message(wf.label)      
      if (water.simu) {
        CW <- create.waterflow(h, MaxIt=MaxIt, hTab=hTab, DMul=DMul)        
        if (is.character(CW)) {
          message(paste("creating swms2d data:", CW))
          return(h)
        }
        
        redraw.horizons(areas)
        if (h$water$mesh) {
          screen(draw.dev, new=FALSE)
          par(new=TRUE)
          plot(Inf, Inf, xlim=draw.lim$xlim, ylim=draw.lim$ylim, axes=FALSE,
               xaxs="i", yaxs="i")
          xx <- t(cbind(CW$KXR[, c(1:4,1)],NaN))
          lines(CW$nCodeM[xx, 3], CW$nCodeM[xx, 4], col=col.mesh)
          points(CW$nCodeM[, 3], CW$nCodeM[, 4], col=col.mesh.pt, cex=0.2)
        } 
        
        h$water.x <- CW$water.x
        h$water.y <- CW$water.y
        h$flux <- CW$flux

        show.res <- function(t, swms2d){
          titl <- paste("t=", format(t, dig=3), "; ", print.list[h$water$print],
                        sep="")
          screen(wf.dev)
          if (is.character(plotWater(list(hQThFlC=swms2d, flux=h$flux,
                                          water.x=h$water.x, water.y=h$water.y),
                                     lim = h$water$lim.swms2d,
                                     what=h$water$print,
                                     col.simu=if (h$water$print==7)
                                     rev(h$col.simu) else h$col.simu,
                                     titl=titl, col.txt=col.txt, zlim=zlim,
                                     col.exception=col.exception)))
            message(paste("t=", format(t, dig=3), ": !any(is.finite(", 
                          print.list[h$water$print], "))", sep=""))
        }

        

        ask.continue <- function(t) {  
          screen(message.dev)
          plot(Inf, Inf, xlim=c(-7,2), ylim=c(0,1), ann=FALSE, axes=FALSE,
            xaxs="i", yaxs="i")
          rect(xleft=c(0.05,1.05), ybottom=rep(0.05, 2), xright=c(0.95, 1.95),
            ytop=rep(0.95, 2), col=col.bg)
          text(x=c(0.5, 1.5), y=0.5, adj=0.5, lab=c("N", "Y"), col=col.mid)
          text(-0.3, 0.5, adj=c(1, 0.5),
               label=paste("time=", format(t,dig=3)," (end=",
                           format(h$water$TPrint, dig=3),
                           ") and iter. > max.iter. Continue?",sep=""),
                           cex=1.2, col=col.mess)
          repeat {
            loc <- Locator(1)
            if (length(loc)>0 && 
               (loc$y>=0) && (loc$y<=1) && (loc$x>=0) && (loc$x<2)) {
            screen(message.dev)
            return(as.logical(as.integer(loc$x)))
            }
          }
        }
        
        result <- swms2d(CW, max.iteration=h$water$max.iteration,
                         iter.print=h$water$iter.print, message=ask.continue,
                         breakpoint=h$water$breakpoint,
                         intermediate.result=show.res)
        if (is.character(result)) {
          message(result)
          return(h)
        }
        h$hQThFlC <- result$hQThFlC 
      }

      screen(wf.dev)
      if (is.character(mess <-
                       plotWater(h, col.txt=col.txt, lim=h$water$lim.swms2d,
                                 zlim=zlim, col.exception=col.exception))) {
        h$hQThFlC <- NULL
        return(h)
      }
      assign("precise", wf.dev==prec.water.dev, envir=ENVIR) ## rather indirect
      ##                         criterion, but it works
    } else {
      if (water.simu) h$hQThFlC <- NULL
    }
    if (!h$rf.complete)
      message(if (is.null(h$RF))
              "random field not simulated yet" else
              "simulation of random fields incomplete: structure definitions missing")
    else message("")    
    return(h)
  } # simulate

                     
  get.horizon <- function() {
    message("Choose horizon or polygon!")
    i <- NA
    screen(draw.dev, new=FALSE)
    par(new=TRUE)
    lx <- length(h$grid.x)
    ly <- length(h$grid.y)
    plot(Inf, Inf, xlim=draw.lim$xlim, ylim=draw.lim$ylim, ann=FALSE, axes=FALSE,
         xaxs="i", yaxs="i")
    if (length(loc <- Locator(1))!=0) {     
      loc.x <- round((loc$x - h$grid.x[1]) / h$step)
      loc.y <- round((loc$y - h$grid.y[1]) / h$step)
      if ((loc.x>=1) && (loc.x<=lx) && (loc.y>=1) && (loc.y<=ly)) {
        i <- (h$idx.rf[loc.x, loc.y] + 1) %% h$max.horizons

        #str(h[[i]]$border)
        #str(c(loc.x, loc.y))

        
        ## 25.12.03
        ## fehler: Error in apply(abs(h[[i]]$border - c(loc.x, loc.y)), 2, sum): 
        ##     dim(X) must have a positive length
        ## 27.12.03 Error found -- error check is kept for a while...
        if (i!=1) {
          if (length(c(loc.x, loc.y))==0) {
            message("E : location unclear")
            return(NA)
          }
          if (is.null(dim(h[[i]]$border))) {
            message("E : border null")
            str(h)
            print( h$idx.rf[loc.x, loc.y] + 1)
         cat(i, loc.x, loc.y)
          return(NA)              
          }
          if (any(dim(h[[i]]$border)==0)) {
            message("E : border dim null")
            return(NA)
          }
        }
        ## end error detection part        

        if ((i==1) ||
            ## if a point on the border is choosen the choice
            ## is ignored; horizon 1 does not have border points      
            all(colSums(abs(h[[i]]$border - c(loc.x, loc.y)))!=0)) {
          redraw.horizons(areas, select=i)      
           message(if (i==1)
                   paste("Start horizon chosen                                 ")
                    else "")
        } else {
          i <- NA
          message("ambiguous choice")
        }
      } else message("out of area")
    }
    return(i)
  }

  updated <- function(x) 
    length(x$.history)==0 || x$.history[[length(x$.history)]][[2]] == "simulate"
  
  evalpar <- function(var, entry, ...) {
    eval.parameters(var=var, entry=entry,
                    h=h, update=update, simulate=simulate, dev=param.dev, 
                    col.rect=col.rect, col.bg=col.bg,
                    col.sep=col.sep, col.left=col.left,
                    col.mid=col.mid, col.right=col.right,
                    col.line = col.line, col.txt=col.txt,
                    cex=cex.eval, cex.i=cex.eval,
                    ...)
  }
 
  if (new) get(getOption("device"))(height=X11.height, width=X11.width)
  else bg.save <- par()$bg
  par(bg="white")

 
  menue.plot <- function() {
    if (par()$col==col.txt) { par(col=col.forbid); col <- col.bg.forbid}
    else { par(col=col.txt); col <- col.bg}
    cbg <- c(col.bg.forbid, col)
    for (i in (1:l.menue)[-menue.title]) {
      rect(0.1, i + 0.1, 0.9, i + 0.9, col=cbg[n.m[i]+1])
    }
    text(0.5, menue.title + 0.5, lab=menues[menue.title],
         col=col.submenue, cex=1.1)
    ctxt <- c(col.forbid, par()$col)[n.m+1]
    text(0.5, ((1:l.menue)[-menue.title]) + 0.6, adj=c(0.5, 0.5),
         menues[-menue.title], cex=1.1, col=ctxt)
  }

  menues.new <- c("horizon", "polygon")
  update.name <- c("updating: no", "updating: yes")
  update.pos <- 2
  water.name <- c("water flow: no", " water flow: yes")
  water.pos <- update.pos + 1
  prec.water.pos <- 4
  menues <- c("end",
              update.name[update+1],
              water.name[waterflow+1],
              "precise waterflow",
              "new simulation",
              "Simulation"
              ,             
              "atmosphere, control",
              "swms2d (chem)", 
              "swms2d (water)",
              "atmosphere, data",
              "root, water uptake",
              "material (chem)",
              "material (phys)",
              "Physical Parameters"
              ,
              "root growth",
              "stones",
              "structure",
              "Stochastic Parameters"
              , 
              "undo", menues.new,"Drawing",
              "postscript")
  menue.title <- c(6, 14, 18, 22)
  l.menue <- length(menues)

  scr <- NULL
  reset.screen <- function() {
    if (!is.null(scr)) close.screen(scr) 
    assign("scr",
           split.screen(figs=rbind(
                          c(0.01, 0.33, 0.01, 0.455),  # simu random field
                          c(0.34, 0.66, 0.01, 0.455),  # simu water flux
                          c(0.01, 0.33, 0.52, 0.98),  # drawing screen
                          c(0.34, 0.66, 0.51, 0.99),  # menue
                          c(0.67, 0.99, 0.01, 0.99),  # parameters (root, stones)
                          c(0.38, 0.87, 0.43, 1),     # show_models
                          c(0.01, 0.66, 0.47, 0.51),  # messages
                          c(0.67, 0.99, 0.01, 0.455),  # simu precise water flux
                          c(0.01, 0.33, 0.51, 0.99),  # invasion percolation
                          c(0.67, 0.99, 0.51, 0.99)  # few parameters (print)
                          ))
           , envir=ENVIR)    
    assign("simu.dev", scr[1], envir=ENVIR)
    assign("draw.dev", scr[3], envir=ENVIR)
    assign("showmodel.dev", scr[6], envir=ENVIR)
    assign("menue.dev", scr[4], envir=ENVIR)
    assign("param.dev", scr[5], envir=ENVIR)
    assign("waterflow.dev", scr[2], envir=ENVIR)
    assign("message.dev", scr[7], envir=ENVIR)
    assign("prec.water.dev", scr[8], envir=ENVIR)
    assign("percolation.dev", scr[9], envir=ENVIR)  ## not used yet 
    assign("few.param.dev", scr[10], envir=ENVIR)

    screen(showmodel.dev)
    par(col.main=col.txt)

    screen(message.dev)
    par(mar=c(0,0,0,0))
    if (old.message!="" && (old.loc<=4 || old.loc>=20)) message(old.message)

    screen(simu.dev)
    par(mar=c(2.2,2.2,2,0.2), col.main=col.forbid)
    plotRF(h, col.txt=col.txt, lim=h$water$lim.rf)

    screen(waterflow.dev)
    par(mar=c(2.2,2.2,2,0.2))
    plotWater(h, col.txt=col.txt, lim=h$water$lim.swms2d, zlim=zlim,
              col.exception=col.exception)

    screen(prec.water.dev)
    par(mar=c(2.2,2.2,2,0.2))
    if (!is.null(precise) && precise) {
      screen(prec.water.dev)
      plotWater(h, col.txt=col.txt, lim=h$water$lim.swms2d, zlim=zlim,
                col.exception=col.exception)
    }
    
    ## not used yet
    screen(percolation.dev)
    par(mar=c(2.2,2.2,3,0.2))

    screen(draw.dev)
    par(mar=c(2.2,2.2,2.5,0.2))
    assign("draw.lim", redraw.horizons(areas), envir=ENVIR)
    
    screen(menue.dev)
    par(mar=rep(0,4), col=col.forbid)
    plot(Inf, Inf, xlim=c(0,1), ylim=c(1, 1 + l.menue), axes=FALSE, ann=FALSE,
         xaxs="i", yaxs="i")
  }
  
  on.exit({close.screen(scr);
           if (!new) par(bg=bg.save);
           if (exists(".dev.orig")) Dev(FALSE);
           options(old.options)})
 
  old.loc <- 999
  n.m <- rep(TRUE, l.menue)
  redraw <- FALSE
  h.prec <- NULL
  old.message <- ""  
  sharpness.pos <- which(sapply(materials.entry, function(x)
                                !is.expression(x$name) && x$name=="sharpness"))
  planttype.pos <- which(sapply(uptake.entry,function(x)
                                !is.expression(x$name) && x$name=="plant type"))
  
  assign("precise", if (is.null(h$hQThFlC)) NULL else FALSE, envir=ENVIR)
  reset.screen()
  if (update) assign("h", simulate(h, 1, 0), envir=ENVIR)
  if (exists(".dev.orig")) Dev(FALSE)
  
  repeat {
    if (frequent.reset && old.loc %in% c(7:18)) reset.screen()
    n.m[(l.menue-2):(l.menue-1)] <- (h$n!=h$max.horizons)    
    n.m[prec.water.pos] <- (h$water$red!=1)
    screen(menue.dev)
    menues[update.pos] <- update.name[update + 1]
    menues[water.pos] <- water.name[waterflow + 1]
    menue.plot()
    Loc <- Locator(1)
    if (redraw) redraw.horizons(areas, all=redraw.all)
    redraw.all <- redraw <- FALSE
    message("")
    screen(menue.dev)
    menue.plot()
    if (length(Loc)==0 || Loc$x<0 || Loc$x>1) next
    loc <- floor(Loc$y)

    if ((loc>=1) && (loc<=l.menue) && (n.m[loc])) {
      switch(loc,
             {
               ## 1 : end
               break
             }, {
               ## 2 : immediate update, negation
               update <- !update
             }, {
               ## 3 :  simulate waterflux?
               if (waterflow <- !waterflow) {
                 h <- simulate(h, what="water")
                 precise <- FALSE
               } else {
                 screen(waterflow.dev)
                 plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), axes=FALSE, ann=FALSE,
                      xaxs="i", yaxs="i")
                 screen(prec.water.dev)
                 plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), axes=FALSE, ann=FALSE,
                      xaxs="i", yaxs="i")
               }             
             }, {              
               ## 4 : precise waterflow
               screen(prec.water.dev)
               plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), axes=FALSE, ann=FALSE,
                    xaxs="i", yaxs="i")
               if ((old.loc==prec.water.pos) || is.null(h$Root.RF) ||
                   is.null(h$hQThFlC)) {
                 h <- simulate(h, what="water")
               }
               h.prec <- h
               h.prec$water$red <- 1
               h.prec$water$max.iteration <-  h$water$max.iter.prec
               h.prec <-
                 simulate(h.prec,
                          what= "water", wf.dev=prec.water.dev,
                          wf.label="water flux (swms_2d) -- may take some time")
               redraw.all <- redraw <- h$water$mesh
             }, {
               ## 5 : new simulation
               h <- simulate(h)
               screen(prec.water.dev)
               plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), axes=FALSE, ann=FALSE,
                    xaxs="i", yaxs="i")
               redraw.all <- redraw <- h$water$mesh && waterflow
             }, {
               ## 6 : menue
             }, {
               ## atmosphere, control
               h <- evalpar(paste("h$atmosphere"),atmosphere.entry,
                            what="water")
               if (!update && !updated(h$atmosphere)) h$hQThFlC <- NULL
             }, {
               ## swms2d, chem
               h <- evalpar(paste("h$chem"), chemical.entry, what="water")
               if (!update && !updated(h$chem)) h$hQThFlC <- NULL
             }, {
               ## swms2d, water
               h <- evalpar(paste("h$water"), water.entry, what="water")
               if (!update && !updated(h$water)) h$hQThFlC <- NULL
             }, {
               ## atmospherical data
               col <- c(col.flash, col.mess)[1 + (diff(c(h$atmosphere$tInit,
                                                         h$atm.data[,1]))<=0)]
               j <-
                 if (nrow(h$atm.data)==1) 1
                 else choice(as.character(format(h$atm.data[,1], dig=3)),
                             txt="instance",col=col, adj=1)
               if (!is.na(j)) {
                 atm.simu <- function(atm=atm){
                   h$atm.data[j,] <- unlist(atm[atm.names])
                   assign("h", simulate(h=h, first=i, what="water"), envir=ENVIR)
                   return(atm)
                 }
                 
                 atm.cur.data <- as.list(h$atm.data[j,])
                 h$atm.data[j,] <-
                   unlist((a <- eval.parameters("atm", atm.data.entry,
                                                atm=atm.cur.data,
                                                update=update, dev=param.dev,
                                                simulate=atm.simu,
                                                col.rect=col.rect,col.bg=col.bg,
                                                col.sep=col.sep,
                                                col.left=col.left,
                                                col.mid=col.mid,
                                                col.right=col.right,
                                                col.line = col.line,
                                                col.txt=col.txt,
                                                cex=cex.eval, cex.i=cex.eval
                                     ))[atm.names])
                 ## h$rf etc is set by side effect of atm.simu
                 if (!update && !updated(a)) h$hQThFlC <- NULL
               }
            }, {
               ## water uptake
               j <-
                 if (plant.types==1) 1 
                 else choice(as.character(unlist(sapply(h$root,
                                                        function(x)x$plant.typ)
                                                 )),
                             txt="plant type", col=col.flash)
               if (!is.na(j)) {
                 h <- evalpar(paste("h$root[[",j,"]]"), uptake.entry,
                              what="root") # "water" is checked by
                 ## function 'simulate'
                 if (!update && !updated(h$root[[j]]))
                   h$hQThFlC <- h$Root.RF <- h$plants <- h$plants.idx <- NULL
               }
             }, {
               ## material definition, chem
               if (((i <- h$n) == 1) || !is.na(i <- get.horizon())) {
                 h <- evalpar(paste("h[[",i,"]]"), chem.material.entry, first=i,
                              what="water")
                 if (!update && !updated(h[[i]])) h$hQThFlC <- NULL
                 redraw <- TRUE
                 redraw.all <- areas
               }               
             }, {
               ## material definition, phys     
               if (((i <- h$n) == 1) || !is.na(i <- get.horizon())) {
                 h <- evalpar(paste("h[[",i,"]]"), materials.entry,
                              first=i, what=if (i==1) "water" else "all")
                 if (!update && !updated(h[[i]])) {
                   if (i==1 || all(lapply(h[[i]]$.history, function(x) x[[1]])
                         != sharpness.pos)) h$hQThFlC <- NULL
                   else h$hQThFlC <- h$RF <- NULL
                 }
                 redraw <- TRUE
                 redraw.all <- areas
               }
             }, {
               ## menue
             }, {
               ## roots
               j <-
                 if (plant.types==1) 1 
                 else choice(as.character(unlist(sapply(h$root,
                                                        function(x)x$plant.type)
                                                 )),
                             txt="plant type", col=col.flash)
               if (!is.na(j)) {
                 h <- evalpar(paste("h$root[[",j,"]]"), root.entry, what="root")
                 if (!update && !updated(h$root[[j]]))
                   h$hQThFlC <- h$Root.RF <- h$plants <- h$plants.idx <- NULL
               }
             }, {
               ## stone definition
               if (((i <- h$n) == 1) || !is.na(i <- get.horizon())) {
                 h <- evalpar(paste("h[[",i,"]]"), stones.entry, first=i,
                              what=c("stone", "root")) # "water" is checked by
                 ## function 'simulate'
                 if (!update && !updated(h[[i]])) h$hQThFlC <- h$Stone.RF <- NULL 
                 redraw <- TRUE
                 redraw.all <- areas
               }
             }, {
               ## structure definition of the random field
               if (((i <- h$n) == 1) || !is.na(i <- get.horizon())) {
                 screen(menue.dev)
                 screen(showmodel.dev)
                 oldmodel <-
                   if (is.null(h[[i]]$model)) model
                   else h[[i]]$model
                # zlim <-
                #   if (is.null(h$Root.RF)) NULL
                 #  else list(range(h$Root.RF, na.rm=TRUE),
                 #            range(m.link(h$Root.RF), na.rm=TRUE))
                 model <- 
                   ShowModels(x=h$grid.x, y=h$grid.y, model=oldmodel,
                              erase=FALSE, x.fraction=0.45,cex.names=cex.eval,
                              link.fct=m.link, covx.default=30,
                              update=update, Col.main="green",
                              # debug=TRUE,
                              col=col.rf, Zlim=Zlim,
                              cex.eval=cex.eval
                              )
                 model <-
                   convert.to.readable(PrepareModel(model), allowed="list")
	         if (is.null(h[[i]]$model) || 
                     !all(as.character(h[[i]]$model)==as.character(model))) {
                   h[[i]]$model <- model   
                   screen(showmodel.dev)  ## delete               
                   h$hQThFlC <- h$RF <- NULL
                   if (update) h <- simulate(h, i)
                 } else screen(showmodel.dev)  ## delete  
               }
               redraw <- TRUE
               redraw.all <- areas
             }, {
               ## menue
             }, {
               ## undo
               h$hQThFlC <- h$RF <- NULL
               if (h$n>1) {
                 ## h[[h$n]] <- NA
                 h$n <- as.integer(h$n - 1)
                 .Call("GetHorizons", h, n=as.integer(c(1,h$n)),
                       PACKAGE="SoPhy") 
                 if (update) h <- simulate(h)
               } else {
                 h[[1]]$model <- NULL
               }
               # if (!frequent.reset)
                 redraw.horizons(areas, all=TRUE)
               redraw <- FALSE
             }, {
               ## new horizon
               message("Draw boundary line of the horizon (Drawing area)!")
               screen(draw.dev, new=FALSE)
               par(new=TRUE)
               plot(Inf, Inf, xlim=draw.lim$xlim, ylim=draw.lim$ylim, axes=FALSE,
                    xaxs="i", yaxs="i")
               pts <- Locator(n=nlocator, type="l",  col=col.draw)
               n.pts <- length(pts$x)
               if (n.pts<=1) {
                 message("definition of a horizon should consist of at least 2 points")
               } else {
                 if (any(pts$x[-1] <= pts$x[1]) ||
                     any(pts$x[-n.pts] >= pts$x[n.pts])) {
                   message(" this does not seem to define a horizon!")
                   if(PrintLevel>0)
                     cat(#"this does not seem to define a horizon",
                       "(first point not the unique most left [",
                       any(pts$x[-1] <= pts$x[1]),
                       "] or last point not the unique most right [",
                       any(pts$x[-n.pts] >= pts$x[n.pts]),"])\n")
                 } else {
                   message("")
                   h$hQThFlC <- h$RF <- NULL
                   ## do not swap the following lines
                   if (pts$x[n.pts]<xlim[2]) {
                     ## do not swap the following lines
                     pts$y <-
                       c(pts$y, linear(pts$x[n.pts-1], pts$y[n.pts-1],
                                       pts$x[n.pts], pts$y[n.pts], xlim[2]))
                     pts$x <- c(pts$x, xlim[2])
                   }
                   if (pts$x[1]>=xlim[1]) {
                     ## GetHorizons considers left open right closed intervals.
                     ## So the left end point must be left from xlim[1]
                     ##
                     ## do not swap the following lines
                     pts$y <- c(linear(pts$x[1],  pts$y[1], pts$x[2], pts$y[2],
                                       xlim[1] - h$step / 2), pts$y)
                     pts$x <- c(xlim[1] - h$step / 2, pts$x)
                   }
                   lines(pts$x, pts$y, col=col.draw)
                   h$n <- as.integer(h$n + 1)
                   h[[h$n]] <-
                     list(type="H", points=pts,
                          stone=if (is.null(h[[h$n]]$stone)) stone else
                                 h[[h$n]]$stone,
                          materials=if (is.null(h[[h$n]]$materials))
                                    h[[h$n-1]]$materials else h[[h$n]]$materials,
                          model=if (!is.null(h[[h$n]]$model)) h[[h$n]]$model
                          )
                   .Call("GetHorizons", h, n=as.integer(rep(h$n, 2)),
                         PACKAGE="SoPhy")
                 }
               }
               # if (!frequent.reset)
                 redraw.horizons(areas, all=areas)
               redraw <- FALSE
            }, {
               ## new polygon
               message("Draw the boundary of a polygon (Drawing area)!")
               ##dev.set(draw.dev)
               screen(draw.dev, new=FALSE)
               plot(Inf, Inf, xlim=draw.lim$xlim, ylim=draw.lim$ylim, axes=FALSE,
                    xaxs="i", yaxs="i")
               pts <- Locator(n=nlocator, type="l", col=col.draw)
               n.pts <- length(pts$x)
               if (n.pts<=2) {
                 message("definition of a polygon should consist of at least 3 points")
               } else {           
                 message("")
                 h$hQThFlC <- h$RF <- NULL
                 ## do not swap the following lines
                 lines(pts$x[c(1, n.pts)], pts$y[c(1, n.pts)], col=col.draw)
                 pts$x <- c(pts$x, pts$x[1])
                 pts$y <- c(pts$y, pts$y[1])
                 #cat("#", h$n,": got it\n")
                 h$n <- as.integer(h$n + 1)
                 h[[h$n]] <-
                   list(type="P", points=pts,
                        stone=if (is.null(h[[h$n]]$stone)) stone else
                                h[[h$n]]$stone,
                        materials=if (is.null(h[[h$n]]$materials))
                               h[[h$n-1]]$materials else h[[h$n]]$materials,
                        model=if (!is.null(h[[h$n]]$model)) h[[h$n]]$model
                        )
                 .Call("GetHorizons", h, n=as.integer(rep(h$n, 2)),
                       PACKAGE="SoPhy")
 #               redraw.horizons(areas)
              }
               redraw.horizons(areas, all=areas)
               redraw <- FALSE
             }, {
               ## menue
             }, {
               ## printing               
               ps.mar <- function(titl) c(rep(2.2, 2), 0.3 + titl * rep(1.2, 2))
               pl.rf <- function(ps, height, titl, legend, what, transf) {
                 ps.width <- height/diff(range(h$grid.y))*diff(range(h$grid.y))
                 ps <- paste(ps, what, transf, "ps", sep=".")
                 Dev(TRUE, TRUE, ps=ps, height=height, width=ps.width)
                 par(mar=ps.mar(titl), cex=0.7)                
                 mess <- plotRF(if (!is.null(precise) && precise) h.prec else h,
                                col.txt=col.txt,
                                lim=h$water$lim.rf,
                                cex=1/0.7, cex.leg=0.8/0.7, quadratic=FALSE,
                                titl=titl, what=what, transf=transf,
                                legend=legend)
                 Dev(FALSE)
                 if (is.character(mess)) message(mess)
                 else message(paste(ps, "done."))
               }
               pl.water <- function(ps, height, titl, legend, what) {
                 ps.width <-height / diff(range(h$grid.y))*diff(range(h$grid.y))
                 for (i in what) {
                   psi <- paste(ps, i, "ps", sep=".")
                   Dev(TRUE, TRUE, ps=psi, height=height, width=ps.width)
                   par(mar=ps.mar(titl), cex=0.7)
                   mess <- plotWater(if (precise) h.prec else h,
                                     col.txt=col.txt, what=i, quadratic=FALSE,
                                     cex=1/0.7, cex.leg=0.8/0.7, titl=titl,
                                     lim=h$water$lim.swms2d, legend=legend,
                                     zlim=zlim, col.exception=col.exception)
                   Dev(FALSE)
                   if (is.character(mess)) message(mess)
                   else message(paste(psi, "done."))
                 }
               }

               val.rf <- "pl.rf(pr.par$ps, pr.par$height, pr.par$titl, pr.par$legend, what="

#               print(precise);  str(h$hQThFlC)
               
               print.par.entry <-
                 c(
                   if (!is.null(h$RF) && any(is.finite(h$RF))) 
                   list(list(name="Gauss random field", var=NULL, col=col.rect, 
                             val=paste(val.rf, "'RF',transf='none')")),
                        list(name=expression("----- " * alpha[K] * " -----"),
                             var=NULL, col=col.rect, 
                             val=paste(val.rf, "'RF',transf='K')"))
                        ),
                   if (!is.null(h$Stone.RF))
                   list(list(name="Gauss rf + stones", var=NULL, col=col.rect,
                             val=paste(val.rf, "'Stone.RF',transf='none')")),
                        list(name=expression("----- "*alpha[K]+stones*" -----"),
                             var=NULL, col=col.rect,
                             val=paste(val.rf, "'Stone.RF',transf='K')"))
                        ),
                   if (!is.null(h$Root.RF))
                   list(list(name="Gauss rf + stones + roots", var=NULL, 
                             col=col.rect,
                             val=paste(val.rf, "'Root.RF',transf='none')")),
                        list(name=expression("----- " * alpha[K] + stones + roots
                            * " -----"), var=NULL, col=col.rect,
                             val=paste(val.rf, "'Root.RF',transf='K')"))
                        ),
                   if (!is.null(precise) && !is.null(h$hQThFlC)) 
                   c(list(list(name=paste(if (precise) "precise" ,"water flow"),
                               var=NULL)), print.all.list),
                   list(list(name="plot title", var="titl", val=TRUE),
                        list(name="legend", var="legend", val=TRUE),
                        list(name="picture height [inch]", var="height",
                             delta=FALSE,
                             val=function(d, v) 1 + 9 * pmin(1, pmax(0, d))),
                        list(name="ps base name", var="ps")
                        )
                   )
               print.par <-
                 eval.parameters("pr.par", print.par.entry, pr.par =print.par,
                                 pl.water=pl.water, pl.rf = pl.rf,
                                 update=FALSE,
                                 dev=if (length(print.par.entry)<10)
                                     few.param.dev else param.dev,
                                 col.rect=col.rect, col.bg=col.bg,
                                 col.sep=col.sep, col.left=col.left,
                                 col.mid=col.mid, col.right=col.right,
                                 col.line = col.line, col.txt=col.txt,
                                 cex=cex.eval, cex.i=cex.eval)
             }
             )
      old.loc <- loc
    }
  }

  invisible(h)
}


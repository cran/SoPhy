## code used in the papers by Schlather and Huwe
################################################
# library(SoPhy); RFparameters(Print=7); source("~/R/SOPHY/SoPhy/R/sh.R"); source("~/R/SOPHY/SoPhy/R/3dplot.R");sh.jch(11, bw=TRUE, readlines=FALSE, dev=TRUE, final=TRUE, Print=0)

sh.jch <- function(input=NULL, dev=2, pspath="./", txt.result.dir = "txt/",
                   simu.path="simu/", standard.method = "fix.m",
                   final=TRUE, PrintLevel=0, bw=TRUE, readlines=TRUE
                   ) {
  ## additional, last input value = Inf to exit safely after simulation
  ## basic definitions, miscellaneous

  keep.figures <- FALSE
  ## keep.figures <- TRUE ## files *.m and *.l are much larger, but contain
  ##                         full information
  simu.pspath <- "./ps"
  if (final) {
    simu.dir <- "./final.simu/"
    tinylambda <- 0.01
    largelambda <- 0.6
    simu.repet <- 100
   # all.methods <- c("fix.m", "optim.m") ## used in all.estimated
    all.methods <- c("fix.m", "optim.m", "ml") ## used in all.estimated
    all.measures <- c("robust", "lsq")   ## dito
    sens.nlines <- 500
    sens.measure <- "robust"
    sens.repet <- 250
    sens.boot.repet <- 250
    delta <- NA  ## automatic determination of area around considered slice
    ##              to avoid side effects
    figures <- c("F04", "F06", "K06")
    data.selected.dist <- 0.75
    data.ppi.par <- 5
    data.measures <- c("robust", "lsq")
    type.select <- c("identical", "unif", "independent", "largelambda",
                    "dependent", "all")
    model.param.select <- c("xiNN", "xiN", "xi0", "xiP", "xiPP")
    select.dist <- c(0, 0.95)
    lambda.F04 <- 1.5
    lambda.F06 <- 1.0
    lambda.K06 <- 1.0
  } else {
    simu.dir <- "./prelim.simu/"
    tinylambda <- 0.01
    largelambda <- 0.01
    simu.repet <- 2
    all.methods <- c(standard.method)   ## used in all.estimated
    all.measures <- c("lsq") ## dito
    sens.nlines <- 10
    sens.measure <- "lsq"
    sens.repet <- 2
    sens.boot.repet <- 2
    delta <- 3
    figures <- c("F04")
    # figures <- c("F04", "F06", "K06")
    data.selected.dist <- c(0.1, 0.15)
    data.ppi.par <- 1
    data.measures <- c("lsq")
    type.select <- c("identical", "independent")
    model.param.select <- c("xiP", "xiPP")
    select.dist <- c(0, 0.05) 
    lambda.F04 <- 0.3
    lambda.F06 <- 0.2
    lambda.K06 <- 0.2
}

  
  
  ## sensitivity analysis as described in the last paragraph of section 4
  sens.simu.name <- "sensitivity.analysis.simu.data"
  sens.bundle.name <- "sensitivity.analysis.bundle.data"
  sens.result.name <- "bundle.smsd"
  sens.psbasename <- "simu2"
  sens.method <- "fix.m"
  sens.bundle <- 1 # sens.bundle <- c(1, 2, 4, 8, 16, 32)
  
  sens.type <- list(F04 =
                    list(type="all",
                         delta.x=40,
                         delta.y=40, 
                         endpoint.tolerance=-50,
                         length= 547 - 56 + 1,
                         lambda=lambda.F04,
                         x.var=  0.16,
                         depth=  (126 - 52),
                         drops=1,
                         xi=-0.72,
                         drop.scale=22,
                         x.h.scale=20, x.v.scale=5, 
                         dd=function(x) (1 - x^(-tt$xi)) * 50
                         ),
                    F06 =
                    list(length=  579 - 51 + 1,
                         lambda= lambda.F06,
                         depth = (407 - 33) ,
                         xi=0.25, 
                         drop.scale=20,
                         x.h.scale=25, #20, 3.9.04
                         x.v.scale=10, #20, 3.9.04
                         dd=function(x) (x^(-tt$xi) - 1) * 250/(0.05^(-tt$xi)-1)
                         ),
                    K06 =
                    list(length= 468  - 118 + 1,
                         lambda=lambda.K06,
                         depth = (350 - 29),
                         x.var= 0.2, #0.1, 3.9.04
                         xi=0.38, 
                         drop.scale=5,
                         x.h.scale=6, x.v.scale=5,
                         dd=function(x) (x^(-tt$xi) - 1)*(150)/(0.05^(-tt$xi)-1)
                         )
                    )
  stopifnot("F04" %in% figures)
  sens.type <- sens.type[figures]
  sens.linesimustep <- 0.1

  ## simulated figures
  pict.name <- "simu"
  pict.several.name <- "real"
  pict.repet <- 6

  ## data
  fig.path <- paste(system.file(package='SoPhy'), 'tracer', sep="/")
  ## fig.path <- "/home/schlather/article/invasion/Pictures/"
  data.ana.ext <- ".l"
  data.estim.ext <- ".m"
  data.dir <- "data/"
  data.par.name <- "param.txt"
  data.risk.name <- "daten.xi.histo.txt"
  data.plot.all <- TRUE
  data.endpoint.tolerance <- -50			       
  data.method <- "fix.m"
  data.mar <- c(4.1, 4.8, 0.3, 0.4)


  ## simulation study
  simu.name <- "parameters"
  result.name <- "smsd"
  delta.name <- "delta.dat"
  ps.name <- "study" ## evaluation of simulation
  med.min <- 0.5
  med.max <- 0.8
  length.for.largelambda <- 500
  depth <- 100
  method <- standard.method
  measure <- "robust"
  simu.eval.variables <- 1 ## subset of c(1,2); 1:xi, 2:sigma
  #continue <- TRUE   ### continue simulation if stopped -- do not restart
  upper.bound <- 0.1 ### ??
  type <- list(identical=
               list(type="identical", delta.y=1, delta.x=1,
                    lambda=tinylambda,
                    length=length.for.largelambda * largelambda / tinylambda,
                    x.h.scale=20, x.v.scale=5, x.var=0,
                    endpoint.tolerance=-50
                    ),
               unif =
               list(type="unif", unif.b=2, delta.y=1, delta.x=2
                    ),
               independent=
               list(type="independent", x.var=0.16, delta.y=delta, delta.x=delta
                    ),
               largelambda=
               list(lambda=largelambda, length=length.for.largelambda
                    ),
               dependent=
               list(type="dependent"
                    ),
               all =
               list(type="all"
                    )
               )
  type <- type[type.select]
  model.param <- list(xiNN=list(xi= -1.5, drops=1, drop.scale=10,
                           dd=function(x) (1 - x^(-m$xi)) * 100
                           ),
                      xiN=list(xi=-1/2),
                      xi0=list(xi=0, 
                           dd=function(x) -log(1 - x) * 20 # 50, 20      
                           ),
                      xiP=list(xi=1/2, drop.scale=5, drops=1,
                           dd=function(x) (x^(-m$xi) - 1) * 100/(0.05^(-m$xi)-1)
                           ),
                      xiPP=list(xi=1.5) #              0.1, 0.005
                      )
  model.param <- model.param[model.param.select]
  hist.bw <- bw
  upper.pos.xi <- 9.5 ## estimates of xi beyond this value are very
  ##                     likely due to numerical or statistical errors
  p.sideeffect <- 0.999  ## 1-p = probability for delta being too small
  ##                        value used in the automatic determination of delta
  ##                        in case delta = NA
  linesimustep <- 0.1

  ## standard modifications of names
  for (pathname in c("pspath", "txt.result.dir", "simu.path")) {
    path <- get(pathname)
    if (substr(path, nchar(path), nchar(path)) != "/")
      assign(pathname, paste(path,  "/", sep=""))
  }
  simu.dir <- paste(simu.path, simu.dir, sep="")
  simu.pspath <- paste(simu.path, simu.pspath, sep="")
  data.dir <- paste(simu.path, data.dir, sep="")
  
  for (name in c(simu.path, simu.dir, simu.pspath, data.dir,
                 pspath, txt.result.dir)) {
    if (PrintLevel>2) print(name)
    if (file.exists(name)) stopifnot(file.info(name)$isdir)
    else stopifnot(dir.create(name))
  }

  for (name in c("sens.simu.name", "sens.bundle.name", "simu.name"))
    assign(name, paste(simu.dir, get(name), sep=""))
  for (name in c("sens.result.name", "result.name",
                 "data.par.name", "data.risk.name"))
    assign(name, paste(txt.result.dir, get(name), sep=""))

  ## further settings
  if (!exists(".Random.seed", envir=.GlobalEnv)) runif(1)
  seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
  close.screen(close.screen())
  Measures <- list(robust=function(x) abs(x), lsq=function(x) x^2)
  Methods <- c("fix.m", "optim.m", "ml") #2: m not optimised, 4:m optimised
  ParIdx <- function(method, measure) {
    which(method==Methods) +
      (which(measure==names(Measures)) - 1) * length(Methods)
  }
  
  selected.rate <- c(med.max, med.min)
  if (exists(".dev.orig")) Dev(FALSE)
  if (is.logical(dev)) {
    par(cex.axis=1.8, cex.lab=1.8)
  } else {
    par(cex.axis=1, cex.lab=1)
  }
  
  if (is.na(delta)) {
    RFparameters(Print=0, pch="")
    ## wieviel raemlicher Puffer wird benoetigt, damit keine
    ## Randeffekte mehr auftreten??
    tt.cmp <- NULL
    filename <- paste(simu.dir, delta.name, sep="")
    if (file.exists(filename)) load(filename)
    else cat("Prelimary calculations. This will take some seconds...\n")
    tt <- list()
    for (typ in 1:length(type)) {
      for (nt in names(type[[typ]])) {
        txt <- paste("tt$", nt, "<-", "type[[typ]]$", nt)
        eval(parse(text=txt))
      }
      if (tt$type!="independent") next
      if (length(tt.cmp) != length(tt) || !all(unlist(tt.cmp)==unlist(tt),
                  na.rm=TRUE) || sum(is.na(unlist(tt.cmp)==unlist(tt)))!=2) { 
        fp <- flowpattern(type=tt$type,
                          x.h.scale=NA, x.v.scale=tt$x.v.scale, x.var=tt$x.var,
                          length=100, lambda=1, width=100, delta.x=0, delta.y=0,
                          method=NULL, simu.method=NULL,  drop.simu.method=NULL,
                          PrintLevel=0)$intermediate[c("xx", "yy")]
        fp$xx <- abs(fp$xx - fp$xx[, 1])
        fp$yy <- abs(fp$yy - fp$yy[, 1])
        delta <- max(quantile(fp$xx, p.sideeffect),
                     quantile(fp$yy, p.sideeffect))
        fp <- NULL
        tt.cmp <- tt
        save(file=filename, tt.cmp, delta)
        tt.cmp <- NULL
        }
      type[[typ]]$delta.y <- type[[typ]]$delta.x <- delta
      break;
    }
  }

 
  rl <- if (readlines) # || .Platform$OS.type!="unix")
    function(x) readline(if (x=="")"press return" else paste(x,": press return"))
  else function(x) { cat(x); sleep.milli(1000); cat("\n")}

  sketches <- function(dev, bw) {
    if (bw) {
      col.1 <- col.2 <- col.3 <- col <- "black"
      lty <- 1:3
    } else {
      col <- c("red", "blue", "green")
      lty <- 1
      col.1 <- "brown"; col.2 <- "blue"; col.3 <- "green" 
    }
    ## paths (Figure 1)
    Pn <- 1 
    x <- seq(0, 20, 0.01)
    enlarge <- 5
    cex <- 2.3
    P <- GaussRF(x=c(range(x), x[2]-x[1]), model="whittle",
               param=c(0, 10, 0, 3.2, 2.5), grid=TRUE, gridtriple=TRUE, n=Pn)
    P <- as.matrix(P)
    Plim <- range(P)
    for (i in 1:Pn) {
      for (n in 0:3) {
        Dev(TRUE, dev, ps=paste(pspath, "Pfad-", i,"-", n, ".eps", sep=""),
            height=5, width=2.5)
        par(mar=c(3.7, 4.5, 0.3, 0.3), cex=1.0)
        farbe <- if (n==0) "black" else "white"
        farbe <- "black"
        plot(Inf, Inf, xlim=Plim, ylim=-c(max(x), min(x)),
             col.axis=col, col.lab=col,
             xlab="x", ylab="depth", cex.axis=cex, cex.lab=cex, yaxs="i")
        par(new=TRUE)
        plot(P[, i], x, type="l", cex.axis=cex, cex.lab=cex,
             ylim=c(max(x), min(x)), lwd=5, col.axis=farbe, col.lab=farbe,
             xlab="x", ylab="depth", yaxs="i", xlim=Plim, frame=TRUE, axes=FALSE)
        if (n>0) {
          Z <- runif(1, min(x), max(x))
          lines(P[x<Z, i], x[x<Z], col="lightblue", lwd=15)
        }
        Dev(FALSE, dev)
        if (is.numeric(dev)) rl(paste("path", n))
      }
    }

    ## pareto density
    x <- seq(0.00001, 4, 0.01)
    height <- 3.5
    width <- 12
    Dev(TRUE, dev, ps=paste(pspath, "paretodensity", sep=""),
        height=height, width=width)
    par(mar=c(4.5, 4.5, 0.2, 0.2), cex=1)
    y <- cbind(pmin(1.2, dpareto(x, -1/2, 1)),dpareto(x, 0, 1),dpareto(x, 1, 1))
    matplot(x, y, type="l", lty=lty, lwd=5, cex.axis=2, cex.lab=2,
            ylab="", xlab="depth", col=col)
    legend(2.8, 1.015, legend=c("h(z ; -0.5, 1)", "h(z ; 0, 1)", "h(z ; 1, 1)"),
           lty=lty, cex=2, lwd=5, col=col)
    Dev(FALSE, dev)
    if (is.numeric(dev)) rl("pareto densities")

    
   ## convolution with uniform distribution 
    faltung <- function(q, xi, s=1, lower.tail=TRUE, a, b) {
      stopifnot(xi!=0)
      res <- (pmax(0, 1 + xi * a * x / s)^(1 - 1/xi) -
              pmax(0, 1 + xi * b * x / s)^(1 - 1/xi)) * s / ((1-xi) * x *(b - a))
      if (lower.tail) res <- 1 - res
      res
    }
    
    x <- seq(0, 1, len=50)
    xi <- -1.5
    a <- 1
    b <- 2
    
    Dev(TRUE, dev, ps=paste(pspath, "uniform", sep=""), height=3.7, width=3.7)
    par(mar=c(4, 4.5, 0.1, 0.1))
    plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), xlab="z", ylab="1-H^*(z)", 
         cex.axis=1.5, cex.lab=1.5)
    lines(x, faltung(x, xi, s=abs(xi), low=FALSE, a=a, b=b), lwd=3)
    lines(x, ppareto(x, xi, s=2/3*abs(xi), low=FALSE), lty=2, lwd=3)
    lines(x, 2 / 3 * ppareto(x, xi/(1-xi), s=abs(xi/(1-xi)), low=FALSE),
          lty=3,lwd=3)
    Dev(FALSE, dev)
    if (is.numeric(dev))
      rl("convolution for variables angles")

    ## block sketch
    Dev(TRUE, dev, ps=paste(pspath, "slice.sketch", sep=""), height=5, width=10)
    par(mar=c(4, 4.5, 2.2, 0.1), cex=1.2)
    laenge <- 200
    breite <- 5 # 3
    tiefe <- 100
    deltax <- 50
    deltay <- 50
    startx <- 5
    starty <- 3
    quader(size=c(laenge + 2* deltax, breite + 2 * deltay, tiefe), #
           bottomleftfront=c(-deltax, -deltay, 0),
           inf=c(1300, 250, 1/800), sun=c(laenge/2, -1500 - deltay, 0), 
           col=grey(c(rep(1, 200), seq(1, 0, -0.001), rep(0, 1000))),
           unit="pixels", cex.axis=2.2, srt=8)
    quader(size=c(laenge, breite, tiefe), bottom=c(0, 0, 0), add=TRUE,
           lty=c(2,2), col.frame=rep("grey20", 2))
    n <- 3
    param <- c(0, sqrt(laenge + 2* deltax) * 4, 0, 5, 2)
    x <- t(c(-startx, laenge / 2, laenge + startx) + 
           t(GaussRF(x=1:100, grid=TRUE, model="whittle", param=param, n=n)))
    
    y <- t(c(-starty, breite/2, breite + starty) + 
           t(GaussRF(x=1:100, grid=TRUE, model="whittle", param=param, n=n)))
    xx <- rbind(as.vector(x), as.vector(y), rep(1:100, 3))
    pp <- pos3D(xx)
    xi <- matrix(pp[1,], ncol=3)
    yi <- matrix(pp[2,], ncol=3)
    p3d.sun <- get(".p3d.sun")
    p3d.col <- get(".p3d.col")
    p3d.colfct <- get(".p3d.colfct")
    for (i in 1:n) {
      for (j in 2:(nrow(xi)-1)) {
        lines(xi[(j-1):j, i],  yi[(j-1):j, i],
              lwd=if ((y[j,i]>=0) & (y[j,i]<= breite) &
                (x[j,i]>=0) & (x[j,i]<=laenge)) 
              8 else 3 ,
              col = p3d.col[p3d.colfct(sqrt(sum((c(x[j,i], y[j,i], j) 
                - p3d.sun)^2)))]
              )
      }
    }
    Dev(FALSE)
  } ## function sketches

  sensitivity.simu.print <- function(dev, bw) {
    risk <- NULL
    load(sens.bundle.name) ## risk is loaded here
    tt <- list()
    if (!is.numeric(dev)) write(file=sens.result.name, "")
    for (typ in 1:length(type)) {
      for (nt in names(type[[typ]])) 
        eval(parse(text=paste("tt$", nt, "<-", "type[[typ]]$", nt)))
      var.value <- tt$xi
      for (b in 1:length(sens.bundle)) {
        med <- risk[[typ]][b, ]
        smsd <- format(sqrt(mean((med - var.value)^2, na.rm=TRUE)), dig=2)
        sd <- sqrt(var(med, na.rm=TRUE))
        max.sd <- sqrt(sens.bundle[1] / sens.bundle[b]) *
          sqrt(var(risk[[typ]][1, ], na.rm=TRUE))
        txt <- paste(formatC(LETTERS[typ], width=1),
                     "   xi=", formatC(var.value, dig=3, width=5),
                     "   n=", formatC(sens.bundle[b], width=2),
                     "   mean=",formatC(mean(med, na.rm=TRUE),
                                        flag="-", dig=3, width=6),
                                        # "   smsd=",formatC(smsd, width=5),
                     "   sd=", formatC(sd, flag="-", dig=3, width=6),
                     "   max.sd=", formatC(max.sd, flag="-", dig=3, width=6),
                     "   %sd=", formatC(100 * (sd-max.sd)/max.sd, flag="-",
                                        dig=3, width=6),
                     sep="")
        if (is.numeric(dev)) cat(txt, "\n")
        else write(file=sens.result.name, txt, append=TRUE)
      }
    } # typ
  }
  
  sensitivity.simu.eval <- function(bootstrap.repet) {
    load(sens.simu.name) #  frequencies, repet, type    
    if (file.exists(sens.bundle.name)) {
      load(sens.bundle.name)
    } else {
      risk <- list()
      for (typ in 1:length(type)) {
        risk[[typ]] <- matrix(nrow=length(sens.bundle), ncol=bootstrap.repet)
      }
    }    
    tt <- list()
    for (typ in 1:length(type)) {
      for (nt in names(type[[typ]])) 
        eval(parse(text=paste("tt$", nt, "<-", "type[[typ]]$", nt)))
      if (PrintLevel>1) cat("type", typ, "\n")
      frequencies[[typ]] <-
        frequencies[[typ]][, 1:sum(colSums(is.na(frequencies[[typ]])) ==0)]
      #####        passt obiges so ??
      maxdist <- nrow(frequencies[[typ]])
      nn <- ncol(frequencies[[typ]])
      dist <- 1:maxdist
      for (repet in 1:bootstrap.repet) {
        cat(formatC(repet, width=4),"")
        for (b in 1:length(sens.bundle)) {
          idx <- as.integer(runif(sens.bundle[b], 1, nn + 1))
          freq <- rowSums(frequencies[[typ]][, idx, drop=FALSE])
          ri <- risk.index(cbind(dist, freq),
                           selected.dist=NULL,
                           PrintLevel=PrintLevel,
                           endpoint.tolerance=tt$endpoint.tolerance,
                           ## weights = dist^2,
                           measure=Measures[sens.measure][[1]],
                           method=sens.method)$risk.index
          risk[[typ]][b, repet] <- ri
          if (PrintLevel>2) cat(formatC(ri, width=6, dig=2), "")
        }
        if (PrintLevel>2) cat("\n")
        if (PrintLevel>0) cat("saving...")
        save(file=sens.bundle.name,  risk, sens.bundle, type)
        if (PrintLevel>0) cat("\n")
      }
    }
  }
  
  sensitivity.simu <- function(final, repet, dev, bw) {
    RFPrintLevel <- PrintLevel
    type <- sens.type
    nn.start <- 1
    
    if (repet>1) {
      cat("R can be interrupted whenever data are not being saved\n",
          "and the simulation can be continued later on. \n",
          "Note that the number of repetetions is memorised.\n",
          "The evaluation algorithm runs on the simulations done up to \n",
          "now if not all the simulations have been performed.\n"
          )
      if (!(start <- !file.exists(sens.simu.name))) {
        load(sens.simu.name)
        assign(".Random.seed", seed, envir=.GlobalEnv)
      }
    } else {
      start <- TRUE
    }
    RFparameters(Print=RFPrintLevel, pch="", # TBM3.every=50,
                 TBM3.lines=sens.nlines,
                 TBM3.linesimufactor = 0.0,
                 TBM2.linesimufactor = 0.0,
                 TBM3.linesimus=sens.linesimustep,
                 TBM2.linesimus=sens.linesimustep)
    split.screen(c(2,2))
    screen(1, new=FALSE)
    tt <- list()
    
    zeit <- system.time(if (nn.start<=repet) for (nn in nn.start:repet) {
      if (PrintLevel>1) cat(nn, "")
      for (typ in 1:length(type)) {
        for (nt in names(type[[typ]])) 
          eval(parse(text=paste("tt$", nt, "<-", "type[[typ]]$", nt)))
        if (PrintLevel>1) cat("\n   ####  type=", typ, " (", tt$type,")", 
                             " (xi=", tt$xi, ")  ####\n", sep="")
        environment(tt$dd) <- environment()        
        estim <- 
          flowpattern(type=tt$type,
                      x.h.scale=tt$x.h.scale, x.v.scale=tt$x.v.scale,
                      x.var=tt$x.var,
                      length=tt$length, depth=as.integer(tt$depth),
                      lambda=tt$lambda, delta.y=tt$delta.y,
                      delta.x=tt$delta.x,
                      endpoint.tolerance=tt$endpoint.tolerance,  
                      drop.distr=tt$dd, drops=tt$drops,
                      drop.sc=tt$drop.scale,
                      select=NULL,
                      method= NULL,
                      measure=NULL,
                      PrintLevel=PrintLevel)
        if (start) {
          start <- FALSE
          frequencies <- list()
          frequencies[[length(type) + 1]] <- NA
        }
        if (is.null(frequencies[[typ]])) {
          frequencies[[typ]] <- array(dim=c(length(estim$freq), repet))
        }
        frequencies[[typ]][, nn] <- estim$freq
        if (repet==1) {
          estim <- cbind(estim$i.x, estim$i.d)
          save(file=paste(sens.simu.name, typ, sep="."), estim)
        }
        gc()
      } #  typ      
      if (repet>1) {
        if (PrintLevel>0) cat("saving...")
        seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
        nn.start <- nn + 1
        save(file=sens.simu.name, frequencies, repet, type, seed, nn.start)
        if (PrintLevel>0) cat("\n")
      }
      gc()
    } # repet
                )[1] # system.time 
    close.screen(close.screen())

    ### plotting the figures
    if (repet==1) {
      cat("\nSingle simulation takes about", zeit,
          "seconds -- remember when\n performing the simulation study \n\n")
      tt <- list()
      for (typ in 1:length(type)) {
        for (nt in names(type[[typ]])) 
          eval(parse(text=paste("tt$", nt, "<-", "type[[typ]]$", nt)))
        load(paste(sens.simu.name, typ, sep="."))
        psname <- paste(pspath, sens.psbasename, ".", typ,  sep="")
        Dev(TRUE, dev, ps=psname, height=tt$depth * 0.015,
            width=tt$length * 0.015)
        par(mar=c(4.4, 4.8, 0.1, 0.3), bg="white")
        plot(estim[, 1], -estim[, 2],
             ylim=c(-tt$depth, -1),
             pch=".", xlab="", ylab="", main="", 
             col=if (bw) "black" else "blue",
             cex=if (is.logical(dev)) 8 else 1,
             cex.axis =if (is.logical(dev)) 1.8 else 1,
             cex.lab = if (is.logical(dev)) 1.8 else 1
             )
        Dev(FALSE)
        if (is.numeric(dev) && typ !=length(type)) rl(paste("type", typ,tt$type))
      }
    }
  } # sensitivity.simu

  simulation <- function(repet,
                         indices=if (i.start<=repet) i.start:repet else NULL) {
    ## indices allows the simulation of certain slots
    split.screen(c(2,2))
    screen(1, new=FALSE)
    xi <- numeric(length(model.param))
    start <- TRUE
    tt <- list()    
    par.idx <- ParIdx(method=method, measure=measure)
    i.start <- 1
    if (file.exists(simu.name)) {
      load(simu.name)
      assign(".Random.seed", seed, envir=.GlobalEnv)
    }
    if (!is.null(indices)) for (i in indices) {
      if (PrintLevel>1) cat(i, "")
      simu.storage <- list()
      for (typ in 1:length(type)) {
        for (nt in names(type[[typ]])) {
          txt <- paste("tt$", nt, "<-", "type[[typ]]$", nt)
          eval(parse(text=txt))
        }
        m <-  list()
        old.simu <- NULL
        for(m.p in 1:length(model.param)) {
          for (n in names(model.param[[m.p]])) {
            txt <- paste("m$", n, "<-", "model.param[[m.p]]$", n) 
            eval(parse(text=txt))
          }
          xi[m.p] <- m$xi
          if (PrintLevel>1)
            cat("\n **** run=", i, "   type=", typ, " (", tt$type,")", 
                "   model=", m.p, " (xi=", m$xi, ") **** \n", sep="")

          environment(m$dd) <- environment()
          estim <-
            flowpattern(type=tt$type,
                        depth=depth,
                        x.h.scale=tt$x.h.scale, x.v.scale=tt$x.v.scale,
                        x.var=tt$x.var,
                        length=tt$length, lambda=tt$lambda,
                        delta.x=tt$delta.x, delta.y=tt$delta.y,
                        endpoint.tolerance=tt$endpoint.tolerance, 
                        unif.b=tt$unif.b,  
                        drop.distr=m$dd, drops=m$drops, drop.sc= m$drop.scale,
                        method=method, select=select.dist,
                        measure=Measures[measure][[1]],
                        old=old.simu,
                        PrintLevel=PrintLevel)
          if (is.null(old.simu)) 
            simu.storage[[typ]] <- old.simu <- estim$intermediate
          if (start) {
            start <- FALSE
            parameters <- list()            
            parameters[[par.idx]] <-
              array(dim=c(dim(estim[[1]]$par), repet,
                      length(type), length(model.param)))
            frequencies <-
              array(dim=c(length(estim$freq), repet, length(type), 
                      length(model.param)))
          }
          frequencies[, i, typ, m.p] <- estim$freq
          parameters[[par.idx]][, , i, typ, m.p] <- estim[[1]]$par
        } ## m.p
      } #  typ      
      seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
      i.start <- i + 1
      sel.dist <-  estim[[1]]$sel.dist
      save(file=simu.name, parameters, frequencies, xi, model.param, sel.dist,
           type, model.param, Methods, Measures, ParIdx, repet, seed, i.start)
      simu.storage <- NULL
    } # i
    close.screen(close.screen())
  }
  
  all.estimate <- function() {
    frequencies <- repet <- NULL
    load(simu.name) ## parameters, frequencies, xi, model.param, 
    ##                 type, model.param, Methods, Measures, ParIdx
    split.screen(c(2,2))
    screen(1, new=FALSE)
    tt <- list()
    for (i in 1:repet) {
      for (typ in 1:length(type)) {
        for (nt in names(type[[typ]])) {
          txt <- paste("tt$", nt, "<-", "type[[typ]]$", nt)
          eval(parse(text=txt))
        }
        m <-  list()
        for(m.p in 1:length(model.param)) {
          for (n in names(model.param[[m.p]])) {
            txt <- paste("m$", n, "<-", "model.param[[m.p]]$", n) 
            eval(parse(text=txt))
          }
          if (PrintLevel>1)
            cat("\n **** run=", i, "   type=", typ, " (", tt$type,")", 
                "   model=", m.p, " (xi=", m$xi, ") **** \n", sep="")
          for (met in all.methods) {
            for (mea in all.measures) {
              par.idx <- ParIdx(method=met, measure=mea)
              if (length(parameters)>=par.idx &&
                  !is.null(parameters[[par.idx]]) &&
                  any(!is.na(parameters[[par.idx]][, , i, typ, m.p]))) {
                if (PrintLevel>1) cat(met, mea, " jumped\n")
                next
              }
              estim <- risk.index(cbind(1:dim(frequencies)[1],
                                        frequencies[, i, typ, m.p]),
                                  selected.dist=select.dist,
                                  selected.rate=NULL,
                                  endpoint.tolerance=tt$endpoint.tolerance,
                                  method=met, 
                                  measure=Measures[mea][[1]],
                                  PrintLevel=PrintLevel)
              if (length(parameters)<par.idx || is.null(parameters[[par.idx]]))
                parameters[[par.idx]] <-
                  array(dim=c(dim(estim$par), repet, length(type),
                          length(model.param)))
              parameters[[par.idx]][, , i, typ, m.p] <- estim$par
            } # mea
          } # met
        } # m.p
        sel.dist <- estim$sel.dist
        save(file=simu.name, parameters, frequencies, m$xi,
             model.param, sel.dist,
             type, model.param, Methods, Measures, ParIdx)
      } # typ
    } # i / repet
    close.screen(close.screen())
  }
  
  eval.simulation <- function(dev) {
    frequencies <- sel.dist <- parameters <- NULL
    load(simu.name)
    if (is.numeric(dev)) split.screen(c(2,2))
    max.freq <- apply(frequencies[, , , ], c(2, 3, 4), max, na.rm=TRUE)
    if (any(idx <- !is.finite(max.freq))) {
      print(max.freq)
      stop("non-finite values in max.freq -- please contact author")
    }
    idx <- sel.dist
    ## 1:(select.dist[2] * depth) ## nicht nachgeprueft, ob definition ok
    ##                                   fuer depth != 100
    distances <- rep(idx, dim(frequencies)[2]) 
    # distances <- rep(1:(select.dist[2] * depth), dim(frequencies)[2]) 
    if (!is.numeric(dev)) write(file=result.name, "")
    tt <- list()
    for (typ in 1:length(type)) {
      for (nt in names(type[[typ]])) 
        eval(parse(text=paste("tt$", nt, "<-", "type[[typ]]$", nt)))
      m <-  list()
      for (m.p in 1:length(model.param)) {
        for (n in names(model.param[[m.p]])) 
          eval(parse(text=paste("m$", n, "<-", "model.param[[m.p]]$", n)))
        if (PrintLevel>1) cat("model=", m.p, " type=", typ, "\n", sep="")
        name.ps <- paste(paste(simu.pspath, ps.name, sep=""),
                         tt$type, tt$lambda * 100,
                         round(tt$x.var * 100), m$xi, sep=".")
        Dev(TRUE, dev, ps=paste(name.ps, "freqhist", sep="."), heigh=5, widt=5)
        par(cex=1, mar=c(2.1, 2.2, 2.1, 0.1)) 
        hist(max.freq[, typ, m.p], 20,
             main=paste(quantile(max.freq[, typ, m.p], na.rm=TRUE),
               collapse=", "))
        Dev(FALSE)
        if (is.numeric(dev)) rl(paste(tt$type, "xi=", m$xi, " histo"))
        for (method in Methods) for (measure in names(Measures)) {
          ii <- ParIdx(method=method, measure=measure)
          if (length(parameters) < ii || is.null(parameters[[ii]])) {
            cat(method, "&", measure, "are not available\n")
            next
          }
          for (variable in simu.eval.variables) {
            if (variable==1) {
              var.value <- m$xi
              var.factor <- 2
              var.upper.allowed <- upper.pos.xi
              var.name <- "xi"
            } else {
              var.value <-
                if (m$xi < 0) -m$xi * 100 else
              if (m$xi==0) 20 else
              if (m$xi > 0) -m$xi * 100 / (0.05^(m$xi) -1) else NA
              var.factor <- 5
              var.upper.allowed <- 100
              var.name <- "s"
            }
            ## evaluate for all repetitions simultaneously
            x.coord <- t(t(frequencies[idx, , typ, m.p]) / max.freq[, typ, m.p])
            xlim <- c(1, 0)
            y.coord <- parameters[[ii]][variable, , , typ, m.p]
            med <- y.coord
            med[x.coord>med.max | x.coord<med.min] <- NA
            med <- apply(med, 2 , median, na.rm=TRUE)
            if (any(is.na(med))) {
              med2 <- y.coord
              med2[x.coord>med.max + 0.1 | x.coord<med.min - 0.1] <- NA
              med2 <- apply(med2, 2 , median, na.rm=TRUE)
              med[is.na(med)] <- med2[is.na(med)]
              
              if (any(is.na(med))) {
                med2 <- y.coord
                med2[x.coord>med.max + 0.2 | x.coord<med.min - 0.2] <- NA
                med2 <- apply(med2, 2 , median, na.rm=TRUE)
                med[is.na(med)] <- med2[is.na(med)]
              }
              if (!final && any(is.na(med))) {
                warning("bad approximation of med")
                med2 <- y.coord
                med2 <- apply(med2, 2 , median, na.rm=TRUE)
                med[is.na(med)] <- med2[is.na(med)]                
              }
              stopifnot(!any(is.na(med)))
            }
            
            if (FALSE) {
              for (r in 1:length(med)) {
                plot(x.coord[, r], y.coord[, r], cex=1, cex.axis=1, cex.lab=1,
                     ylim=c(-3,0))
                lines(c(0.5, 0.8), rep(med[r], 2))
                print(med[r])
                rl("")
              }
            }
            
            x.coord <- as.vector(x.coord)
            y.coord <- as.vector(y.coord)
            ylim <- var.value + c(-1, 1) * var.factor
            allowed <- is.finite(y.coord) & y.coord < var.upper.allowed
            P <- y.coord <- y.coord[allowed] 
            x.coord <- x.coord[allowed]
            pch <- rep(16, length(y.coord))
            outside <- y.coord < ylim[1] | y.coord > ylim[2]
            if (any(outside)) {
              pch[outside] <- 4
              y.coord[outside] <- pmin(pmax(ylim[1], y.coord[outside]), ylim[2])
            } 
            fr.name <- paste(var.name, ".relfreq", sep="")
            xlab <- "m(D)/m(0)"
            ylab <- expression(xi(D))
            for (r in 0:2) { # 0:2
              par.name <- paste(name.ps, method, measure, fr.name, sep=".")
              Dev(TRUE, dev, ps=par.name, height=5, width=5, quiet=TRUE)
              par(cex=1, mar=c(2.1, 2.2, 0.1, 0.1))
              plot(x.coord, y.coord, xlab=xlab, ylab=ylab,
                   pch=pch, xlim=xlim, ylim=ylim)
              if (!final) points(rep((med.min + med.max)/2, length(med)), med,
                                 pch=pch, col="red")
              ks <- ksmooth(x.coord, P, ker="box", band=0.1)
              ki <- is.finite(ks$y)
              lines(ks$x[ki], ks$y[ki], col="red", lty=1)
              if(!is.logical(dev)) rl("")
              Dev(FALSE)
              if (is.numeric(dev)) rl(paste(tt$type, "xi=", m$xi, fr.name))
              Dev(TRUE, dev, ps=paste(par.name, "var", sep="."),
                  height=5, width=5, quiet=TRUE)
              ks <- ksmooth(x.coord, (P - var.value)^2, ker="box", band=0.1)
              ki <- is.finite(ks$y)
              plot(ks$x[ki], pmin(upper.bound, ks$y[ki]), type="l",
                   xlab="", ylab="",
                   col="black", xlim=xlim, lty=1,
                   ylim=c(0, upper.bound))
              if (r==0) {
                points((med.min + med.max)/2,
                       sum((med - var.value)^2)/ (length(med) -1 ),
                       pch=pch, col="red")
                smsd <- format(sqrt(mean((med -var.value)^2, na.rm=TRUE)), dig=2)
                txt <- paste(formatC(method, width=7),
                             formatC(measure, width=5),
                             formatC(LETTERS[typ], width=1),
                             formatC(m$xi, width=5),
                             formatC(var.name, width=3), "=",
                             formatC(var.value, dig=3, width=5),
                             "   mean=",formatC(mean(med, na.rm=TRUE),
                                                dig=3, width=8),
                             "   smsd=",formatC(smsd, width=7),
                             "   sd=", formatC(sqrt(var(med, na.rm=TRUE)),
                                               dig=3, width=7)
                             )
                if (is.numeric(dev)) cat(txt,"\n")
                else write(file=result.name, append=TRUE, txt)
              }
              if (!is.logical(dev)) rl("")
              Dev(FALSE)
              if (is.numeric(dev) && r<2)
                rl(paste(tt$type, "xi=", m$xi, fr.name))
              if (r==0) {
                fr.name <- paste(var.name, ".D", sep="")
                x.coord <- distances[allowed]
                xlim <- range(x.coord, na.rm=TRUE)
                xlab <- "D"
              } else {
                ## unused 
                fr.name <-  paste(var.name, "-all.D", sep="")
                y.coord <- P
                ylim <- range(P, na.rm=TRUE)
                pch <- 16  
              } # r==0
            } # for r in 1:2
          } # variable
        } # method
      } # j.p
    } # typ
    if (is.numeric(dev)) close.screen(close.screen())
  } # eval.simu

  data.print <- function(figures, dev=dev, bw=bw, plot.all) {
    height <- 5
    par(bg="white")
    write(file=data.par.name, paste("\\begin{tabular}{c",
            paste(rep("rl", 6), collapse=""), "}"))
    write(file=data.par.name, append=TRUE, nc=100,
          c("&\\multicolumn{4}{c}{threshold parameters}", "&",
            "\\multicolumn{8}{c}{position of the edges of the rectangle}",
            "\\\\"))
    write(file=data.par.name, append=TRUE, nc=100,
          c("profile &",
            paste(paste("\\multicolumn{2}{c}{",
                        c("$c_0$", "$c_2$",  "left",  "right", "bottom", "top"),
                        "}"), collapse="&"),
            "\\\\\\hline"))
    write(file=data.risk.name,
          paste("estimation of xi using ranges of values and",
                "subsequently the median of all xi; see run.R, daten,druck"))
    write(file=data.risk.name, append=TRUE, nc=10,
          c("file      ", "   measure", " raw.idx", "   index", "raw-sens",
            "idx-sens"))
    for (k in 1:length(figures)) {
      if (!is.logical(dev) && k!=1) rl("")
      load(paste(data.dir, figures[k], data.estim.ext, sep=""))
      if (is.null(m[[1]]$picture)) {
        stopifnot(!is.null(m[[1]]$name))
        m[[1]]$picture <- read.picture(m[[1]]$name, Pr=PrintLevel)
       }
      dp <- dim(m[[1]]$picture)
      psname <- paste(pspath, figures[k], sep="")
      if (PrintLevel>1) cat(figures[k], "\n")
      width <- height / dp[2] * dp[1]
      par(mar=data.mar)
      freq <- m[[1]]$r.i$data[, 2]
      dist <- m[[1]]$r.i$data[, 1]
      max.freq <- m[[1]]$r.i$max.freq
      if (plot.all) {
        Dev(TRUE, dev, ps=paste(psname, "lines", sep="-"), hei=height, wid=width)
        plotRGB(m[[1]]$picture, xlab="x [pixels]", ylab="depth [pixels]",
                cex.axis=par()$cex.axis)
        Dev(FALSE)
        if (is.numeric(dev)) rl(paste(figures[k], "lines"))
        else if (PrintLevel>1) cat("lines...")
        
        origfile <- paste(fig.path, "/", figures[k], sep="")
        Dev(TRUE, dev, ps=paste(psname, "orig", sep="-"), hei=height, wid=width)
        plotRGB(read.picture(origfile), xlab="x [pixels]", ylab="depth [pixels]",
                cex.axis=par()$cex.axis)
        Dev(FALSE)
        if (is.numeric(dev)) rl(paste(figures[k], "orig"))
        else  if (PrintLevel>1) cat("orig...")
                
        stained <- m[[1]]$fct(m[[1]]$picture, m[[1]]$param)
        Dev(TRUE, dev, ps=paste(psname, ".stained", sep=""), hei=height,
            wid=width)
        par(mar=c(4.4, 4.8, 0.2, 0.4))
        image(1:nrow(stained), -ncol(stained):-1, stained, xaxs="i", yaxs="i", 
              xlab="x [pixels]", ylab="depth [pixels]",
              col=c("white", if (bw) "black" else "blue"),
              axes=FALSE,  frame.plot=TRUE)
        axis(1)
        pry <- pretty(c(-ncol(stained), -1))
        axis(2, at=pry, labels=-pry)     
        loc <- m[[1]]$loc[[1]]
        loc$x <- round(loc$x) + c(-1, 1) 
        loc$y <- round(loc$y) + c(-1, 1)
        lcol <- if (bw) "slategrey" else "red"
        lines(rep(loc$x[1], 2), loc$y - ncol(stained), lwd=3, lty=5, col=lcol)
        lines(rep(loc$x[2], 2), loc$y - ncol(stained), lwd=3, lty=5, col=lcol)
        lines(loc$x, rep(loc$y[1] - ncol(stained), 2), lwd=3, lty=5, col=lcol)
        lines(loc$x, rep(loc$y[2] - ncol(stained), 2), lwd=3, lty=5, col=lcol)
        Dev(FALSE)
        if (is.numeric(dev)) rl(paste(figures[k], "stained"))
        else if (PrintLevel>1) cat("stained...")
        
        Dev(TRUE, dev, ps=paste(psname, ".pixels", sep=""), hei=height,
            wid=height)
        par(mar=data.mar)
        plot(freq, dist + ncol(stained) - m[[1]]$loc[[1]]$y[2],
             xaxs="i", yaxs="i", pch=16, ylim=c(dp[2], 1), 
             xlab="number of stained pixels", ylab="depth [pixels]",
             col=if (bw) "black" else "blue")
        Dev(FALSE)
        if (is.numeric(dev)) rl(paste(figures[k], "pixels"))
        else if (PrintLevel>1) cat("pixels...")
      } ## plot.all
      idx <- m[[1]]$r.i$selected.dist
      for (ms in 1:length(m)) {
        x.coord <- freq[idx] / max.freq
        xlim <- c(1, 0)
        xlab <- expression(m(D)/m(0))
        y.coord  <- m[[ms]]$r.i$par[1, idx]
        fr.name <- "xi.relfreq"
        col <- rep(if (bw) grey(0.5) else "darkgreen", length(x.coord))
        pch <- rep(if (bw) 20 else 16, length(x.coord))
        Idx <- x.coord<=med.max & x.coord>=med.min
        col[Idx] <- if (bw) "black" else "red"
        pch[Idx] <- 16
        for (r in 0:1) {
          if (PrintLevel>1) cat("r=", r, "")
          par.name <- paste(psname, fr.name, ms, sep=".")
          Dev(TRUE, dev, ps=par.name, height=5, width=5)
          par(mar=data.mar)
          plot(x.coord, pmax(-5, pmin(5, y.coord)), col=col, pch=pch, xlim=xlim,
               xlab=xlab, ylab=expression(xi(D)))
          Dev(FALSE) 
          if (is.numeric(dev)) rl(paste(figures[k], fr.name))
          fr.name <- "xi.D"
          if (!is.null(m[[1]]$r.i$selected.dist)) {          
            x.coord <- m[[1]]$r.i$selected.dist
          } else {
            x.coord <- m[[1]]$contamination.risk$selected.dist
          }
         xlim <- range(x.coord)
          xlab <- "threshold depth D"
          col <- if (bw) "black" else "darkgreen"
          pch <- 16
        } #r
        
        if (is.numeric(dev)) {
          split.screen(c(2,2))
          screen(1)
          plotRGB(m[[1]]$picture, xlab="", ylab="", cex.axis=par()$cex.axis)
          screen(3)
        } else {
          ## raw.hist: histogram contains values close to the border of
          ## the maximisation domaine
          Dev(TRUE, dev, ps=paste(psname, "raw.hist", ms, sep="."))
        }
        # cat("hist...")
        hist(m[[ms]]$raw.risk.index, 30, main="",
             col=if (hist.bw) "grey" else "green",
             xlab=parse(text=paste("'estimated ' * xi")), ylab="frequency")
        median.raw <- median(m[[ms]]$raw.risk.index)
        points(m[[ms]]$r.i$raw.risk, 0, col=if (hist.bw) "black" else "red")
        points(median.raw, 0, pch=20, col=if (hist.bw) "black" else "blue")
        if (is.numeric(dev)) screen(4) else {
          Dev(FALSE)
          Dev(TRUE, dev, ps=paste(psname, "hist", ms, sep="."))
        }
        # cat("hist2...")
        hist(m[[ms]]$risk.index, 30, main="",
             col=if (hist.bw) "lightgrey" else "green",
             xlab=parse(text=paste("'estimated ' * xi")), ylab="frequency")
        median.index <- median(m[[ms]]$risk.index)
        points(m[[ms]]$r.i$risk, 0, pch=2, cex=2,
               col=if (hist.bw) "black" else "red")
        points(median.index, 0, pch=16, cex=2,
               col=if (hist.bw) "black" else "blue")
        if (is.numeric(dev)) {
          if (k!=length(figures) || ms!=length(m)) rl("")
        } else Dev(FALSE)
        
        write(file=data.risk.name, append=TRUE, nc=10,
              c(formatC(figures[k], width=10),
                formatC(paste(as.character(as.list(m[[ms]]$r.i$measure))[-1],
                              collapse=";"), width=10),
                formatC(m[[ms]]$r.i$raw.risk, dig=4, width=8),
                formatC(m[[ms]]$r.i$risk, dig=4, width=8),
                formatC(median.raw, dig=4, width=8),
                formatC(median.index, dig=4, width=8))
              )
      } # ms
      
      plot.interval <- function(rg)
        if (diff(rg)==0) "[ -- ]" else paste("[",rg[1], ", ", rg[2],"]", sep="")
      
      ## reversed for y -- document it!!!!
      interval <- function(xy, i) {
        loc <- m[[1]]$loc
        rg <- round(range(loc[[2]][xy][[1]][i], loc[[3]][xy][[1]][i]))
        lc <- loc[[1]][xy][[1]][i]
        if (xy=="y") {
          rg <- rev(dim(m[[1]]$picture)[2] - rg + 1)
          lc <- dim(m[[1]]$picture)[2] - lc
        }
        c(round(lc), "& \\kern-3mm", plot.interval(rg))       
      }
      
      interval.par <- function(par) {
        rg <-  round(range(m[[1]]$lower[par][[1]] ,m[[1]]$upper[par][[1]]))
        c(round(m[[1]]$param[par][[1]]), "&\\kern-3mm", plot.interval(rg))
      }
      
      write(file=data.par.name, append=TRUE, nc=100,
            c(letters[k], " %", figures[k], "\n&", 
              interval.par("minR"), "&", interval.par("maxGB"), "&",
              interval("x", 1), "&", interval("x", 2), "&",
              interval("y", 1), "&", interval("y", 2), "\\\\"))
    } ## figures
    write(file=data.par.name, append=TRUE, nc=100, "\\end{tabular}")
    if (is.numeric(dev)) close.screen(close.screen())
  } # data.print


 pictures <- function(dev, repet=1, types, models, name, estimate=FALSE){
   assign(".Random.seed", seed, envir=.GlobalEnv)
   if (is.numeric(dev)) {
     split.screen(c(2,2))
     screen(1, new=FALSE)
   }
   xi <- numeric(length(model.param))
   tt <- list()
   for (typ in 1:length(type)) {
     est <- list()
     zaehler <- 0
     for (nt in names(type[[typ]])) 
       eval(parse(text=paste("tt$", nt, "<-", "type[[typ]]$", nt)))
     if (!any(typ %in% types)) next
     m <-  list()
     for (w in 1:repet) {
       old.simu <- NULL
       for(m.p in 1:length(model.param)) {
         for (n in names(model.param[[m.p]])) 
           eval(parse(text=paste("m$", n, "<-", "model.param[[m.p]]$", n)))
         if (!any(m.p %in% models)) next
         environment(m$dd) <- environment()        
         psname <- paste(paste(pspath, name, sep=""), tt$type, tt$lambda * 100, 
                         round(tt$x.var * 100), m$xi, sep=".")
         if (repet>1) psname <- paste(psname, w, sep=".")
         for (ms in 1:length(Measures)) {
           est <- flowpattern(type=tt$type,
                              x.h.scale=tt$x.h.scale, x.v.scale=tt$x.v.scale,
                              x.var=tt$x.var,
                              length=tt$length, lambda=tt$lambda,
                              delta.x=tt$delta.x, delta.y=tt$delta.y,
                              endpoint.tolerance=tt$endpoint.tolerance,  
                              unif.b=tt$unif.b,
                              drop.distr=m$dd, drops=m$drops,
                              drop.sc=m$drop.scale,
                              method=if (estimate) Methods else NULL, #depend
                              select=select.dist,
                              measure=Measures[[ms]],
                              old=old.simu,
                              PrintLevel=PrintLevel)
           old.simu <- est$intermediate
           est$intermediate <- est$input$old.paths <- NULL

           ## print
           if (ms==1) {
             max.freq <- max(est$freq)
             Dev(TRUE, dev, ps=psname, height=2.5, width=5)
             par(cex=1, mar=c(4.1, 4.5, 0.3, 1.4))
             plot(est$i.x, # * if (final) 1 else 100,
                  -est$i.d, ylim=c(-est$input$depth, -1),
                  pch=".", xlab="x [pixels]", ylab="depth [pixels]",
                  main="", cex=1/2, col=if (bw) "black" else "blue")
             Dev(FALSE)
             if (is.numeric(dev)) rl(psname)
             
             Dev(TRUE, dev, ps=paste(psname, "freq", sep="."),
                 height=5, width=5)
             par(cex=1, mar=c(4.1, 4.5, 0.3, 0.5))
             plot(est$freq, -1:-est$input$depth, ylim=c(-est$input$depth, -1),
                  pch=".", xlab="# pixels", ylab="depth [pixels]",
                  main="", cex=1/2, col=if (bw) "black" else "blue")
             Dev(FALSE)
             if (is.numeric(dev)) rl(psname)
            }
           
           if (length(est$input$method)>0) {
             mea <- names(Measures)[ms]
             for (me in 1:length(est$input$method)) {
               idx <- est[[me]]$selected.dist
               param <- est[[me]]$par[1, idx]         
               x.coord <- est$freq[idx]  / max.freq
               xlim <- c(1, 0)
               y.coord <- param
               ylim <- m$xi + c(-1, 1) * 2
               pch <- rep(16, length(y.coord))
               outside <- (!is.finite(y.coord) | y.coord < ylim[1] |
                           (y.coord > ylim[2] & y.coord < upper.pos.xi))
               if (any(outside, na.rm=TRUE)) {
                 pch[outside] <- 4
                 y.coord[outside] <-
                   pmin(pmax(ylim[1], y.coord[outside]), ylim[2])
               } 
               fr.name <- "xi.relfreq"
               xlab <- "m(D)/m(0)"
               ylab <- expression(xi(D))
               col <- rep(if (bw) "black" else "darkgreen", length(x.coord))
               col[x.coord<=med.max & x.coord>=med.min] <-
                 if (bw) "black" else "red"
               for (r in 0:2) { 
                 par.name <-
                   paste(psname, fr.name, est$input$method[me], mea, sep=".")
                 Dev(TRUE, dev, ps=par.name, height=3, width=3)
                 par(cex=1, mar=c(4.1, 4.5, 0.3, 0.5))
                 plot(x.coord, y.coord, xlab=xlab, ylab=ylab, pch=pch, 
                      xlim=xlim, ylim=ylim, col=col)
                 Dev(FALSE)
                 if (is.numeric(dev))
                   rl(paste(mea, est$input$method[me], fr.name))
                 col <- if (bw) "black" else "darkgreen"
                 if (r==0) {
                   fr.name <- "xi.D"
                   x.coord <- est$dist[idx]
                   xlim <- range(x.coord)
                   xlab <- "D"
                 } else {
                   fr.name <- "xi-all.D"
                   y.coord <- param
                   ylim <- range(y.coord[y.coord < upper.pos.xi], na.rm=TRUE)
                   pch <- 16  
                 }
               } # r
             } #me + length>0
           }
           x.coord <- y.coord <- e <- NULL
           gc()
         } # ms
       } ## m.p
     } # repet
   } #  typ
   close.screen(close.screen())
 }


  sorry <- function() {
    cat("this functionality is currently only available on Linux/Unix systems\n")
  }

#######################################################################
  items <-
    c("Fig. 1,2,4,5 (sketches)",
      "Fig. 13 -- it also gives you a time estimate for the next menu point",
      "Simulation of profiles for last paragraph of Section 4", #3
      "Risk indices of profiles for last paragraph of Section 4",
      "(Re-)Evaluation of indices of profiles for last paragr., Sec.4",
      "simulation study",    #6
      "use all estimator variants in the simulation study (time demanding!)",
      "evaluation of simul. study ('simulation study' must be called first) (Tab.1)",
      "data image analysis", #9
      "data risk estimation ('data image analysis' must have been called first)",
      "data print ('data risk estimation' must have been called first) (Fig. 10, Table 2, Fig. 11,12)",
      "Fig. 3",              #12
      "Fig. 6,8,9",
      "Fig. 7"
               )
  while (TRUE) {
    RFparameters(Print=PrintLevel, pch="", TBM3.linesimus=linesimustep,
                  TBM2.linesimus=linesimustep, TBM3.linesimufactor = 0.0,
                 TBM2.linesimufactor = 0.0)
    if (length(input)==0) input <- menu(items)
    if (input[1] > length(items) || input[1]==0) break
    if (PrintLevel>1) cat(input[1], ":", items[input[1]], "\n")

    ut <- unix.time(
    switch(input[1],
             {
               sketches(dev=dev, bw=bw)
             }, {
               sensitivity.simu(final=final, repet=1, dev=dev, bw=bw)
             }, { #3
               if (length(input)>1 ||
                   (wiederhol <-
                    readline(paste("number of repetitions (",
                                   sens.repet, ") -- Press return"))) == "")
                 wiederhol <- sens.repet
               wiederhol <- as.integer(wiederhol)
               sensitivity.simu(final=final, repet=wiederhol, dev=dev, bw=bw)
             }, {
               if (length(input)>1 ||
                   (wiederhol <-
                    readline(paste("number of repetitions (",
                                   sens.boot.repet, ") -- Press return"))) == "")
                 wiederhol <- sens.boot.repet
               wiederhol <- as.integer(wiederhol)
               sensitivity.simu.eval(bootstrap.repet=wiederhol)
               input <- c(NA, 5, input[-1])
             }, {
               sensitivity.simu.print(dev=dev, bw=bw)
             }, { #6
               if (length(input)>1 ||
                   (wiederhol <-
                    readline(paste("number of repetitions (",
                                   simu.repet, ") -- Press return"))) == "")
                 wiederhol <- simu.repet
               simulation(repet=wiederhol)
#               simulation(repet=wiederhol, 1:8)
             }, {
               all.estimate()
             }, {
               eval.simulation(dev=dev)
             }, { #9
               # if (.Platform$OS.type!="unix") {sorry(); next} 
               for (i in figures) {
                 data(list=i, envir=environment(), package="SoPhy")
                 l <- get(i)
                 l$name <- paste(fig.path, "/", i, ".G", sep="")
                 l <- analyse.profile(l, estimate.all=NULL, tit=i,
                                      PrintLevel=PrintLevel,
                                      endpoint.tolerance=data.endpoint.tolerance)
                 if (!keep.figures) l$picture <- NULL
                 save(file=paste(data.dir, i, data.ana.ext, sep=""), l)
                 dev.off()
               }
             }, {
               # if (.Platform$OS.type!="unix") {sorry(); next} 
               for (k in figures) {
                 load(paste(data.dir, k, data.ana.ext, sep=""));
                 m <- list()
                 for (mi in 1:length(data.measures)) {
                   if (PrintLevel>1) cat(k, data.measures[mi], "\n")
                   m[[mi]] <-
                     analyse.profile(l, estimate.all=TRUE, interactive=FALSE,
                                     selected.dist=data.selected.dist,
                                     selected.rate=selected.rate,
                                     Print=PrintLevel,
                                     ppi.par=data.ppi.par,
                                     measure=Measures[data.measures[mi]][[1]],
                                     endpoint.tolerance=data.endpoint.tolerance,
                                     method=data.method
                                     )
                   if (!keep.figures) m[[mi]]$picture <- NULL
                   save(file=paste(data.dir, k, data.estim.ext,sep=""),
                        m)
                 }
               }
             }, {
               # if (.Platform$OS.type!="unix") {sorry(); next} 
               data.print(figures, dev=dev, bw=bw, plot.all=data.plot.all)
             }, { #12
               pictures(dev=dev, types=1:6, models=1, name=pict.name)
             }, {
               pictures(dev=dev, types=5, models=1:5, name=pict.name, estim=TRUE)
             }, {
               pictures(dev=dev, repet=pict.repet, types=6, models=5,
                        name=pict.several.name)
             }
             ))
    if (!final) cat("\nsystem time for #", input[1], ":", ut,"\n\n")
    input <- input[-1]
  } # while (true)
}


######################################################################
######################################################################
######################################################################
  
sh.jh <- function(input=NULL, dev=2, pspath="./", final=TRUE, PrintLevel=0,
                  readlines=TRUE, low.resolution=TRUE, bw=TRUE
                  ) {
  blue <- if (bw) 1 else "#0000FF"
  hex <- function(x) {
    h <- c(0:9, LETTERS[1:6])
    paste(h[1 + x / 16], h[1 + x %% 16], sep="")
  }
  default.first <- FALSE
  #  default.first <- TRUE
  ## basic definitions
  if (final) {
    lambda <- 0.7         ## in case of simulation on Poisson points
    len.x <- len.y <- 101 ## in case of simulation on grid
    nlines <- 250         ## TBM3D method in GaussRF
  } else {
    lambda <- 0.002       ## in case of simulation on Poisson points
    len.x <- len.y <- 10  ## in case of simulation on grid
    nlines <- 10          ## TBM3D method in GaussRF
  }
  length.profile <- 100
  drop.scale <- 15
  default.distrib <- 4
  radius <- 2
  fullsize <- FALSE
  delta <- 10
  linesimustep <- 0.1

  nu1 <- 0.3 # to investigate the influence of the parameter nu 
  nu2 <- 1
  nu3 <- 8
  f1 <- GetPracticalRange("whittle", nu1) / GetPracticalRange("whittle", 2)
  f2 <- GetPracticalRange("whittle", nu2) / GetPracticalRange("whittle", 2)
  f3 <- GetPracticalRange("whittle", nu3) / GetPracticalRange("whittle", 2)


  drop.distr <- list(x2 = function(x) x^0.25 * 100,
                     x05 = function(x) x^2 * 120,
                     exp2 = function(x) (-log(1-x))^(1/2) * 50,  ##
                     exp1 = function(x) -log(1-x) * 25,  ##
                     xM2 = function(x) (1-x)^(-0.5) * 15,
                     xM1 = function(x) (1-x)^(-1) * 3,
                     fix = function(x) 0 * x + 400
                     )

  default <- list(type="all", 
                  length.profile=length.profile,
                  width.slice=100,
                  delta.y=delta, delta.x=delta, 
                  lambda=lambda,
                  depth=100,
                  x.var=0.25, y.var=0.25,
                  x.h.scale=20, x.v.scale=5, x.kappa=2,  
                  drop.scale=drop.scale,
                  drop.kappa=2,
                  drop.distr=drop.distr[default.distrib],
                  n.balls=1,      ## "pixel" resolution
                  Horizont="no", ## "", ""
                  Profiles=1, radius=radius,
                  grid=FALSE, raw=FALSE,
                  unit="cm", unit.scale=1,
                  nxt=-1
                  )

  crestana.resolution <- 5
  crestana.length <- 30
  crestana.width <- 1
  crestana.delta <- 0
  schwartz.resolution <- 5
  schwartz.length <- 100
  schwartz.width <- 6
  schwartz.delta <- 2 * schwartz.resolution
  type <- list(default, #1
               standard = list(Profiles=9,
                 drop.distr=drop.distr[if (default.first)
                   c(default.distrib, 1, 2, 6) else
                   c(1, 2, default.distrib, 6)]
                 ),          
               default, #3
               standardgrid = list(raw=TRUE, grid=TRUE, Profiles=1), 
               default, #5
               dropkappa1 = list(drop.kappa=nu1, drop.scale=drop.scale / f1),
               dropkappa2 = list(drop.kappa=nu2, drop.scale=drop.scale / f2),
               dropkappa3 = list(drop.kappa=nu3, drop.scale=drop.scale / f3),
               default, #9
               xkappa1 = list(x.kappa=nu1, x.h.scale=20 / f1, x.v.scale=5 / f1),
               xkappa2 = list(x.kappa=nu2, x.h.scale=20 / f2, x.v.scale=5 / f2),
               xkappa3 = list(x.kappa=nu3, x.h.scale=20 / f3, x.v.scale=5 / f3),
               default, #13
               absorbierend= list(Horizont="absorbing", n.balls=5),
               durchbruch = list(Horizont="breakthrough", n.balls=5),
               default, #16
               dropscale1 = list(drop.scale=5), 
               dropscale2 = list(drop.scale=30),          
               default, #19
               vscale1 = list(x.v.scale=2),
               vscale2 = list(x.v.scale=15),
               default, #22
               hscale1 = list(x.h.scale=5),
               hscale2 = list(x.h.scale=40),
               default, #25
               xvar1 = list(x.var=0.1, y.var=0.1),
               xvar2 = list(x.var=1, y.var=1),
               xvar0 = list(x.var=10^{-12}, y.var=10^{-12}), #
               xvar00grid = list(raw=TRUE, grid=TRUE),
               default, #30
               xvar0grid = list(y.var=10^{-8}, raw=TRUE, grid=TRUE),
               default, #32
               crestana = list(
                 length.profile=crestana.length * crestana.resolution,
                 width.slice=crestana.width * crestana.resolution,
                 delta.y=crestana.delta * crestana.resolution,
                 delta.x=crestana.delta * crestana.resolution, 
                 lambda= (7 / (crestana.length + 2 * crestana.delta) /
                          (crestana.width + 2 * crestana.delta) /
                          crestana.resolution^2),    
                 depth=100 * crestana.resolution,
                 x.var=0.0028 * crestana.resolution^2, y.var=0,
                 x.h.scale=0.4 * crestana.resolution,
                 drop.scale=1 * crestana.resolution,
                 drop.distr= list(
                   crestana=function(x) x * 100 * crestana.resolution),
                 radius=0.5,
                 unit="cm", unit.scale=1 / crestana.resolution
                 ),
               default, #34
               schwartz = list(
                 length.profile=schwartz.length * schwartz.resolution,
                 width.slice = schwartz.width * schwartz.resolution,
                 delta.y=schwartz.delta * schwartz.resolution,
                 delta.x=schwartz.delta * schwartz.resolution,  # in pixels
                 lambda = (10000 / (schwartz.length + 2 * schwartz.delta)
                           / (schwartz.width + 2 * schwartz.delta) /
                           schwartz.resolution^2),
                 depth=50 * schwartz.resolution,
                 x.var=0.06 * schwartz.resolution^2,
                 x.h.scale=4 * schwartz.resolution,
                 x.v.scale=1 * schwartz.resolution,
                 drop.scale=10 * schwartz.resolution,
                 drop.kappa=2,
                 ##drop.distr= function(x) (-log(1-x))^0.7 * 170,
                 drop.distr= list(schwartz=function(x) (-log(1-x))^0.6 * 30 *
                   schwartz.resolution),
                 radius=0.2, ## in true units
                 unit="cm", unit.scale=1/schwartz.resolution),
               schwartz2 = list(
                 drop.distr= list(schwartz2=function(x) (1-x)^(-0.3) * 50),
                 drop.scale=2 * schwartz.resolution,
                 x.var=0.004 * schwartz.resolution^2,
                 x.h.scale=2 * schwartz.resolution,
                 x.v.scale=20 * schwartz.resolution
                 ),
               default #36
               )

  slct <- c(1, 3, 5, 9, 13, 16, 19, 22, 25, 30, 32, 34)
  inf <- c(length.profile * 1.6, -400, 1/500) # x,y, dehnungsfaktor
  sun <- c(200, -30000, 50)
  cex.axis <- 2
  mai <- c(1.2, 1.2, 0.2, 0.3)
  profileheight <- 4
  totalheight <- profileheight + mai[1] + mai[3]
  colour <- if (bw) grey(seq(1, 0, len=100)) else rainbow(100)
  reverse <- TRUE ## labelling of the y axis
  quiet <- PrintLevel>1
  
  ##selected <-  FALSE # selected for historical reasons
  selected <- TRUE
 
  close.screen(close.screen())
  par(bg="white")
  if (!exists(".Random.seed",  envir=.GlobalEnv)) runif(1)
  seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
  RFparameters(Print=1, pch="", PracticalRange=11,
               TBM2.num=TRUE, TBMCE.forc=TRUE,
               TBM3.lines=nlines, CE.useprimes=TRUE,
               TBM3.linesimufactor = 0.0,
               TBM2.linesimufactor = 0.0,
               TBM3.linesimus=linesimustep, 
               TBM2.linesimus=linesimustep)
  for (name in c(pspath)) {
    if (PrintLevel>2) cat("path", name, "\n")
    if (file.exists(name)) stopifnot(file.info(name)$isdir)
    else {
      if (.Platform$OS.type=="unix") stopifnot(dir.create(name))
      else try(dir.create(name))
    }
  }
  if (substr(pspath, nchar(pspath), nchar(pspath)) != "/")
    pspath <- paste(pspath,  "/", sep="")
  rl <- if (readlines) # || .Platform$OS.type!="unix")
    function(x) if (x=="") readline("press return")
    else readline(paste(x, ": press return"))
  else function(x) { cat(x); sleep.milli(1000); cat("\n")}

  sketches <- function(dev) {
    x <- seq(0, 100, 10)
    incr <- GaussRF(x, model="exp", param=c(0,4000,0,5), grid=TRUE)
    p <- cumsum(incr)
    Dev(TRUE, dev=dev, ps=paste(pspath, "path.ps", sep=""))
    par(mar=c(4, 4.5, 0.2, 0.2))
    plot(p, x, type="l", axes=FALSE, frame=TRUE, cex=2,
         ylab="depth", xlab="x [units]",
         cex.lab=2, lwd=3, col=if (bw) 1 else "brown")
    axis(1, cex.axis=2)
    axis(2, at=seq(0, 100, 20), labels=seq(100, 0, -20), cex.axis=2)
    Dev(FALSE)
    if (is.numeric(dev)) rl("first sketch")

    Dev(TRUE, dev=dev, ps=paste(pspath, "intregr.path.ps", sep=""))
    par(mar=c(4, 4.5, 0.2, 1.3))
    ip <- cumsum(sqrt(diff(x)^2 + diff(p)^2))
    plot(c(0, x[-1]), c(0,ip), type="l", axes=TRUE, frame=TRUE, cex=2,
         xlab="depth", ylab="integrated path length",
         cex.lab=2, cex.axis=2, lwd=3, xaxs="i", yaxs="i",
         col=if (bw) 1 else "brown")
    D <- 450
    i <- sum(ip <= D) + 1
    xD <- x[i] + (D-ip[i-1])/(ip[i] - ip[i-1]) * (x[i+1] - x[i])
    points(xD, D, col="red")
    text(-1, D, "D", adj=c(1,0.35), xpd=TRUE, cex=2)
    lines(c(0, xD), c(D, D), lty=2)
    lines(c(xD, xD), c(0, D), lty=2)
    text(xD, -1, "z", adj=c(0.85, 1), xpd=TRUE, cex=2)
    Dev(FALSE)
  }
  
#######################################################################

  while (TRUE) {
    if (length(input)==0) {
      input <-
        menu(c("Fig. 4/5",
               "Fig. 2b/c/f/g/h, 1b", #
               "Fig. 10",
               "nu (xy)",
               "Fig. 11",
               "Fig. 6",
               "Fig. 8",
               "Fig. 9",
               "Fig. 7 u 3a", #
               "Fig. 2d/e", #
               "Crestana & Posadas (1998, Fig. 8)", #
               "Schwartz et al. (1999, Fig. 4)", # 
               "Fig. 3 (sketches)",#
               )
             )
    }
    if (PrintLevel>1) cat("input", input[1], "\n")
    if (input[1]>length(slct)) {
      switch(input[1]-length(slct),
             {
               sketches(dev=dev)
               input <- input[-1]
               next
             }
             )
    }
    if (input[1] > length(slct) || input[1] < 1) break
    tt <- list()
    nxt <- 0
    firstplot <- TRUE
    for (typ in slct[input[1]]:length(type)) {
      ## length(type) is above is dummy nxt=1 breaks, see below and
      ## the definition of type above
      if (PrintLevel>1) cat("typ",typ, "\n")
      for (nt in names(type[[typ]])) {
        txt <- paste("tt$", nt, "<-", "type[[typ]]$", nt)
        eval(parse(text=txt))
      }
      if (!is.null(type[[typ]]$nxt)) {
        nxt <- nxt + 1
        if (nxt==1) next else break
      }
      if (firstplot) firstplot <- FALSE
      else if (is.numeric(dev)) rl("") 
      
      if (is.numeric(dev)) while(!is.null(dev.list())) dev.off()
      assign(".Random.seed", seed, envir=.GlobalEnv)
      res <- flowpattern(type=tt$type,
                         x.h.scale=tt$x.h.scale,
                         x.v.scale=tt$x.v.scale,
                         depth=tt$depth,
                         x.var=tt$x.var,
                         y.var=tt$y.var,
                         x.kappa=tt$x.kappa,
                         length.profile=tt$length.profile,
                         lambda=tt$lambda,
                         delta.x=tt$delta.x,
                         delta.y=tt$delta.y,
                         width.slice=tt$width.slice,
                         endpoint.tolerance=NA,  
                         unif.b=NA,
                         drop.distr=function(x) x, drops=1,
                         drop.sc=tt$drop.scale,
                         drop.kappa=tt$drop.kappa,
                         grid=tt$grid,
                         len.x = len.x, len.y=len.y,
                         method=NULL,
                         simu.method = "TBM3",
                         drop.simu.method = "spectral",
                         PrintLevel=PrintLevel,
                         raw=tt$raw
                         )
      if (PrintLevel>1) cat("paths simulated\n")
      DeleteAllRegisters()
      inp <- res$input
      idx <- res$interm$x >= 0 & res$interm$x <= inp$length.profile
      idy <- res$interm$y >= 0 & res$interm$y <= inp$width.slice
      iy <- sum(!idx) / 2 +
        if (tt$Profiles==1) (len.y + 1) / 2 else seq(1, len.y, len=tt$Profiles)
      totalwidth <- profileheight / inp$depth * inp$length.profile +
        mai[2] + mai[4]
      if (tt$raw) {
        stopifnot(tt$Profiles >= 1, tt$grid)
        zlim <- range(res$interm$raw.xx)
        for (cs in c(FALSE, TRUE)) {
          for (ii in 1:tt$Profiles) {
            if (PrintLevel>1) cat("cumsum=", cs, "profile=", ii, "\n")
            i <- iy[ii]
            base.name <- paste("2d.pattern.xz", names(type[typ]), ii, sep=".")
            Dev(TRUE, dev, ps=paste(pspath, base.name, ".", cs, sep=""),
                height=totalheight, width=totalwidth, quiet=quiet)
            par(mai=mai)
            cum <- t(apply(res$interm$raw.xx[idx, i, ], 1, cumsum))
            ## t() da apply(, 1, .) zeilen und spalten vertauscht
            cszlim <- range(cum)
            image(x=seq(0, inp$length.profile, len=len.x),
                  y=seq(1, inp$depth, 1), xaxs="i", yaxs="i",
                  if (cs) cum[, inp$depth:1]
                  else res$interm$raw.xx[idx, i, inp$depth:1],
                  col=colour, zlim=if (cs) cszlim else zlim,
                  axes=FALSE, frame=TRUE, xlab="x [cm]", ylab="z [cm]",
                  cex.lab=cex.axis)
            lab <- pretty(c(0, inp$depth))
            at <- if (reverse) {inp$depth - lab} else {lab}
            axis(1, cex.axis=cex.axis)
            axis(2, at=at , lab=lab, cex.axis=cex.axis)
            my.legend(lb.x=0, lb.y=2, col=colour, zlim=if (cs) cszlim else zlim,
                      cex=1.5)
            Dev(FALSE)
            if (is.numeric(dev)) rl(base.name)
          }
        }
      }
      
      for (dd in 1:length(tt$drop.distr)) {
        base.name <- paste("3d", names(type[typ]),
                           names(tt$drop.distr[dd]), sep=".")
        if (PrintLevel>1)
          cat("drop distr:", as.character(as.list(tt$drop.distr[[dd]])),"\n")
        x <- plotFlow3d(res, tt$Horizont, n.balls=tt$n.balls,
                        drop.distr=tt$drop.distr[[dd]],
                        ps=paste(pspath, base.name, sep=""), dev=dev,
                        pointradius=tt$radius,
                        ps.background=FALSE,
                        inf=inf, sun=sun, rl=rl,
                        unit=tt$unit, unit.scale=tt$unit.scale,
                        low.resolution = low.resolution,
                        col = if (bw) grey(pmin(1, pmax(0, seq(0.95,0,-0.001))))
                        else paste("#0000", hex(seq(0, 255,
                          length=100)),sep="")
                        )
        if (is.null(tt$Profile) || tt$Profiles<1) {
          if (PrintLevel>1) cat("profiles jumped \n")
          next;
        }
        plotFlow2d(x, ps=paste(pspath, base.name, ".profile.", sep=""),
                   dev=dev, Profiles=tt$Profiles,
                   pointradius=tt$radius, unit=tt$unit, 
                   full.size=fullsize, rl=rl, col=blue)
        if (is.numeric(dev)) rl(paste(base.name, "profile"))
        
        if (tt$raw) {
          plotFlow2d(x, ps=paste(pspath, base.name, ".profile.raw.", sep=""),
                     dev=dev, Profiles=tt$Profiles, pointradius=tt$radius,
                     full.size=FALSE, slice=0.1, rl=rl,
                     col=blue)
          if (is.numeric(dev)) rl(paste(base.name, "profile.raw"))
        }
        lab <- pretty(c(0, inp$depth))
        at <- if (reverse) {inp$depth - lab} else {lab}
        stopifnot(tt$Profiles > 0)
        Dev(TRUE, dev, ps=paste(pspath, "2d.drop.", names(type[typ]), ".",
                         names(tt$drop.distr[dd]), sep=""),
            height=totalheight, width=totalwidth, quiet=quiet)
        par(mai=mai)
        z <- apply(tt$drop.distr[[dd]](res$interm$LEN), 1, max)
        zlim.drop <- range(z)
        if (tt$grid) {        
          z <- matrix(z, nrow=inp$len.x)[idx, idy]
           image(x=seq(0, inp$length.profile, len=len.x),
                y=seq(0, inp$width.slice, len=len.y),
                z, col=colour, zlim=zlim.drop, axes=TRUE, frame=TRUE,
                cex.axis=cex.axis, cex.lab=cex.axis,
                xaxs="i", yaxs="i", xlab="x [cm]", ylab="y [cm]")
        } else {
          idxy <- (res$interm$xx[,1] >= 0 &
                   res$interm$xx[,1] <= inp$length.profile &
                   res$interm$yy[,1] >= 0 &
                   res$interm$yy[,1] <= inp$width.slice)
          plot(x=res$interm$xx[idxy, 1], xlim=c(0, inp$length.profile),
               y=res$interm$yy[idxy, 1], ylim=c(0, inp$width.slice), 
               col=colour[(z[idxy] - zlim.drop[1]) / diff(zlim.drop) *
                 (length(colour) - 1) + 1],
               axes=TRUE, frame=TRUE, pch=20,
               cex.axis=cex.axis, cex.lab=cex.axis,
               xaxs="i", yaxs="i", xlab="x [cm]", ylab="y [cm]")
        }
        my.legend(lb.x=0, lb.y=2, col=colour, zlim=zlim.drop, cex=1.5)
        Dev(FALSE)
        if (is.numeric(dev)) rl(paste(base.name, "2d.drop"))
         
        if (tt$raw) {
          iy <- if (tt$Profiles==1) (len.y+1)/2 else seq(1,len.y,len=tt$Profiles)
          for (ii in 1:tt$Profiles) {
            i <- iy[ii]
            Dev(TRUE, dev, ps=paste(pspath, "2d.drop.xz.", names(type[typ]), ".",
                             names(tt$drop.distr[dd]), ".", ii, sep=""),
                height=totalheight, width=totalwidth, quiet=quiet)
            par(mai=mai)
            plot(Inf, Inf, xlim=c(0, inp$length.profile),
                 ylim=c(-inp$depth, 0), axes=FALSE, frame=TRUE,
                 xaxs="i", yaxs="i", cex.axis=cex.axis,
                 xlab="x [cm]", ylab="z [cm]", cex.lab=cex.axis)
            lab <- pretty(c(0, inp$depth))
            at <- -lab
            axis(1, cex.axis=cex.axis)
            axis(2, at=at , lab=lab, cex.axis=cex.axis)
            polygon(x=c(seq(0, inp$length.profile, len=len.x),
                      rep(inp$length.profile * 1.1, 2),
                      rep(-5, 2)),
                    y=c(-z[,i], -z[nrow(z), i], 5, 5, -z[1,i]),
                    col=blue)
            Dev(FALSE)
            if (is.numeric(dev)) rl(paste(base.name, "raw", ii))
          } ## ii
        } ## tt$raw
      } # length(drop.distr)
    } # typ
    input <- input[-1]
  } # while (true)
}

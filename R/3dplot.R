pos3D <- function(x, inf) {
  if (missing(inf)) inf<- get(".p3d.inf", env=.ENV)
  x <- as.matrix(x)
  stopifnot(nrow(x)==3)
  p <- rbind(x[1, ] + (inf[1] - x[1, ]) * (1 - exp(-inf[3] * x[2, ])),
             x[3, ] + (inf[2] - x[3, ]) * (1 - exp(-inf[3] * x[2, ])))
  return(p)
}

# initperspective <- function() {}

quader <- function(size, inf, sun,
                   bottomleftfront= c(0,0,0),  #bottom front left corner of
                                        # the brick
                   col=grey(seq(1,0,-0.001)),
                   col.frame=c('grey', 'black'), # invisible, visible
                   lty=c(2,1), # invisible, visible
                   cex.axis=1.5, reverse=TRUE,
                   unit = "cm",
                   add=FALSE, plot=TRUE,
                   font = c("sans serif", "bold"), srt=NULL
                   ){
  lab.line <- 0.9
  cex.axis <- cex.axis * 0.6
  lty <- rep(lty, len=2)
  if (add) { ## add to existing quader plot 
    inf <- get(".p3d.inf", env=.ENV)
    if (missing(size)) size <- get(".p3d.profile", env=.ENV)
    if (missing(bottomleftfront))
      bottomleftfront <- get(".p3d.bottomleftfront", env=.ENV)
  }

  blf <- bottomleftfront
  ## qq1: vordere 4 ecken 
  qq1 <- rbind(blf[1] + c(0,size[1], size[1], 0),
              blf[2],
              blf[3] + c(0,0,size[3], size[3]))
  q1 <- pos3D(qq1, inf)
  ## qq2: hintere 4 ecken 
  qq2 <- rbind(qq1[1,], blf[2] + size[2], qq1[3,])
  q2 <- pos3D(qq2, inf)
  xlim <- range(q1[1, ], q2[1, ])
  ylim <- range(q1[2, ], q2[2, ])
  
  if (plot) {
    signx <- 2 * (min(q1[1, ]) < min(q2[1, ])) - 1
    signy <- signz <- 1 - 2 * (min(q1[2, ]) < min(q2[2, ]))
    
    if (!add) {
      plot(Inf, Inf,  xlim = xlim, ylim = ylim, xlab="", ylab="",
           cex.lab=cex.axis, xaxs="i", yaxs="i",
           frame.plot=FALSE, axes=FALSE)
      
      ## x Achse, Beschriftung
      if (min(q1[2, ]) < min(q2[2, ])) {     
        at <-  range(q1[1, ])
        xlab.pos <- pos3D(c(blf[1] + size[1] / 2, blf[2], blf[3]), inf)
      } else {
        at <- range(q2[1, ])
        xlab.pos <- pos3D(c(blf[1] + size[1] / 2, blf[2] + size[2], blf[3]), inf)
      }
      lab <- pretty(blf[1] + size[1] * 0:1)
      at <- at[1] + diff(at) * (lab - blf[1]) / size[1]
      axis(1, at=at, lab=rep("", length(at)))
      text(at, xlab.pos[2] - 1 * par()$cxy[2], lab,
           adj=c(0.5, 1), xpd=TRUE, cex=cex.axis, vfont=font)
      if (unit!="") {
         text(xlab.pos[1], xlab.pos[2] - 3 * lab.line * par()$cxy[2],
              paste("x [", unit, "]", sep=""),
              adj=c(0.5, 1), xpd=TRUE, cex=cex.axis, vfont=font)
      }
      
      ## y Achse, Beschriftung
      lab <- pretty(blf[2] + size[2] * c(0,1) , n=4)
      size.units <- (lab - blf[2]) / size[2]
      idx <- size.units>0.03& size.units<0.97
      lab <- lab[idx]
      size.units <- size.units[idx]
      #print(blf[2] + size[2] * c(0.03,0.97)); print(cbind(lab, size.units))
      
      at <- rbind(signx * min(signx * qq1[1, ]) , lab,
                  signz * min(signz * qq1[3, ]))
      at <- pos3D(at, inf)
      ypart <- 40     
      atx <- at[1, ] - signx * size[1]/ypart
      aty <- at[2, ] - signy * size[3]/ypart
      adj <- c(0.5, 0.5 + signy * 0.6) 
       
      matplot(rbind(at[1, ], atx), rbind(at[2, ], aty),
              type="l", add=TRUE, col=1, lty=1, lwd=1)
      if (is.null(srt)) {
        srt <- 180 / pi * atan(diff(at[2, 1:2])/ diff(at[1, 1:2]))
      }
      text(atx, aty, adj=adj, lab=lab, cex=cex.axis, xpd=TRUE, vfont=font,
           srt=srt) 
      idx <- as.integer(length(atx) / 2)
      adj <- c(0.5, 0.5 + signy * 2.3)
      if (unit!="")
        text(mean(atx[idx:(idx+1)]), mean(aty[idx:(idx+1)]),
             adj=adj, lab=paste("y [", unit, "]", sep=""),
             cex=cex.axis,
             xpd=TRUE, vfont=font, srt=srt) 

      
      ## z Achse, Beschriftung
      ylab.pos <- pos3D(c(blf[1], blf[2], blf[3] + size[3] / 2), inf)
      lab <- pretty(blf[3] + size[3] * 0:1)
      at <- range(q1[2, ])
      
      at <- at[1] + diff(at) *
        (if (reverse) size[3] - lab + blf[3] else lab - blf[3]) / size[3]
      axis(2, at=at, lab=rep("", length(at)))
      text(ylab.pos[1] - 1.5 * par()$cxy[1], at, lab, srt=90,
           adj=c(0.5, 0), xpd=TRUE, cex=cex.axis, vfont=font)
      if (unit!="") {
        text(ylab.pos[1] - 4.5 * lab.line * par()$cxy[1], ylab.pos[2],
             paste("z [", unit, "]", sep=""),
             adj=c(0.5, 0), xpd=TRUE, cex=cex.axis,  vfont=font, srt=90)
      }     
    }
    
    ## vordere 4 Kanten
    polygon(q1[1,], q1[2,], lty=lty[2], border=col.frame[2])
    
    ## hintere 4 Kanten
    dq2 <- apply((q2-inf[1:2])^2, 2, sum)
    covered <- which(dq2==max(dq2))
    for (i in 1:4) {
      index <- c(i, i %% 4 + 1)
      visible <- 1 + all(covered!=index)
      lines(q2[1, index], q2[2,index], lty=lty[visible], col=col.frame[visible])
    }
    
    ## Kanten, die nach hinten laufen 
    for (i in 1:4) {
      visible <- 1 + (i!=covered)
      lines(c(q1[1,i], q2[1,i]), c(q1[2,i], q2[2,i]), lty=lty[visible],
            col=col.frame[visible])
    }
  } # plot
  
  if (add) return(invisible())
  
  assign(".p3d.range.dist",
         range(sqrt(apply((cbind(qq1, qq2) - sun)^2, 2, sum))), envir=.ENV)

  assign(".p3d.colfct", function(d)
         pmin(length(col),
              pmax(1, (d - .p3d.range.dist[1]) / diff(.p3d.range.dist) *
                   length(col))), envir=.ENV)
  environment(.p3d.colfct) <- NULL

  assign(".p3d.profile", size, envir=.ENV)
  assign(".p3d.inf", inf, envir=.ENV)
  assign(".p3d.sun", sun, envir=.ENV)
  assign(".p3d.col", col, envir=.ENV)
  assign(".p3d.bottomleftfront", bottomleftfront, envir=.ENV)
  assign(".p3d.xlim", xlim, envir=.ENV)
  assign(".p3d.ylim", ylim, envir=.ENV)
  invisible()
}


belichtet <- function(x, radius=1, fix.scale=TRUE, correction=1.2) {
  ## 1.2, da breite des kringels kleiner ist als die
  ##       fuer ihn als Buchstabe reservierte breite
  ##       (Wert fuer eps-Dateien, nicht X11!)
  x <- as.matrix(x)
  x <- x[, order(x[2,], decreasing=TRUE)]
  
  col <- get(".p3d.colfct",env=.ENV)(sqrt(rep(1, nrow(x)) %*%
                             (x - get(".p3d.sun", env=.ENV))^2))

  if (fix.scale) {
    radius = radius / get(".p3d.profile", env=.ENV)[1] * par()$pin[1] /
      par()$cin[1] #par()$csi
   
   points(t(pos3D(x, get(".p3d.inf", env=.ENV))), pch=16,
           col=get(".p3d.col", env=.ENV)[col],
           #col=rainbow(50),
           cex= (2 * correction)*
           ##           * 2, von radius zu durchmesser
           radius * exp(-get(".p3d.inf", env=.ENV)[3] * x[2, ]))
    
  } else {
    x <- t(pos3D(x, get(".p3d.inf", env=.ENV)))
    alpha <- seq(0, 2 * pi, len=15)[-1]
    xcircle <- cos(alpha) * radius
    ycircle <- sin(alpha) * radius
    for (i in 1:nrow(x)) {
      radius = exp(-get(".p3d.inf", env=.ENV)[3] * x[2, ])
      polygon(x = x[i, 1] + xcircle * radius, y = x[i, 2] + ycircle * radius, 
              col=get(".p3d.col", env=.ENV)[col], pch=16)
    }
  }
}


flowpattern <- 
  function(type=c('identical', 'unif', 'independent', 
             'dependent', 'all'),
           length.profile = 200, #200,
           depth = 100, ## of paths
           width.slice = 1,      ## visible one
           delta.x = depth * sqrt(x.var), ## simulated outside length.profil
           delta.y = depth * sqrt(y.var), ## simulated outside width.slice

           lambda.path=0.1,## number of i.i.d ?? paths
           len.x, len.y, # if grid=TRUE alternative to n.path
                                        # overwrites n.path
           grid=FALSE,
           
           x.name='whittle',
           x.var=1,
           x.v.scale=10,
           x.kappa=2,
           x.h.scale=1,#only relevant if dependent is not 'identical'
           ##             or 'independent'
           y.name=x.name,
           y.var=x.var,
           y.v.scale=x.v.scale,
           y.kappa=x.kappa,
           y.h.scale=x.h.scale,
           unif.b=2,
                      
           drop.distr=function(x) x * 80,
           drop.name='whittle',  ## only used if type=all
           drop.scale=1,         ## only used if type=all
           drop.kappa=2,         ## only used if type=all
           drops=1,

           ## risk index calculation
           selected.dist = 2/3,
           front.factor=2, ## should not be changed
           method = NULL, # 'fix.m',
           endpoint.tolerance=0,
           measure=function(x) x^2,
           
           simu.method=if (type %in% c("dependent", "all")) "TBM3",
           drop.simu.method=NULL,
           register = c(1,2),
           old.paths=NULL,
           PrintLevel=RFparameters()$Print,
           wait=FALSE,
           raw=FALSE,
           compress=TRUE,
           max.points = 5500000
           ) {
  raw.xx <- raw.yy <- len.dx <- len.dy <- NULL
  register <- rep(register, len=2)
  if (missing(old.paths)) old.paths <- NULL
  if (is.list(old.paths))
    eval(parse(text=paste(names(old.paths), "<- old.paths$", names(old.paths)))) 
  type <- match.arg(type)
  
  if (any(unlist(RFparameters()[c("TBM2.simulinefactor","TBM2.simulinefactor")])
          != 0.0) && (is.null(simu.method) || simu.method="TBM3" ||
                  simu.method=="TBM2")) {
    stop("simulinefactor must be zero for this function; see RFparameters for this global parameter")
  }
  stopifnot(!type %in% c("dependent", "all") ||
            length.profile == round(length.profile),
            depth==round(depth) )
 
# library(SoPhy, lib="~/TMP"); source("~/R/SOPHY/SoPhy/R/3dplot.R"); example(flowpattern);
  
  dependent.path <- function(startx, starty, depth, model, grid, reg=0) {

#    print(startx)
#    print(starty)
#    str(depth)
#    str(model)
#    str(grid)
    
    InitGaussRF(x=c(-delta.x, length.profile + delta.x),
                y=c(-delta.y, width.slice + delta.y),
                z=c(1, depth),
                grid=FALSE, model=model, reg=reg, method=simu.method)

#str(GetRegisterInfo(reg), vec=20)
#str(RFparameters())
    
    mem <- GetRegisterInfo(reg)$method[[1]]$mem
    if (grid) {
      t(matrix(GaussRF(x=startx,
                       y=starty,
                       z=1:depth,
                       grid=grid, model=model, reg=reg, method=simu.method,
                       TBM2.linesimufactor=0.0, TBM2.linesimustep=0.0,
                       TBM3.linesimufactor=0.0, TBM3.linesimustep=0.0,
                       TBM.points=length(mem$l),
                       TBM.center=if (!is.null(mem$aniso)) solve(mem$aniso,
                         mem$center) + 0.0 else 0.0
                       ),
               ## +0.0: force double value for TBM.center...
               ncol=depth))
    } else {
      matrix(GaussRF(x=rep(startx, each=depth),
                     y=rep(starty, each=depth),
                     z=rep(1:depth, length(startx)),
                     grid=grid, model=model, method=simu.method, reg=reg,
                     TBM2.linesimufactor=0.0, TBM2.linesimustep=0.0,
                     TBM3.linesimufactor=0.0, TBM3.linesimustep=0.0,
                     TBM.points=length(mem$l),
                     TBM.center=if (!is.null(mem$aniso)) solve(mem$aniso,
                       mem$center) + 0.0 else 0.0
                     ),
             nrow=depth)
    }
  }
  
  if (!is.list(old.paths)) {
    x.model <-
      switch(type,
             identical = list(
               list(model=x.name, variance=x.var,
                    scale=x.v.scale, kappa=x.kappa)
               ),
             independent = list(
               list(model=x.name, variance=x.var,
                    scale=x.v.scale, kappa=x.kappa),
               ),
             dependent = list(
               list(model=x.name, variance=x.var, kappa=x.kappa,
                    aniso=diag(1.0 / c(x.h.scale, x.h.scale, x.v.scale))),
               ),
             all = list(
               list(model=x.name, variance=x.var, kappa=x.kappa,
                    aniso=diag(1.0 / c(x.h.scale, x.h.scale, x.v.scale))),
               )
             )
    y.model <-
      switch(type,
             identical = list(
               list(model=y.name, variance=y.var,
                    scale=y.v.scale, kappa=y.kappa)
               ),
             independent = list(
               list(model=y.name, variance=y.var,
                    scale=y.v.scale, kappa=y.kappa),
               ),
             dependent = list(
               list(model=y.name, variance=y.var,
                    aniso=diag(1.0 / c(y.h.scale,
                      y.h.scale, y.v.scale)), kappa=y.kappa),
               ),
             all = list(
               list(model=y.name, variance=y.var,
                    aniso=diag(1.0 / c(y.h.scale,
                      y.h.scale, y.v.scale)), kappa=y.kappa),
               )
             )

    rp <- rpois(1, lambda.path *
                (length.profile + 2*delta.x) * (width.slice + 2*delta.y))
    if (PrintLevel>2) {
      if (rp>100000) {
        print(c(lambda.path *
                (length.profile + 2*delta.x) * (width.slice + 2*delta.y),
                lambda.path,
                length.profile,
                delta.x,
                width.slice,
                delta.y, rp)
              )
      } else print(rp)
    }
    
    sx <- runif(rp, -delta.x, length.profile + delta.x)
    sy <- runif(rp, -delta.y, width.slice + delta.y)
    if (grid) {
      if (missing(len.x) || missing(len.y)) {
        PointsPerSegm <- sqrt(rp / (length.profile * width.slice))
        len.x <- round(length.profile * PointsPerSegm)
        len.y <- round(width.slice * PointsPerSegm)
      }
      startx <- seq(0, length.profile, len=len.x)
      
      dx <- floor(delta.x / length.profile * len.x)
      dx <- if (dx==0) NULL else 1:dx * startx[2]
      startx <- c(-1 * rev(dx), startx, length.profile + dx)# -1*, da dx u.U. NULL
      len.x <- length(startx)
      starty <- seq(0, width.slice, len=len.y)
      dy <- floor(delta.y / length.profile * len.y)
      dy <- if (dy==0) NULL else 1:dy * starty[2]
      starty <- c(-1 * rev(dy), starty, width.slice + dy)
      len.y <- length(starty)
      rp <- len.x * len.y     
    } else {
      ALLSTARTX <- sort(sx)
      ALLSTARTY <- sy
      idxlength <- if (compress && any(type %in% c("dependent", "all")))
          as.integer(max.points / depth) else Inf

      if (PrintLevel>2)
        cat("idxlength=", idxlength, "    true length=", length(ALLSTARTX), "\n")
      idxmax <- min(idxlength, length(ALLSTARTX))
      startx <- ALLSTARTX[1:idxmax]
      starty <- ALLSTARTY[1:idxmax]
      seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
      YY <- STARTX <- STARTY <- NULL
    }
    
    sx <- sy <- NULL
    repeat {
      yy <- switch(type,
                   identical = GaussRF(x=1:depth, grid=TRUE, model=y.model,
                     n=rp, reg=register[1], method=simu.method, pch=""),
                   unif=matrix(0, nrow=depth, ncol=rp),
                   independent = GaussRF(x=1:depth, grid=TRUE, model=y.model,
                     n=rp, reg=register[1], method=simu.method, pch=""
                     ),
                   dependent=dependent.path(startx, starty, depth, y.model,
                     grid=grid, reg=register[1]),
                   all = dependent.path(startx, starty, depth, y.model,
                     grid=grid, reg=register[1]) 
                   )      
      stopifnot(is.numeric(yy))
      ## y path is "horizontal" !
      if (PrintLevel>2) {
        print(quantile(starty, c(0, 0.005, 0.01, 0.02, 0.05, 0.95,
                                 0.98, 0.99, 0.995, 1)))
        if (PrintLevel>3)
          print(quantile(apply(yy, 2, cumsum), 
               c(0, 0.005, 0.01, 0.02, 0.05, 0.95, 0.98, 0.99, 0.995, 1)))
      }
      if (raw && grid) raw.yy <- array(t(yy), dim=c(len.x, len.y, depth))
      
      if (PrintLevel>2) cat("size of yy ", length(yy),"\n")
     
      yy <- (if (grid) rep(starty, ea=len.x) else starty) +
        t(.Call("matrixCumsum", yy, 2, PACKAGE="SoPhy"))

      if (compress && !grid) {
        idx <- integer(nrow(yy))
        if (PrintLevel>1) cat("compressing y: ")
        .C("anyinside", yy, as.integer(nrow(yy)), as.integer(ncol(yy)),
           as.double(0), as.double(width.slice), idx,
           PACKAGE="SoPhy", DUP=FALSE)
        idx <- as.logical(idx)
        if (PrintLevel>1) cat(length(idx), "->", sum(idx), "\n")
        YY <- rbind(YY, yy[idx, ])
        yy <- NULL
        STARTX <- c(STARTX, startx[idx])
        STARTY <- c(STARTY, starty[idx])
        ALLSTARTX <- ALLSTARTX[-1:-idxmax]
        ALLSTARTY <- ALLSTARTY[-1:-idxmax]
        gc()
        if (length(ALLSTARTX)==0) {
          startx <- STARTX
          starty <- STARTY
          rp <- length(startx)
          yy <- YY
          remove("STARTX", "STARTY", "YY")
          break
        } else {
          idxmax <- min(idxlength, length(ALLSTARTX))
          startx <- ALLSTARTX[1:idxmax]
          starty <- ALLSTARTY[1:idxmax]
          assign(".Random.seed", seed, envir=.GlobalEnv)
        }
      } else break
    }
    gc()
    
    xx <-
      switch(type, 
             identical = (matrix(GaussRF(x=1:depth, grid=TRUE, model=x.model,
                                         method=simu.method),
                                      nrow=depth, ncol=rp)),
             unif=matrix(sqrt(runif(rp, 1, unif.b)^2 - 1), nrow=depth, 
                         ncol=rp, byrow=TRUE) * 
                  (2 * rbinom(depth * rp, 1, 0.5) - 1),
             independent = GaussRF(x=1:depth, grid=TRUE, model=x.model, n=rp,
               method=simu.method),
             dependent=dependent.path(startx, starty, depth, x.model, grid=grid),
             all = dependent.path(startx, starty, depth, x.model, grid=grid)
             ) 
    if (raw && grid) raw.xx <- array(t(xx), dim=c(len.x, len.y, depth))
    if (PrintLevel>3) {
      screen(1)
      image(1:ncol(xx), -nrow(xx):-1, t(xx)[, nrow(xx):1, drop=FALSE],
            cex=1, cex.axis=1, cex.lab=1)
    }
    
    ## x path is "horizontal" !
    xx <- (if (grid) rep(startx, len.y) else startx) +
      t(.Call("matrixCumsum", xx, 2, PACKAGE="SoPhy"))
    gc()
    if (compress && !grid) {
      idx <- integer(nrow(xx))
      .C("anyinside", xx, as.integer(nrow(xx)), as.integer(ncol(xx)),
         as.double(0), as.double(length.profile), idx,
         PACKAGE="SoPhy", DUP=FALSE)
      idx <- as.logical(idx)
      if (PrintLevel>1) cat("compressing x:",length(idx), "->", sum(idx), "\n")
      xx <- xx[idx, , drop=FALSE]
      yy <- yy[idx, , drop=FALSE]
      rp <- sum(idx)
      startx <- startx[idx]
      starty <- starty[idx]
    }
    if (PrintLevel>3) {
      screen(2)
      image(1:nrow(xx), -ncol(xx):-1, xx[, ncol(xx):1, drop=FALSE] - startx,
            cex=1, cex.axis=1, cex.lab=1)
    }
    gc()
  } # !is.list(old.paths)
     
  drop.model <-
    if (type=="all") list(list(model=drop.name, variance=1, scale=drop.scale,
          kappa=drop.kappa))
    else list(list(model="nugget", variance=1, scale=1))
  if (!is.list(old.paths)) {
    if (PrintLevel>2) cat("drops\n")
    LEN <- pnorm(matrix(GaussRF(startx,
                                 starty,
                                 grid=grid, model=drop.model, n=drops,
                                 reg=register[2], method=drop.simu.method),
                         ncol=drops))
  }
  len <- apply(drop.distr(LEN), 1, max)
 
  nc <- ncol(xx)
  if (PrintLevel>2) cat("matrixcumsum\n");
  IDX <- .Call("matrixCumsum",
               sqrt((xx[, -1, drop=FALSE] - xx[, -nc, drop=FALSE])^2 +
                    (yy[, -1, drop=FALSE] - yy[, -nc, drop=FALSE])^2 + 1),
               1, PACKAGE="SoPhy") <= len
  ## entweder liegt der Punkt in der Schicht, oder aber
  ## bei zwei aufeinanderfolgenden Tiefen liegt ein Punkt vor, der andere
  ## hinter der Schicht.
  idx <- t(  (yy > 0 | cbind(yy[, -1, drop=FALSE], -Inf) > width.slice)
           & (yy < width.slice | cbind(yy[, -1, drop=FALSE], Inf) < 0)
           & cbind(TRUE, IDX) # travelling distance not exceeded
           & xx >= 1 &
            xx < length.profile+1)
  i.x <- as.integer(t(xx)[idx])
  i.d <- 1 + (which(idx) - 1) %% depth
    
  picture <- matrix(FALSE, nrow=depth, ncol=length.profile)
  picture[cbind(i.d, i.x)] <- TRUE
    
  freq <- apply(picture, 1, sum)
  dist <- 1:nrow(picture)
  
  if (PrintLevel>2) {
    if (PrintLevel>3) cat("\n max.freq=", max(freq))
    screen(3)
    plot(i.x, -i.d, ylim=c(-depth, -1), pch=".",
         cex=1, cex.axis=1, cex.lab=1)
    if (wait) readline()
  }
  
  c.r <- list()
  if (length(method)>0) for (m in 1:length(method)) {
    if (PrintLevel>2) cat("\n method=", method[m], "\n", sep="")
    c.r[[m]] <- risk.index(cbind(dist, freq),
                           selected.dist=selected.dist,
                           selected.rate=NULL,
                           front.factor=front.factor,
                           PrintLevel=PrintLevel,
                           endpoint.tolerance=endpoint.tolerance,
                           measure = measure,
                           method=method[m])
  }
  
  nn <- names(formals())
  if (missing(len.x) || missing(len.y)) len.x <- len.y <- NULL

  return(c(c.r,
           list(intermediate=list(xx=xx, yy=yy, LEN=LEN,
                  raw.xx=raw.xx, raw.yy=raw.yy, x=if (grid) startx,
                  y=if (grid) starty)),
           list(dist=dist, freq=freq, i.x=i.x, i.d=i.d),
	   list(input=eval(parse(text=paste("list(",
                                   paste(paste(nn, nn, sep="="),
                                   collapse=", "), ")"))))
	   ))
}


plotFlow3d <- function(paths, horizons=c("no", "absorbing", "breakthrough"),
                       drop.distr, n.balls=1, pointradius=1,
                       dev=1, ps="3d.dye.pattern", ps.background=FALSE,
                       profileheight=4, unit="cm", unit.scale=1, inf, sun,
                       rl = function(x) readline(paste(x, ": press return")),
                       low.resolution = TRUE,
                       col=grey(pmin(1, pmax(0, seq(0.95,0,-0.001))))
                       ){  
  if (missing(drop.distr)) drop.distr <- paths$input$drop.distr
  horizons <- match.arg(horizons)
  xx <- paths$intermediate$xx
  yy <- paths$intermediate$yy
  LEN <-  paths$intermediate$LEN
  paths$intermediate <- NULL
  x1 <- xx[,1]
  y1 <- yy[,1]
    
  hor <- rep(1, (ncol(xx)-1))
  if (horizons!="no") {
    # warning("only temptative approach")
    hor[50:60] <- runif(11, -1, 1) * 50 ## oder fix to 50
  }
  if (horizons=="absorbing") {
    xx <- matrix(rep(t((xx[,-1] - xx[, -ncol(xx)]) / n.balls) * hor,
                     each=n.balls), nrow=nrow(xx), byrow=TRUE)
    yy <- matrix(rep(t((yy[,-1] - yy[, -ncol(yy)]) / n.balls) * hor,
                     each=n.balls), nrow=nrow(yy), byrow=TRUE)
    hor <- 1
  } else {
    xx <- matrix(rep(t((xx[,-1] - xx[, -ncol(xx)]) / n.balls),
                     each=n.balls), nrow=nrow(xx), byrow=TRUE)
    yy <- matrix(rep(t((yy[,-1] - yy[, -ncol(yy)]) / n.balls),
                     each=n.balls), nrow=nrow(yy), byrow=TRUE)
    hor <- matrix(rep(hor, each=n.balls), ncol=ncol(xx),
                  nrow=nrow(xx), byrow=TRUE)
  }
    
  ## Reihenfolge beibehalten!!
  path.length <- .Call("matrixCumsum", sqrt(xx^2 + yy^2 + n.balls^(-2)), 1,
                       PACKAGE="SoPhy")  
  xx <- x1 + .Call("matrixCumsum", cbind(0, xx * hor), 1, PACKAGE="SoPhy")
  yy <- y1 + .Call("matrixCumsum", cbind(0, yy * hor), 1, PACKAGE="SoPhy")
  hor <- NULL
      
  within.slice <- (xx>=0 & xx<paths$inp$length.profile &
                   yy>=0 & yy<paths$inp$width.slice)
      
  ## recalculated only because horizon can be != "no"
  ## and for efficiency as drop.distr might be a vector
  IDX <-(cbind(TRUE, path.length <= apply(drop.distr(LEN),
                       1, max)) & within.slice)
  x <- rbind(as.vector(xx[IDX]), as.vector(yy[IDX]),
             as.vector(matrix((ncol(xx):1)/n.balls, nrow=nrow(xx),
                              ncol=ncol(xx), byrow=TRUE)[IDX])) * unit.scale

  qq <- function(plot)
      quader(size=unit.scale * c(paths$inp$length.profile, paths$inp$width.slice,
               paths$inp$depth),
             inf=inf, sun=sun,
             cex.axis=1.8, # 1.2
             col=col,
             plot=plot,
             unit=unit
             )

  qq(FALSE)
  if (!is.null(dev)) {
    mai <- c(rep(0.6, 2) + is.logical(dev) * 0.6 , 0.1, 0.1)
    totalheight <- profileheight + mai[1] + mai[3]
    totalwidth <- profileheight / diff(get(".p3d.ylim", env=.ENV)) *
      diff(get(".p3d.xlim", env=.ENV)) + mai[2] + mai[4]
    Dev(TRUE, dev, ps=ps, height=totalheight, width=totalwidth)
    par(mai=mai, lwd=1 + is.logical(dev)) 
    qq(TRUE)
    belichtet(x, radius = pointradius)
    quader(add=TRUE, lty=c(0,1), col=c("transparent", "black"))
    if (is.numeric(dev) && dev!=1) rl(ps)
    Dev(FALSE)
  }
  
  ## tiff.file <- ps
  if (is.logical(dev) && dev && low.resolution && .Platform$OS.type=="unix") {
    exts <- c("ps", "eps")
    if (any(l <- (splt <- rev(strsplit(ps,"\\.")[[1]]))[1] == exts) &&
        length(splt)>1) {
      ps <- paste(rev(splt[-1]), collapse=".")
      ext <- exts[l]
    } else ext <- "eps"
    tiff.file <- paste(ps, "...dummy.tiff",  sep="")
    txt <- paste("nice -n 19 convert ", ps, ".", ext, " ", tiff.file, " && ",
                 "nice -n 19 convert ", tiff.file, " eps:", ps, ".", ext, " && ",
                 "rm ", tiff.file, if (ps.background) "&" else "",
                 sep="")
    system(txt)
  }
  invisible(t(x))
}



plotFlow2d <- function(coord, 
                       pointradius=1, slice=2 * pointradius, full.size=TRUE,
                       Profiles=1, dev=1, ps="",  
                       height=4,  unit="cm", cex=2, correction=1.2, col=1,
                       rl = function(x) readline(paste(x, ": press return"))
)
{
  stopifnot(Profiles>=1)
  size <- get(".p3d.profile", envir=.ENV)
  blf <- get(".p3d.bottomleftfront", env=.ENV)
  m <- blf[2] + if (Profiles==1) size[2]/2  else seq(0, size[2], len=Profiles)
  Pro <- m - 0.5 * slice
  dig <- nchar(format(Profiles)) - 1 # -1 ist mir unlogisch aber
  ## notwendig, sonst eine fuehrende Null zu viel
  a <- 0.25 / slice^2
  mai <- c(rep(0.8, 2) + 0.2 * is.logical(dev), 0.1, 0.3)
  totalheight <- height + mai[1] + mai[3]
  totalwidth <- height / size[3] * size[1] + mai[2] + mai[4]
  
  for (i in 1:Profiles) {
    idx <- coord[, 2] >= Pro[i] & coord[, 2] <= Pro[i] + slice
    psname <- paste(ps, formatC(i, flag="0", digits=dig), sep="")
    Dev(TRUE, dev, ps=psname, height=totalheight, width=totalwidth)
    par(mai=mai)
    plot(Inf, Inf,
         xlim=(c(0, size[1]) + blf[1]), ylim=(c(0, size[3]) + blf[3]),
         xlab=if (unit=="") "" else paste("x [", unit, "]", sep=""),
         ylab=if (unit=="") "" else paste("z [", unit, "]", sep=""),
         xaxs="i", yaxs="i", axes=FALSE, frame=TRUE,
         cex.axis=cex, cex.lab=cex)
    lab <- pretty(c(0, size[3]))
    reverse <- TRUE
    at <- if (reverse) size[3] - lab else lab
    axis(1, cex.axis=cex)
    axis(2, at=at , lab=lab, cex.axis=cex)

    radius = (2 * correction) * pointradius / size[1] *
      par()$pin[1] / par()$cin[1]
    points(coord[idx, 1], coord[idx, 3], pch=16,
           cex = (if (full.size) radius else radius *
                  sqrt(1 - a * (coord[idx, 2] - m[i])^2)),
           col=col
           )
    Dev(FALSE)
    if (is.numeric(dev) && dev!=1 && i<Profiles) rl(psname)
  }
  invisible(NULL)
}

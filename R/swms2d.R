.ENV <- environment()
.ENV <- .GlobalEnv

.onLoad <- function (lib, pkg) {
 # library("RandomFields", lib=lib)
#  library.dynam("SoPhy", pkg, lib) 

  if (file.exists("/home/schlather/bef/x")) {
    ## to do list -- since my surname is rare, the message should 
    ## appear only on computers I have a login
    cat("To-Do List\n==========\n")

    print("package grid\n s4: useful for sophy print functions, etc.")
    print("wurzeln, parameter wie depth zuschlag funktionieren nicht mehr")
    print("convert create.stones into c programme")
    print("convert create.roots into c programme")
    print("improve cut[][]-area!! -- add lines at the very end that check how much the cut area can be reduced.")
    print("rf.link not tested yet -- improve standard definition of rf.link!")
  }
}


##AtmInF=AtmInF: is.null(atmosphere)
##MaxAL = nrow(atmosphere)
##lChem = (!is.null(cBound) && (!is.na(cBound[1])))
##DrainF= (length(d$KElDr)!=0) && !is.na(d$KElDr[1])
##SeepF = (length(d$NP)!=0) && !is.na(d$NP[1])
swms2d <- function(d,
                   max.iteration = 1e+05,
                   iter.print = max.iteration,
                   ShortF = TRUE, #information printed only at selected instances
                   message = NULL,
                   breakpoint = 1e+10,
                   intermediate.result = NULL
#                   NumNPD=NumNP,
#                   NumElD=NumEl,
#                   NumBPD=NumBP,
#                   NMatD =NMat
                   ){
  h <- NULL
  if (is.null(d$nCodeM)) {
    if (d$H1$type=="Start") {
      h <- d
      d <- create.waterflow(d)
      if (is.character(d)) stop(d)
    } else stop("neither a swms2d input list nor a profile definition")
  }
  FluxF <- TRUE  ## hence output of Q, v, conc!
  NSeepD <- 2    #
  NumSPD <- 50
  NDrD <- 2
  NElDrD <- 8
  NumKD <- 6   ## always >= 5 !!
  NObsD <- 4
  MBandD <- 20
  NTabD <- 100
  MNorth <- 4

  len <- function(x, l, m, value=0) {
    ##  print(match.call())
    if (is.list(x)) {
      stop("??")
      m <- max(unlist(lapply(x, length)))
      x <- t(sapply(x,function(e) c(e, rep(NA,m))[1:m]))
    } else if (is.matrix(x)) {
      m <- matrix(value, nrow=l, ncol=m)
      m[1:nrow(x),1:ncol(x)] <- x
      return(m)
    }
    return(c(as.vector(x), rep(value, l))[1:l])
  }

  ###################### printing information  
  TPrint <- d$TPrint
  if (any(diff(TPrint)<=0)) return("TPrint not an ordered sequence")
  if ((NObs <- length(d$Node)) > NObsD) return("too many observation points")

  ######################## nodal information
  ###### row for each node
  ## if only 11 columns then assumed first column is missing 
  if ((NumNP <- nrow(d$nCodeM)) < 2) return("too few points")
  ##  if (NumNP > NumNPD) return("too many points")
  if (ncol(d$nCodeM)==11) d$nCodeM <- cbind(NumNP, d$nCodeM) else
  if (ncol(d$nCodeM)!=12) return("ncol(nCodeM)!=12")

  D1double <- len(as.matrix(d$nCodeM[,c(3:7,9:12)]), NumNP, 9)
  nCodeM <- len(as.matrix(d$nCodeM[,-c(3:7,9:12)]), NumNP, 3)
  if (any(nCodeM!=as.integer(nCodeM)))
    return("integers expected in columns 1, 2, and 8 (of 12 columns)")
  
  ######################## boundary information

  KXB <- d$nCodeM[d$nCodeM[,2] != 0, 1]
  NumBP <- length(KXB)
  ##if (NumBP > NumBPD) {
  ##  return(paste("too many boundary points (len(KXB)=", length(KXB),
  ##               "> NumBPD=", NumBP,")"))
  ##}
  if (length(d$Width) == 1) d$Width <- rep(d$Width, NumBP) else
  if (length(d$Width) != NumBP) return("length(Width) != NumBP")


  ######################## element information
  ###### row for each element
  ## if only 8 columns then assumed first column is missing 
  if ((NumEl <- nrow(d$KXR)) < 2) return("too few elements")
  ## if (NumEl > NumElD) return("too many elements")
  if (ncol(d$KXR)==8) d$KXR <- cbind(1:NumEl, d$KXR) else
  if (ncol(d$KXR)!=9) return("ncol(KXR)!=9")
  ConAX <- len(as.matrix(d$KXR[,6:8]), NumEl, 3)
  KXR <- len(as.matrix(d$KXR[,-6:-8]), NumEl, 6) 
  
  ######################## material
  ###### row for each material
  NPar <- 9  ## only van Genuchten model used
  if (is.vector(d$Par)) d$Par <- t(d$Par)
  if (ncol(d$Par)!=NPar)return("only van Genuchten model allowed (9 parameters)")
  NMat <- nrow(d$Par)
  ## if (NMat > NMatD) return("too many materials")
  ## to do: check grid definition, material numbers!  

  ######################### atmospherical input
  ####### row for each time step
  root <- rep(NA, 7)
  
  if (AtmInF <- !is.null(d$atmosphere) && !is.na(d$atmosphere[1])) {
    if (ncol(d$atmosphere)!=10) return("ncol(atmosphere)!=10")
    MaxAL <- nrow(d$atmosphere) ## no constraints since data are
    ##                             read in step by step within TIME2.f
    tMax <- d$atmosphere[nrow(d$atmosphere),1]
    if (tMax < (TPrintMax <- TPrint[length(TPrint)])) {
      warning(paste(tMax,"= tMax < max(TPrint)=",TPrintMax,
                    " -- TPrint is shortened", sep=""))
      TPrint <- c(TPrint[TPrint < tMax], tMax)
    } else {
      ## tmax is automatically shortened !!!
      tMax <- TPrintMax
    }
    ######################## root uptake/sink, only if atmospherical information
    if (SinkF <- d$SinkF) {
      if ((length(d$POptm)==0) || is.na(d$POptm[1])) 
         return("SinkF, but invalid POptm")
      ## P0, P2H, P2L, P3, r2H, r2L
      if (is.vector(d$root)) {
        if (length(d$root)!=6) return("length(root)!=6")
        if (length(d$POptm)!=NMat) return("length(POptm)!= # of materials")
        root <- cbind(t(matrix(d$root,nrow=6,ncol=NumNP)),
                      d$POptm[d$nCodeM[,8]])
      } else {
        if (!is.matrix(d$root) || (ncol(d$root)!=6) || (nrow(d$root)!=NumNP))
          return("root not a matrix of size NumNP x 6")
        if (length(d$POptm)!=nrow(d$nCodeM)) return("length(POptm)!= # of nodes")
        root <- cbind(d$root, as.vector(d$POptm))
      }
      ## siehe INPUT2.f !
      root[,1:4] <- -abs(root[,1:4])
    } 
  } else {
    ## no atmospheric input
    SinkF <- d$qGWLF <- FALSE
    d$atmosphere <- tMax <- d$GWL0L <- d$Aqh <- d$Bqh <-  d$hCritS <- NA
    d$tInit <- 0 ## the default value in SWMWS_2d, necessary to check
    ##              whether the 1st TPrint value is far enough away from
    ##              tInit -- otherwise swms2d crashes!
    MaxAL <- 0
  }
  if (TPrint[1] < d$tInit + d$dt) {
    ## unclear which condition is indeed necessary; seems that
    ## this one works fine -- if 1st value is too small the system
    ## crashes
    init.TP <- d$tInit + d$dt
    TPrint <- c(init.TP, TPrint[TPrint > init.TP])
    warning(paste("TPrint[1] is too small. -- increased to", init.TP))    
  }
  
  ######################### chemical input
  NChPar <- 9
  if (lChem <- (!is.null(d$cBound) && (!is.na(d$cBound[1])))) {
   if (length(d$cBound) != 6) return("length of cBound incorrect")
   if (length(d$KodCB) != NumBP) return("length(KodCB) != NumBP")
   if (is.vector(d$ChPar)) d$ChPar <- t(d$ChPar)   
   if (nrow(d$ChPar) != NMat) return("nrow(ChPar) != NMat")
   if (ncol(d$ChPar) != NChPar) return("ncol(ChPar) != 9")       
   } else {
     ## no need to give the parameters
     d$Epsi <- d$lUpW <-d$lArtD <- d$PeCr <- d$ChPar <- d$KodCB <- d$tPulse <- NA
     cBound <- rep(NA, 6)
  }

  ######################### drainage 
  if (DrainF <- (length(d$KElDr)!=0) && !is.na(d$KElDr[1])) {    
    KElDr <- matrix(nrow=NDrD, ncol=NElDrD) 
    if (is.list(d$KElDr)) {
      NDr <- length(d$KElDr)
      NED <- integer(NDr)
      for (i in 1:NDr) {
        if ( (l <- NED[i] <- length(d$KElDr[[i]])) > NElDrD)
          return(paste("KElDr[[", i, "]] has too many elements"))
        KElDr[i, 1:l] <- d$KElDr[[i]]
      }
    } else {
      ## information given in form of a matrix
      if (!is.matrix(d$KElDr)) return("KElDr must be a list or a matrix")
      NDr <- nrow(d$KElDr)
      NED <- integer(NDr)
      if ((l <- ncol(d$KElDr)) > NElDrD) return(paste("KElDr: too many columns"))
      for (i in 1:NDr) {
        NED[i] <- sum(is.finite(d$KElDr[i, 1:l]))
        if(any(!is.finite(KElDr[i, 1:l] <- d$KElDr[i, 1:l])))
          return("KElDr: non finite elements are not allowed")
      }
    }
    if (NDr > NDrD) return("NDr too large")
    if (length(ND) != NDr) return("length(ND) != NDr")
    if ((NDr==1) && is.vector(d$EfDim)) d$EfDim <- t(d$EfDim)
    if (!is.matrix(d$EfDim) || (ncol(d$EfDim)!=2) || (nrow(d$EFDim!=NDr)))
      return("EfDim must be a matrix of 2 columns and NDr rows")
  } else {
    ## no drain
    NDr <- DrCorr <- ND <- NED <- EfDim <- KElDr <- NA    
  }

  ######################### seepage 
  if (SeepF <- (length(d$NP)!=0) && !is.na(d$NP[1])) {
    NP <- matrix(nrow=NSeepD, ncol=NumSPD)
    if (is.list(d$NP)) {
      NSeep <- length(d$NP)
      NSP <- integer(NSeep)
      for (i in 1:NSeep) {
        if ((l <- NSP[i] <- length(d$NP[[i]])) > NumSPD)
          return(paste("NP[[", i, "]] has too many elements"))
        NP[i, 1:l] <- d$NP[[i]]
      }
    } else {
      if (!is.matrix(d$NP)) return("NP must be a list or a matrix")
      NSeep <- nrow(d$NP)
      NSP <- integer(NSeep)
      if ((l <- ncol(d$NP)) > NumSPD) return(paste("NP has too many columns"))
      for (i in 1:NSeep) {
        NSP[i] <- sum(is.finite(d$NP[i, 1:l]))
        if (any(!is.finite(d$NP[i, 1:NSP[i]])))
          return("NP: non finite elements are not allowed")
        NP[i, 1:l] <- d$NP[i, 1:l]
      }
   }
    if (NSeep>NSeepD) return("NSeep too large")
  } else {
    ## no seepage
     NSP <- NP <- NSeep <- NA
  }
  if (length(d$hTab)!=2) return("length(hTab)!=2")
  if (length(d$dtMinMax)!=2) return("length(dtMinMax)!=2")
  if (length(d$DMul)!=2) return("length(DMul)!=2")
  
  NInt <- 25
  NDbl <- 27
  
  IntVec <-
    c(Kat=d$Kat,     MaxIt=d$MaxIt, lWat=d$lWat,  lChem=lChem,   ShortF=ShortF,
      FluxF=FluxF,   AtmInF=AtmInF, SeepF=SeepF,  FreeD=d$FreeD, DrainF=DrainF,
      NLay=d$NLay,   NSeep=NSeep,  NDr=NDr,       NObs=NObs,     SinkF=SinkF,
      qGWLF=d$qGWLF, lUpW=d$lUpW,  lArtD=d$lArtD, iter.print=iter.print,
                                                                 start = TRUE,
      maxLoop = max.iteration, breakpoint=breakpoint, brk.tprint=NA, ii=0,
                                                                 iiprint = NA
      )
  names.IntVec <- names(IntVec)
  
  DblVec <-
    c(TolTh=d$TolTh,        TolH=d$TolH,  hTab=d$hTab,  dt=d$dt,
      dtMinMax=d$dtMinMax,  DMul=d$DMul,  DrCorr=d$DrCorr,
      rLen=d$rLen,  GWL0L=d$GWL0L,  Aqh=d$Aqh,  Bqh=d$Bqh, tInit=d$tInit,
      hCritS=d$hCritS, tMax=tMax,  Epsi=d$Epsi,  PeCr=d$PeCr,
      cBound=d$cBound, # 6 variables
      tPulse=d$tPulse, time=NA)
  names.DblVec <- names(DblVec)

  if (length(IntVec)!=NInt || length(DblVec)!=NDbl) return("unexpected nulls")
  
  MPL <- length(TPrint)
  hQThFlC <- array(0, dim=c(6, NumNP, MPL+2))

  TlSolObs.ncol <- 10 + 3 * (NumKD + 1) + 3 + 2 * NumKD + 3 * NObs
  TlSolObs.nrow <- max(1 - d$lWat, if (ShortF) MPL + 1 else max.iteration)
  TlSolObs <- matrix(nrow=TlSolObs.nrow, ncol=TlSolObs.ncol)
  ## note interpretation of the variables depends on lWat, see TLInf
  
  atmOut <-  matrix(nrow=MaxAL, ncol=9)

  boundary.n <- 2 + if(!d$lWat && !lChem) max.iteration else MPL #!!
  boundary <- array(0, dim=c(NumBP, 11,  boundary.n))

  balance.ncol <- 5 + (4 + 2 * lChem) * (d$NLay + 1)
  balance <- matrix(0, nrow=boundary.n, ncol=balance.ncol)

  if (AtmInF && (!is.matrix(d$atmosphere) || (ncol(d$atmosphere)!=10)))
      return("`atmosphere' not compatible with AtmInF")

  warn.orig <- options()$warn
  options(warn=-1)
  storage.mode(IntVec) <- "integer"
  storage.mode(DblVec) <- "double"
  storage.mode(NMat) <- "integer"
  ## first t() then len(), since d$Par is not matrix
  ## but will be transformed to a matrix by t()
  Par <- as.double(len(t(d$Par), 10, NMat))
  storage.mode(MPL) <- "integer"
  ## see INPUT2.f, time point added!
  TPrint <- as.double(len(TPrint, MPL + 1))
  NSP <- as.integer(len(NSP, NSeepD))
  storage.mode(NP) <- "integer"
  ND <- as.integer(len(d$ND, NDrD))
  NED <- as.integer(len(NED, NDrD))
  EfDim <- as.double(len(d$EfDim, NDrD, 2))
  storage.mode(KElDr) <- "integer"
  storage.mode(NumNP) <- "integer"
  storage.mode(NumEl) <- "integer"
  storage.mode(NumBP) <- "integer"
  storage.mode(nCodeM) <- "integer"
  storage.mode(D1double) <- "double"
  storage.mode(KXR) <- "integer"
  storage.mode(ConAX) <- "double"
  storage.mode(KXB) <- "integer"
  Width <- as.double(d$Width)
  Node <- as.integer(d$Node)
  storage.mode(MaxAL) <- "integer"
  ChPar <- as.double(len(t(d$ChPar), 10, NMat))
  KodCB <- as.integer(d$KodCB)# no len(), see SWMS_2D !
  atmosphere <- as.double(d$atmosphere)
  storage.mode(root) <- "double"
  storage.mode(TlSolObs.nrow) <- "integer"
  storage.mode(TlSolObs.ncol) <- "integer"
  storage.mode(boundary.n) <- "integer"
  storage.mode(balance.ncol) <- "integer"
  storage.mode(hQThFlC) <- "double"
  storage.mode(TlSolObs) <- "double"
  storage.mode(atmOut) <- "double"
  storage.mode(balance) <- "double"
  storage.mode(boundary) <- "double"
  error <- integer(1)
#  storage.mode() <- ""
  options(warn=warn.orig)
  

  n.intpck <- NumEl * 5 + NumNP * (MBandD + 5) + 21 + 15 + 1
  intpck <- integer(n.intpck)
  n.dblpck <- (NumNP * (2 * MBandD + 2 * MNorth + 31) + NMat * 24 +
               NTabD * (1 + 3 * NMat) + NumEl * 11 + NumKD * 6 +
               NumBP * 3 + 6 + MNorth + 51
               + 1)
  dblpck <- double(n.dblpck)
  
  repeat {
    .Fortran("swms2d",
             ##input
             IntVec=IntVec,
             DblVec=DblVec,
             NMat=NMat,
             Par=Par,
             MPL=MPL,
             TPrint=TPrint, 
             NSP=NSP,
             NP=NP,
             ND=ND,
             NED=NED,
             EfDim=EfDim,
             KElDr=KElDr,
             NumNP = NumNP, 
             NumEl = NumEl,   
             NumBP = NumBP,   
             nCodeM= nCodeM,
             D1double=D1double,
             KXR=KXR,
             ConAX=ConAX,
             KXB=KXB,
             Width= Width,
             Node= Node,
             MaxAL=MaxAL,
             ChPar=ChPar,
             KodCB=KodCB,
             atmosphere=atmosphere,
             root = root,
             ## output
             TlSolObs.nrow = TlSolObs.nrow,
             TlSolObs.ncol = TlSolObs.ncol,
             boundary.n = boundary.n,
             balance.ncol = balance.ncol, 
             hQThFlC=hQThFlC,
             TlSolObs=TlSolObs,
             atmOut = atmOut,
             balance = balance,
             boundary = boundary,
             error=error,
             intpck = intpck, 
             dblpck = dblpck,
             PACKAGE="SoPhy", NAOK=TRUE, DUP=FALSE)
 
    ## false only if programming error, cf. PCK.f
    stopifnot(intpck[n.intpck]==999999, dblpck[n.dblpck]==999999.25)
    
    dbl <- DblVec
    names(dbl) <- names.DblVec
    names(IntVec) <- names.IntVec
    IntVec <- as.list(IntVec)

    if (((error==7) || (error==6)) && !is.null(intermediate.result)) {
      intermediate.result(as.list(dbl)$time,
                          array(hQThFlC, dim=c(6, NumNP, 2 + MPL),
                                dimnames=list(c("h","Q","th","vx","vz","Conc"),
                           NULL,NULL))[,,IntVec$brk.tprint, drop=FALSE])
      ## PLevel, see PCK.f
      if (error==6) error <- 13
    }

    if (error!=7) {
      if (error==13 && !is.null(message)) {
        if (!message(as.list(dbl)$time)) {
          error <- 20
          break
        }
      } else break
    }
    
    IntVec$start <- FALSE
    if (error==13) {
      IntVec$ii <- 0
      IntVec$iiprint <- IntVec$iter.print
    }
    IntVec <- as.integer(unlist(IntVec))
    error <- as.integer(0)
  }

  if (error!=0) {
    txt <- switch(error,
                  "NSeep exceeded",
                  "NumSPD exceeded",
                  "NDrD exceeded",
                  "NElDrD exceeded",
                  "jj -- should never happen", #5
                  "breakpoint && maxLoop -- should never get here",
                  "breakpoint -- should never get here",
                  "NObsD exceeded",
                  "NumKD exceeded",
                  "NodInf failed", #10
                  "No steady state solution found",
                  "too many iterations (orthomin)",
                  "iteration > max. iter., see swms2d (water)",
                  "nPLvl > MPL + 2",
                  "nbalnc > balaR", #15
                  "",
                  "",
                  "",
                  "",
                  "stopped", #20
                  )
    return(if (error<20) paste("error", error, "in swms_2d occured:", txt) else
           paste("swms_2d:", txt) 
           )
  }

  hQThFlC <- array(hQThFlC, dim=c(6, NumNP, 2 + MPL),
                          dimnames=list(c("h","Q","th","vx","vz","Conc"),
                            NULL,NULL))[, , -2-MPL, drop=FALSE] # -2-MPL,13.11.04
  # hQThFlC <- hQThFlC[, , !apply(is.na(hQThFlC), 3, all)] # 13.11.04

  nx <-  function(x) paste(x, c("Atm", "Root", 3, 1, "Seep", 5:NumKD), sep="")
  TlSolObs <-
    matrix(TlSolObs, nrow=TlSolObs.nrow, ncol=TlSolObs.ncol,
           dim=list(NULL, c("Time", "rAtm", "rRoot", nx("vK"), nx("hK"),  #TLInf
             "CumQAP", "CumQRP", nx("CumQ"),
             "dt", "Iter", "ItCum", "Peclet", "Courant",
             "CumCh0", "CumCh1", "CumChR", paste("ChemS",1:NumKD, sep=""),#SolInf
             paste("SMean", 1:NumKD, sep=""), if (NObs>0) #obsnod
             c(paste("hNew", 1:NObs, sep=""), paste("ThNew", 1:NObs, sep=""),
               paste("ConcNew", 1:NObs, sep=""))
             ))
           )
  TlSolObs <- TlSolObs[!apply(is.na(TlSolObs), 1, all), ]

  atmOut <-
    matrix(atmOut, nrow=MaxAL, ncol=9,
           dim=list(NULL, c("AtmTime", "CumQAP", "CumQRP","CumQA",
             "CumQR", "CumQ3", "hAtm","hRoot","hKode3")))
 
  name <- c("time","WatBalT", "WatBalR", "CncBalT", "CncBalR")
  name.pref <- c("Area","Volume","InFlow","hMean")
  if (lChem) name.pref <- c(name.pref,"ConcVol","cMean")
  for (i in 1:length(name.pref))
    name <- c(name, paste(name.pref[i], c("tot",1:d$NLay),sep=""))

  balance <-
    matrix(balance, nrow=boundary.n, ncol=balance.ncol,
           dim=list(NULL, name))
  balance <- balance[!apply(is.na(balance), 1, all), ]
  
  boundary <- array(boundary, dim=c(NumBP, 11, boundary.n),
                           dimnames=list(NULL,c("i","n","x","z","Code","Q",
                             "v","h","th","Conc","time"),NULL))
  boundary <- boundary[, ,!apply(is.na(boundary), 3, all)]

  return(c(h, list(hQThFlC=hQThFlC, TlSolObs=TlSolObs, atmOut=atmOut,
              balance=balance, boundary=boundary,
              flux=d$flux, water.x=d$water.x, water.y=d$water.y
              )))
}



read.swms2d.table <- function(path,
                              selct.in = "SELECTOR.IN",
                              grid.in = "GRID.IN",
                              atm.in = "ATMOSPH.IN"
                              ){
  NSeepD <- 2
  NumSPD <- 50
  NDrD <- 2
  NElDrD <- 8
  NTabD <- 100
  NumKD <- 6   ## always >= 5 !!
  NObsD <- 4
  MNorth <- 4
  
  ENVIR = environment()
  READLINES <- 0
  
  read <- function(file, skip, nrows=1, com="(", header=TRUE,
                   fct=function(x)x, solve=header, mode="any",
                   ndata=NA) {    
    assign("READLINES", 0, envir=ENVIR)
    stopifnot(((nrows==1) && !header && !solve) || is.na(ndata))
    x <- read.table(file=file, skip=skip, header=header, nrows=nrows, com=com)
    if (solve)
      for (i in names(x)) {
        if (any(i==ls(envir=ENVIR))) stop(paste(i,"exists"))
        txt <- paste("assign('",i,"',fct(x$",i,"),envir=ENVIR)",sep="")
        eval(parse(text=txt))
      }
    else {
      if (is.na(ndata)) return(x)
      x <- as.numeric(x)
      assign("READLINES", 1, envir=ENVIR)
      while (length(x) < ndata) {
        x <- c(x, as.numeric(read.table(file=file, skip=skip+READLINES,
                                        header=header, nrows=1, com=com)))
        assign("READLINES", READLINES + 1, envir=ENVIR)
      }
      if (length(x) > ndata) warning(paste("too many data in",match.call()))
      return(x[1:ndata])
    }
  }
  
  slct <- paste(path,selct.in,sep="/")
  grid <- paste(path,grid.in,sep="/")
  atmos <- paste(path,atm.in,sep="/")

  read(slct, sk=3, fct=as.character)  # units
  if (!exists("BUnit")) BUnit <- NULL

  read(slct, sk=5)  # kat
  read(slct, sk=7)  # maxit, tolthX
  read(slct, sk=9, fct=function(x)x=="t") #lWat lChem CheckF ... 
  read(grid, sk=1)  # NumNP     NumEl       IJ      NumBP     NObs
  nCodeM <- read(grid, sk=3, nrows=NumNP, solve=FALSE)
  ## n Code x z h Conc Q M B Axz Bxz Dxz
  
  KXR <- read(grid, sk=NumNP+5, nrows=NumEl, solve=FALSE)
  stopifnot(ncol(KXR)==9)
  gskip <- NumNP+NumEl+8
  KXB <- read(grid, sk=gskip, h=FALSE, ndata=NumBP)
  if ((length(kxb <- nCodeM[ nCodeM[,2]!=0, 1]) != length(KXB)) ||
      (any(KXB!=kxb)))
    warning("KXB does not match Kode information")
  gskip <- gskip + 1 + READLINES
  Width <- read(grid, sk=gskip, h=FALSE, ndata=NumBP)
  gskip <- gskip + 1 + READLINES
  rLen <- read(grid, sk=gskip, h=FALSE)
  ## last of grid !!
  Node <- if (NObs>0) read(grid, sk=gskip+2, h=FALSE, ndata=NObs) else NULL
  read(slct, sk=12) # NMat    NLay    hTab1   hTabN   NPar
  hTab <- c(hTab1, hTabN)
  
  Par <- read(slct, sk=14, nrow=NMat, solve=FALSE)
                                        #thr ths tha  thm  Alfa n Ks Kk thk
  if (AtmInF) {
    read(atmos, sk=3) # SinkF   qGWLF
    read(atmos, sk=5) # GWL0L   Aqh     Bqh
    read(atmos, sk=7) # tInit   MaxAL
    read(atmos, sk=9) # hCritS
    atmosph <- read(atmos, sk=11, nrows=MaxAL, solve=FALSE)
    stopifnot(ncol(atmosph)>=8,ncol(atmosph)<=10,
              nrow(atmosph)==MaxAL)    
    atmosphere <- matrix(ncol=10, nrow=MaxAL)
    atmosphere[, 1:ncol(atmosph)] <- as.matrix(atmosph)
  } else {
    SinkF <- qGWLF <- FALSE
    atmosphere <- GWL0L <- Aqh <- Bqh <- tInit <- hCritS <- NA
    MaxAL <- 0
    }
  skip <- 16 + NMat
  read(slct, sk=skip) # dt  dtMin dtMax DMul DMul2 MPL
  dtMinMax <- c(dtMin, dtMax)
  DMul <- c(DMul, DMul2)
  skip <- skip+3
  TPrint <- read(slct, sk=skip, h=FALSE, ndata=MPL)
  skip <- skip + READLINES    
  if (SinkF) {
    root <- read(slct, sk=skip+2, h=FALSE, ndata=6)
    skip <- skip + READLINES + 1 + 2
    POptm <- read(slct, sk=skip, h=FALSE, ndata=NMat)
    skip <- skip + READLINES
  } else {
    root <- POptm <- NA
  }
  if (SeepF) {
    ## only NP will be used further, other information can be reconstructed
    read(slct, sk=skip+1) # NSeep
    NSP <- read(slct, sk=skip+4, h=FALSE, ndata=NSeep)
    NP <- matrix(NA, nrow=NSeep, ncol=NumSPD) 
    select <- NULL
    for (i in 1:NSeep) select <- rbind(select, cbind(i,1:NSP[i]))
    NP[select] <- read(slct, sk=skip+6, h=FALSE, ndata=sum(NSP, na.rm=TRUE))
    skip <- skip + 6 + READLINES
  } else {
    NSP <- NP <- NSeep <- NA
  }
  
  if (DrainF) {
    ## only EfDim & KElDr will be used further, other information can be
    ## reconstructed
    read(slct, sk=skip+1) # NDr, DrCorr
    skip <- skip + 4
    ND <- read(slct, sk=skip, h=FALSE)
    skip <- skip + READLINES + 1
    NED <- read(slct, sk=skip, h=FALSE)
    skip <- skip  + READLINES + 1
    EfDim <- matrix(read(slct, sk=skip, nrows=NDr, h=FALSE), ncol=2)
    skip <- skip + 3
    KElDr <- matrix(nrow=NDr, ncol=NElDrD) 
    select <-  NULL
    for (i in 1:NDr) select <- rbind(select, cbind(i,1:NED[i]))
    KElDr[select] <- read(slct, sk=skip, h=FALSE, ndata=sum(NED, na.rm=TRUE))
    skip <- skip + NDr
  } else {
    NDr <- DrCorr <- ND <- NED <- EfDim <- KElDr <- NA
  }
  
  if (lChem) {
    read(slct, sk=skip+1) # Epsi,lUpW,lArtD,PeCr
    lUpW <- lUpW=="t"
    lArtD <- lArtD=="t"
    ChPar <- read(slct, sk=skip+4, h=FALSE, nrow=NMat)
    skip <- skip + 5 + NMat
    KodCB <- read(slct, sk=skip, h=FALSE, ndata=NumBP)
    skip <- skip + READLINES + 1
    cBound<- as.numeric(read(slct, sk=skip, h=FALSE))
    stopifnot(length(cBound)==6)
    read(slct, sk=skip+1) # tPulse
    skip <- skip + 3      
  } else {
    Epsi <- lUpW <- lArtD <- PeCr <- ChPar <- KodCB <- tPulse <- NA
    cBound <- rep(NA, 6)
  }
 
  return(list(Units=c(LUnits=LUnit, TUnit=TUnit, MUnit=MUnit, BUnit=BUnit),
              Kat=Kat,
              MaxIt=MaxIt,
              TolTh=TolTh,
              TolH=TolH,
              lWat=lWat,     
              FreeD=FreeD,

              ## material
              NLay=NLay, 
              hTab=hTab, # c(hTab1, hTabN)
              Par=Par, # material information
              
              dt=dt,
              dtMinMax=dtMinMax,  # dtmin, dtmax
              DMul=DMul, # 
              TPrint=TPrint,

              ## root(sink information): see after atmoshere
              
              ## seepage
              NP=NP,

              ## drains
              DrCorr=DrCorr,
              ND=ND,
              EfDim=EfDim,
              KElDr=KElDr,

              ## chemical
              Epsi=Epsi,
              lUpW=lUpW,
              lArtD=lArtD,
              PeCr=PeCr,
              ChPar=ChPar,
              KodCB=KodCB,
              cBound=cBound,
              tPulse=tPulse,
              
              ## nodes
              nCodeM=nCodeM,

              ## elements
              KXR=KXR,

              ## boundary
              Width=Width,
              rLen=rLen,              
              Node=Node,

              ## atmosphere
	      SinkF=SinkF,
              qGWLF=qGWLF,
              GWL0L=GWL0L,
              Aqh=Aqh,
              Bqh=Bqh,
              tInit=tInit,
              hCritS=hCritS,
              atmosphere=atmosphere,

              ## Block D table 8.4, logically after atmosphere
              ## (sink information)
              root=root, # 6 variables
              POptm=POptm,
              
              ## additional information
              path = path,
              selct.in = selct.in,
              grid.in = grid.in,
              atm.in = atm.in
              ))
}








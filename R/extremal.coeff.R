

rpareto <- function(n, xi, s=1) {
  if (length(n)>1) n <- length(n)
  xi <- rep(xi, len=n)
  s <- rep(s, len=n)
  res <- rep(NA, n)
  allowed <- s>0
  if ((sidx <- sum( idx <- xi==0 & allowed)) > 0) {
    res[idx] <- rexp(sidx, 1/s)
  }
  if ((sidx <- sum( idx <- xi>0 & allowed)) > 0) {
    res[idx] <-  ((1.0 - runif(sidx))^(-xi) - 1) * s / xi
  }
  if ((sidx <- sum( idx <- xi<0 & allowed)) > 0) {
    res[idx] <- (1.0 - (1.0 - runif(sidx))^(-xi)) * s / xi
  }
  return(res)
}

ppareto <- function(q, xi, s=1, lower.tail=TRUE) {
  l <- max(length(q), length(xi), length(s))
  q <- rep(q, len=l)
  xi <- rep(xi, len=l)
  s <- rep(s, len=l)
  res <- rep(1, l)
  res[s<=0] <- NA
  res[q / s * xi <= -1] <- 0
  allowed <- q>0 & s>0 & (xi>=0 | q / s * xi > -1)
  idx <- xi==0 & allowed
  res[idx] <- exp(-q[idx] / s[idx])
  idx <- xi!=0 & allowed
  xi <- xi[idx]
  res[idx] <- (1 + xi * q[idx] / s[idx])^(-1/xi)
  if (lower.tail) res <- 1 - res
  return(res)
}


dpareto <- function(x, xi, s=1) {
  l <- max(length(x), length(xi), length(s))
  x <- rep(x, len=l)
  xi <- rep(xi, len=l)
  s <- rep(s, len=l)
  res <- rep(0, l) 
  res[s<=0] <- NA
  allowed <- x>0 & s>0 & (xi>=0 | x / s * xi > -1)
  idx <- xi==0 & allowed
  res[idx] <- exp(-x[idx] / s[idx]) / s[idx]
  idx <- xi!=0 & allowed
  xi <- xi[idx]
  s <- s[idx]
  res[idx] <- (1 + xi * x[idx] / s[idx])^(-1/xi-1) / s[idx]
  return(res)
}

qpareto <- function(p, xi, s=1, lower.tail=TRUE) {
  if (lower.tail) p <- 1-p
  l <- max(length(x), length(xi), length(s))
  p <- rep(p, len=l)
  xi <- rep(xi, len=l)
  s <- rep(s, len=l)
  res <- rep(NaN, l)
  allowed <- p>=0 & p<=1 & s>0
  idx <- xi==0 & allowed
  res[idx] <- -s * log(p)
  idx <- xi!=0 & allowed
  xi <- xi[idx]
  s <- s[idx]
  res[idx] <-  s / xi * (p^(-xi) - 1)
  return(res)
}


## ausfuehrliche Dokumentation wie data neu sortiert wird und was
## zurueckgegeben wird!!
risk.index <- 
       function(data, # matrix(col1=distance to soil surface, col2=frequency)
                selected.dist=0.95,#less than 1 then interpreted as
                ##              first ... last part
		##               -- distances for which xi is calculated
                selected.rate=cbind(c(0.5, 0.8), c(0.4, 0.9), c(0.3, 1.0)),
                weights=1, ## weights if NA have been in the data, and
		##           data have been corrected accordingly
                ##           In a sequence of weights values < max(weight)
		##           should be the exeption ! -- otherwise the 
                ##           algorithm will not work well
                measure=function(x) x^2, ## distance function for least "squares"
                method=c('fix.m', 'optim.m', 'ml'),
                min.neg.xi = -10, # used for fix.m and optim.m
                max.neg.xi = -0.1,# used for fix.m and optim.m
                max.pos.xi = 10,  # used for fix.m and optim.m
                endpoint.tolerance = 0, # if pos then tol w.r.t. to freq
                ##                       if neg then tol w.r.t. to scale
		front.factor=2, ## 1.2 is too small; results for 2 and
                ##                 20 are not distinguishable
                ##                 (preliminary examination) 
                max.no.paths=10 * max(data[, 2]),
                PrintLevel=RFparameters()$Print,
                max.rate = TRUE
                )
{

#  str(list(data=data,selected.dist=selected.dist, selected.rate=selected.rate,
#           weights=weights, measure= measure,method=method,
#           min.neg.xi=min.neg.xi, max.neg.xi=max.neg.xi, max.pos.xi=max.pos.xi,
#           endpoint.tolerance=endpoint.tolerance, front.factor=front.factor,
#           min.no.paths, max.no.paths, PrintLevel))
 
  ord <- order(data[, 1], decreasing=FALSE)     
  data <- data[ord,  ] ###!!!!!!!!
  weights <- rep(weights, len=nrow(data))
  weights <- weights[ord]

  stopifnot(is.character(method),
            is.numeric(endpoint.tolerance), length(endpoint.tolerance)==1,
            is.numeric(min.neg.xi), length(min.neg.xi)==1,
            is.numeric(max.neg.xi), length(max.neg.xi)==1,
            is.numeric(max.pos.xi), length(max.pos.xi)==1,
            is.numeric(front.factor), length(front.factor)==1,
            is.function(measure))
  method <- match.arg(method)

  stopifnot(ncol(data)==2)
  max.dist <- max(data[, 1])
  mxfreq <- max(data[, 2], na.rm=TRUE)

  ## transform the required selected.dist and the selected.rate
  ## into explicit indices
  if (!is.null(selected.dist)) {
    if (any(selected.dist<1)) {
      stopifnot(all(selected.dist<=1), 
                length(selected.dist)==1 || length(selected.dist) %% 2 == 0)
      if (length(selected.dist)==1)
        selected.dist <- 1 + 0:(round(nrow(data)-1) * selected.dist)
      else 
        selected.dist <- 
          apply(matrix(selected.dist, nrow=2), 2,
                function(i) 1 + round((nrow(data)-1) * i[1]) :
                round((nrow(data)-1) * i[2]))
      selected.dist <- sort(unique(as.vector(selected.dist)))
    }
  } else stopifnot(!is.null(selected.rate))

  
  sel.dist <- selected.dist
  if (!is.null(selected.rate)) {
    selected.rate <- as.matrix(selected.rate)
    stopifnot(nrow(selected.rate)==2, all(selected.rate <= 1),
              all(selected.rate >= 0))
    for (i in 1:ncol(selected.rate)) {
      sel.rate <- range(selected.rate[, i]) * mxfreq
      dummy <- data[, 2]
      if (max.rate) {
        dummy[is.na(dummy)] <- -1
        dummy <- rev(cummax(rev(dummy)))
      }
      idx <- !is.na(data[, 2]) & (dummy >= sel.rate[1]) & (dummy <= sel.rate[2])
      if (any(idx)) break;
    }
    if (!any(idx))
      stop("no value found for the given intervals of selected.rate")
    sel.rate <- which(idx)
  } else {
    sel.rate <- NULL
    raw.risk.index <- risk.index <- NA
  }
  
  sel.dist <- c(sel.dist, sel.rate)
  if (is.null(sel.dist)) {
    total.par <- NA
    values <- NA
  } else {
    sel.dist <- sort(unique(sel.dist))

    ## set target functions and initialisation functions according to method
    if (method=="optim.m") {
      targetLS <- function(eta) {
        ## 1E-50 avoid error messages in case of positive eta and non-working
        ## optim (parameters out of bounds)
        if (eta[1] > max.pos.xi) {
          theo <-
            eta[3] * (pmax(0, 1 + max.pos.xi * dist / eta[2]))^(-1/max.pos.xi)
          res <- sum(w * measure(theo - freq)) + 100 * (eta[1] - max.pos.xi)^2
        } else {
          theo <- eta[3] * (pmax(0, 1 + eta[1] * dist / eta[2]))^(-1/eta[1])
          res <- sum(w * measure(theo - freq)) 
          if (is.finite(res) && res < VALUE) {
            assign("VALUE", res, envir=ENVIR)
            assign("PARAM", eta, envir=ENVIR)
          }
        }
        res
      }
      
      targetLSneg <- function(eta) {
        ## 1E-50 avoid error messages in case of positive eta and non-working
        ## optim (parameters out of bounds)
        theo <- eta[3] * (pmax(0, 1 - dist / eta[2]))^(-1/eta[1])
        res <- sum(w * measure(theo - freq))
        if (is.finite(res) && res < VALUE) {
          assign("VALUE", res, envir=ENVIR)
          assign("PARAM", eta, envir=ENVIR)
        }
        res
      }
      
      target <- list(targetLSexp=function(eta) { sum(w * measure(eta[2]
                       * exp(- dist / eta[1]) - freq))}, # exp
                     targetLS, # pos
                     targetLSneg) # neg
      
      front <- data[max(1, which(data[,2] > 0)), 1]
      nk <- c(1, 10, 1)
      ini <- function(k) {
        ## different scales, for initial xi values of 0, 1, -1
        scale <- switch(k[1],
                        max.dist / nk[k[1]] * k[2],
                        max.dist / nk[k[1]] * k[2],
                        max(1, dist[which(freq>0)]))
        c(switch(k[1], scale, c(1, scale), c(-1, scale)), 1.1*max(0, freq))
      }
      
      ## bounds, depending on xi=0, pos, neg
      lower <- 0.0001
      low <- function(k) {
        c(switch(k, lower, c(0.0001, lower),
                 c(min.neg.xi, 
                   if (endpoint.tolerance <= 0) {
                     max(1, dist[which(freq>0)]) + endpoint.tolerance
                   } else max(1, dist[which(freq > endpoint.tolerance)])
                   )),
          max(0, freq)
          )
      }   
      upper <- Inf
      up <- function(k) {
        c(switch(k, upper, c(max.pos.xi, upper),
                 c(max.neg.xi, front.factor * front)), 
          5 * mxfreq)
      }
      
      ## after optimisation, transform the optimised parameters to
      ## readable, standardised parameters
      param <- function(k, par) {
        if (k==1) par <- c(0, par) else  ## par[2] <- par[2] * par[1]
        if (k==3) par[2] <- - par[2] * par[1]
        return(par)
      }     
    } else if (method=="fix.m") {
      targetLS <- function(eta) {
        ## 1E-50 avoid error messages in case of positive eta and non-working
        ## optim (parameters out of bounds)

        ## theo <-min.no.paths *(pmax(0, 1 + eta[1] * dist / eta[2]))^(-1/eta[1])
        ## max.freq <- max(freq[w==max(w, na.rm=TRUE)], freq[1])
        max.freq <- max(0, freq)
        if (eta[1] > max.pos.xi) {
          theo <- 
            max.freq * (pmax(0, 1 + max.pos.xi * dist / eta[2]))^(-1/max.pos.xi)
          res <- sum(w * measure(theo - freq)) + 100 * (eta[1] - max.pos.xi)^2
        } else {
          theo <- max.freq *(pmax(0, 1 + eta[1] * dist / eta[2]))^(-1/eta[1])
          res <- sum(w * measure(theo - freq)) 
          if (is.finite(res) && res < VALUE) {
            assign("VALUE", res, envir=ENVIR)
            assign("PARAM", eta, envir=ENVIR)
            assign("MAX.FREQ", max.freq, envir=ENVIR)
          }
        }
        res
      }
      targetLSneg <- function(eta) {
        ## 1E-50 avoid error messages in case of positive eta and non-working
        ## optim (parameters out of bounds)
        ## theo <- freq[1] * (pmax(0, 1 - dist / eta[2]))^(-1/eta[1])
        theo <- max(0, freq) * (pmax(0, 1 - dist / eta[2]))^(-1/eta[1])
        res <- sum(w * measure(theo - freq))
        if (is.finite(res) && res < VALUE) {
          assign("VALUE", res, envir=ENVIR)
          assign("PARAM", eta, envir=ENVIR)
        }
        res
      }
      
      target <- list(targetLSexp=function(eta) { sum(w * measure(max(0, freq) *
                       exp(- dist / eta[1]) - freq))}, # exp
                     targetLS, # pos
                     targetLSneg) # neg
      
      front <- data[max(1, which(data[,2]>0)), 1]
      nk <- c(1, 10, 1)
      ini <- function(k) {
        ## different scales, for initial xi values of 0, 1, -1
        scale <- switch(k[1],
                        max.dist / nk[k[1]] * k[2],
                        {
                          basis <- 1.6
                          kneg <- 0.6 * nk[k[1]] 
                          max.dist * basis^(k[2] - kneg)
                          ##  max.dist / nk[k[1]] * k[2],
                        },
                        max(1, dist[which(freq>0)]))
                                        # scale <- diff(range(dist)) / nk * k.m
        switch(k[1], scale, c(1, scale), c(-1, scale))
      }
      
      ## bounds, depending on xi=0, pos, neg
      lower <- 0.0001
      low <- function(k) {       
        switch(k, lower, c(0.0001, lower),
               c(min.neg.xi, 
                 if (endpoint.tolerance <= 0) {
                     max(1, dist[which(freq>0)]) + endpoint.tolerance
                 } else max(1, dist[which(freq > endpoint.tolerance)])
                 ))
      }
      
      upper <- Inf
      up <- function(k) 
        switch(k, upper, c(max.pos.xi, upper),
               c(max.neg.xi, front.factor * front))
      
      ## after optimisation, transform the optimised parameters to
      ## readable, standardised parameters
      param <- function(k, par) {
        if (k==1) par <- c(0, par) else # par[2] <- par[2] * par[1]
        if (k==3) par[2] <- - par[2] * par[1]
        return(par)
      }
      
    } else if (method=="ml") {
      target <- function(x) {
        zero <- 1e-30
        xi <- x[1]
        s <- x[2]
        n <- x[3]
        - sum(lgamma(n + 1) - lgamma(n - freq + 1)
              - freq / xi * log(pmax(zero, 1 + xi * ML.DIST / s))
              + (n - freq) * log(1 - pmax(zero, 1 + xi * ML.DIST / s)^(-1 / xi)))
      }
      target0 <- function(x) {
        zero <- 1e-30
        s <- x[1]
        n <- x[2]
        - sum(lgamma(n + 1) - lgamma(n - freq + 1) - freq * ML.DIST / s +
              (n - freq) * log(1 - exp(-ML.DIST / s)))
      }
      
      target <- list(target0, # exp
                     target, # pos
                     target) # neg
    
      front <- data[max(1, which(data[,2]>0)), 1]
      nk <- c(1, 1, 1)
      ini <- function(k)
         c(switch(k[1], NULL, 1, -1),
           switch(k[1], max.dist/5, max.dist/5, max(0, ML.DIST[which(freq>0)])),
           max(0, freq) + 1
           )
      
      ## bounds, depending on xi=0, pos, neg
      lower <- 0.0001
      low <- function(k) c(switch(k, NULL, lower, -Inf), lower,
                           max(0, freq)+0.001)
      up <- function(k) c(switch(k, NULL, Inf, -lower), Inf, Inf)
    
      ## after optimisation, transform the optimised parameters to
      ## readable, standardised parameters
      param <- function(k, par) return(if (k==1) c(0, par) else par)
      
    } else stop(paste("unknown method:", method))
    
    total.par <- matrix(nrow=3, ncol=length(sel.dist)) ## total.par
    ## stored the optimal parameter for each distance, above which data
    ## are used for fitting

    lsq <- list()  ## collection of the initial optimisation results
    factor <- 100 / length(sel.dist)
    values <- numeric(length(sel.dist)) ## optim()$value
    
    mess <- c("exp", "pos", "neg")
    nkcum <- cumsum(c(0, nk))
    ENVIR <- environment()
    MAX.FREQ <- NULL
    for (i in 1:length(sel.dist)) {
      ## select the data, whose distance to the soil surface is >= i
      idx <-  sel.dist[i] : nrow(data)
      idx <- idx[is.finite(data[idx, 2])]
      dist <- data[idx, 1] - min(data[idx, 1])
      ML.DIST <- dist + 1
      freq <- data[idx, 2]
      if (sum(freq>0) <= 1) {
        for (signum in 1:3) { # 0, +, -
          segm <- nkcum[signum]
          for (k in 1:nk[signum]) {
            lsq[[segm + k]] <- list(value=NA) 
          }
        }
      } else {    
        w <- weights[idx]
        if (PrintLevel>1) {
          cat("\b\b\b\b", format(round(i * factor), dig=2),"%", sep="")
        }
        
        for (signum in 1:3) { # 0, +, -
          segm <- nkcum[signum]
          for (k in 1:nk[signum]) {
            if (PrintLevel>6) cat(mess[[signum]], " ", k, "\n")
            ## optimisation for starting value of xi=0,1,-1 and different scales
            ## (( the choice of the scale is critical -- optim may run into
            ##    local minima -- observed previously; xi does not seem to be
            ##    that crititcal))
            PARAM <- ini(c(signum, k))
            VALUE <- target[[signum]](PARAM)
            init <- ini(c(signum, k))
            zaehler <- 0
            # print(init)
#            print(PARAM)
#            print(signum)
#            print(target[[signum]])
#            print(low(signum))
#            print(up(signum))
            
            while
            (!is.list(lsq[[segm + k]] <-
                      #try
                      (optim(PARAM, fn=target[[signum]], lower=low(signum), 
                                upper=up(signum), meth="L-BFGS-B",
                                control=list(fnscale=c(if (signum!=1) 1,
                                               max.dist/3, if (method!="fix.m")
                                               max(0, data[, 2])))
                                ))))
            {
              ## it may happen that the optim algorithm fails without any real
              ## reason, so the optimisation is restarted using the currently
              ## best parameters as the initial ones.
              ## retry it only 4 times before giving up. 
              if (PrintLevel>3)
                cat("repeat: i=", i,
                    " threshold d=", sel.dist[i],
                    " sign=", signum,
                    " k=", k,
                    "\n", sep="")
              if ((zaehler <- zaehler + 1) > 4) break;
            }
 
            lsq[[segm + k]]$fixed.m <- MAX.FREQ
            if (!is.list(lsq[[segm + k]])) {
              lsq[[segm + k]] <- list()
              if (PrintLevel>3) 
                cat(" optim failed partially for sel.distance=", i, 
                    "; xi:", mess[signum], "; ini=", ini(c(signum, k)), "\n")
            }
            if (VALUE<Inf &&
                (length(lsq[[segm + k]])==0 || lsq[[segm + k]]$value>VALUE)) {
              if (PrintLevel>3) cat("  undetected optimum\n")
              lsq[[segm + k]]$value <- VALUE
              lsq[[segm + k]]$par <- PARAM
             }
          }  # for k
        } # for signum
      } # else any freq>0

      value <- sapply(lsq, function(x) if (is.null(x$value)) NA else x$value)

      ## extract best parameter out of the three cases -,0,+
      if (all(is.na(value))) {
        eta <- total.par[, i] <- NA     
      } else {
        idx <- which(value==(values[i] <- min(value, na.rm=TRUE)))[1]
        idx.i <- sum(idx > nkcum)
         total.par[, i] <-
          c(param(idx.i, lsq[[idx]]$par), lsq[[segm + k]]$fixed.m)        
        if (PrintLevel>2) {
          
          #for (jj in 2:11) cat(param(2, lsq[[jj]]$par)[1], "")
          #cat("\n")
          
          eta <- total.par[, i]
          print(log(value))
          print(log(value/values[i]))
          print(c(idx, NA, nkcum, NA, idx.i))
          print(eta)
          if (PrintLevel>4) {
            if (!is.logical(close.screen())) screen(max(close.screen()))
            plot(dist, freq)
            if (FALSE)
            print(cbind(dist, switch(method, fix.m=max(0, freq), optim.m=eta[3])
                        *(pmax(0, 1 + eta[1] * dist /
                               eta[2]))^(-1/eta[1]))
                  )
            lines(dist, switch(method, fix.m=max(0, freq), optim.m=eta[3],
                               ml=dist)
                  *(pmax(0, 1 + eta[1] * dist /
                         eta[2]))^(-1/eta[1]), col="red")
            lines(dist, switch(method, fix.m=freq[1], optim.m=eta[3], ml=dist) * 
                  ppareto(dist, xi=eta[1], s=eta[2], l=FALSE),
                  lty=2, col="blue")
            if (PrintLevel>5) readline()
          }
        }
      }
    }

    ## extract risk.index 
    if (is.null(selected.rate)) risk.index <- NULL
    else {
      idx <- (sel.dist %in% sel.rate) & is.finite(total.par[1,])
      raw.risk.index <-
        if (any(idx)) median(total.par[1,idx], na.rm=TRUE) else NA
      idx <- (idx &
              (total.par[1,] < max.pos.xi * 0.999) &
              (total.par[1,] > min.neg.xi * 0.999))
      risk.index <-
        if (any(idx)) median(total.par[1,idx] , na.rm=TRUE) else NA  
    }
    total.par <- rbind(total.par, sel.dist)
    dimnames(total.par) <- list(c("xi", "s", "m", "dist (pixels)"), NULL)
    
  } # not is.null(sel.dist)
  
  return(list(par=total.par,
         #lsq, ## lsq is several times overwritten -- only last run for largest
         ##   distance is returned; lsq is only for debugging and developing
         ##   checks
              data=data, weights=weights, 
              selected.dist=selected.dist,
              selected.rate=selected.rate,
              sel.rate=sel.rate,
              sel.dist=sel.dist, ## the distances used (Klartext und als
              ##              Startindex fuer frequency Vektor)
              max.freq = mxfreq,
              values=values,
              method = method,
              measure = measure,
              raw.risk.index = raw.risk.index,              
              risk.index = risk.index
              ))
}


analyse.profile <-
  function(picture, fct, param, lower, upper, loc,          
           estimate.all=NULL, selected.dist=0.95, selected.rate=c(0.5, 0.8),
           measure=function(x) x^2, method=c("fix.m","optim.m", "ml"),
	   endpoint.tolerance=-10,
           ppi.par = 5, ## point per interval
           ##                                    names important, not sequence
           ppi.xy = c(xlow=2, xup=2, ylow=2, yup=2), ##' sequ important not name
          
           interactive=TRUE, PrintLevel=RFparameters()$Print,
           extensions=c("tif", "tiff", "gif", "jpeg"),
           
           X11.width=13, new=TRUE,
           col.thresh=c("white", "yellow", "black", "blue"),
           col.rect="red", col.bg="yellow", col.sep="gray",
           col.left="red", col.mid="darkgreen", col.right="darkgreen",
           col.line = "red", col.txt="black",
           reverse=TRUE, title
          )
{
  debug <- FALSE # debug <- TRUE
  
  line.tolerance <- 2 ## this parameter might be put into the parameter
  ##                       list of the function
  bDa <- 0.15 ## relation of the angles to switch between corner and line
  
  
  markemptylines <- function(picture) { 
    ll <- apply(picture[,,1]==255 & picture[,,2]==255 &
                picture[,,3]==255, 2, sum)
      linie <- ll > dim(picture)[2] * 0.7
    picture[,linie,1] <- NA
    picture[,linie,2] <- NA
    picture[,linie,3] <- NA
    return(picture)
  }
 
### reading in picture  
  if (debug) cat("reading in picture\n")
#  if (.Platform$OS.type!="unix")
#    stop("'analyse.profile' only available on unix systems")
 
  if (is.list(picture)) {
    stopifnot(missing(fct) || is.null(fct), missing(param) || is.null(param),
              missing(upper) || is.null(upper), missing(lower) || is.null(lower),
              missing(loc) || is.null(loc))
    stopifnot(!is.null(picture$fct),
              !is.null(picture$param), !is.null(picture$lower),
              !is.null(picture$upper), !is.null(picture$loc))
    if (is.null(picture$picture)) {
      stopifnot(!is.null(picture$name))
      picture$picture <- read.picture(picture$name, ex=extensions, Pr=PrintLevel)
      picture$picture <- markemptylines(picture$picture)
    }
    dp <- dim(picture$picture)
    fct <- picture$fct
  } else {
### define threshold function for extracting blue area
    if (debug) cat("define threshold function for extracting blue area\n")
    if (missing(fct) || !is.function(fct)) {
      stopifnot(missing(fct) || is.null(fct))
      fct <-
        function(i, p) {
          if (missing(p)) list(minR=160, maxGB=200)
          else {
            gb <- i[,,2] + i[,,3]
            (gb >= i[,,1] * p$minR / 100) & (80<=gb) & (p$maxGB>=gb)
          }
        }
      environment(fct) <- NULL
    }
    
### check input of param, lower and upper
    if (debug) cat("check input of param, lower and upper\n")
    if (PrintLevel>7) cat("param\n")
    if (missing(param)) param <- fct()
    stopifnot(is.list(param))
    if (is.null(names(param)))
      names(param) <- paste("parameter", 1:length(param), sep="")
    if (missing(lower)) lower <- param
    else {
      stopifnot(is.list(lower))
      names(lower) <- names(param)
    }
    if (missing(upper)) upper <- param
    else {
      stopifnot(is.list(upper))
      names(upper) <- names(param)
    }
    
### for compatibility to older versions only -- delete later on
    if (debug)
      cat("for compatibility to older versions only -- delete later on\n")
    if (!is.null(param$minGB)) {
      param$minGB <- NULL
      upper$minGB <- NULL
      lower$minGB <- NULL
    }
    
### read in picture    
    if (debug) cat("read in picture\n")
    if (PrintLevel>7) cat("picture\n")
    if (is.character(picture)) {
      name <- paste(path.expand(dirname(picture)), basename(picture), sep="/")
      picture <- read.picture(name, extensions=extensions, PrintLevel=PrintLevel)
    } else {
      stopifnot(is.array(picture))
      name <- "unnamed array"
    }
    dp <- dim(picture)
    if (length(dp)==0) stop("invalid picture format")
    else stopifnot(length(dp)==3, dp[3] %in% c(3,4))
    picture <- markemptylines(picture)
        
### define list of data
    if (debug) cat("define list of data\n")
    picture <- list(picture=picture,
                    name=name,
                    fct=fct,               
                    param = if (missing(param)) fct() else param,
                    lower = if (missing(lower)) fct() else lower,
                    upper = if (missing(lower)) fct() else upper,
                    loc = if (!missing(loc) && !is.null(loc)) {
                      stopifnot(is.list(loc))
                      if (is.list(loc[[1]])) loc else list(loc)
                    } else NULL
                    )
  } # else picture not list

  if (!is.null(ppi.par))
     picture$ppi.par <-
       if (is.list(ppi.par)) ppi.par
       else lapply(picture$param, function(x) ppi.par[1])
  if (!is.null(ppi.xy)) picture$ppi.xy <- ppi.xy


### extracting blue area
  if (debug) cat("extracting blue area\n")
  ENVIR <- environment()
  assign("old.picture", picture$fct(picture$picture, picture$par), envir=ENVIR)

### interactive plot

  if (debug) cat("interactive plot\n")
  if (debug) cat("\n")
  if (interactive) {
    if (PrintLevel>7) cat("interactive\n")
    assign("old.param", picture[c("param", "lower", "upper")], envir=ENVIR)
    X11.height <-
      X11.width / nrow(picture$picture) * ncol(picture$picture) * 40 / 99
    if (new) get(getOption("device"))(height=X11.height, width=X11.width)
    else bg.save <- par()$bg

    figs <- rbind(c(0.01,0.40,0.01,0.99),
                  c(0.41,0.80,0.01,0.99),
                  c(0.81,0.99,0.01,0.99)
                  )
    scr <- split.screen(figs=figs) # , erase=FALSE) ##29.12.03 : erase=-FALSE !!
    orig.dev <- scr[1]
    thresh.dev <- scr[2]
    menue.dev <- scr[3]

    old.options <- options()["locatorBell"]
    options(locatorBell=FALSE)
    on.exit({
      close.screen(rev(close.screen())[1:3]);
      # if (new) dev.off() else par(bg=bg.save);
      options(old.options)
    })

    quadratic <- function(d, v, a, mini=0, maxi=Inf) {
      d <- pmin(1, pmax(0, d)) - 0.5
      d <- ((d>0) * 2 - 1) * d^2 * a * 4
      if (missing(v)) d else pmax(mini, pmin(maxi, v + d))
    }
   
    linear <- function(d, v, mini=1, slope=8) {
      d <- (pmin(1, pmax(0, d)) - 0.5) * slope
      if (missing(v)) d else pmax(mini, v + d)
    }

    y.coord <- if (reverse) dp[2]:1 else 1:dp[2]
    if (y1Gy2 <- y.coord[1]>y.coord[2]) y.coord <- -y.coord
    Axes <- function(){
      axis(1, xaxs="i")
      pry <- pretty(y.coord)
      axis(2, at=pry, labels=-pry, yaxs="i")
    }
    thresh <- function(fig, what=NULL) {
      ## check first which parameter has been changed?
      min.parameter <- any(unlist(fig$lower)!=unlist(old.param$lower))
      max.parameter <- any(unlist(fig$upper)!=unlist(old.param$upper))
      result <- fct(fig$picture, if (min.parameter) fig$lower else
                  if (max.parameter) fig$upper else fig$par)
      if (what=="all") {
        close.screen(c(orig.dev, thresh.dev, menue.dev))
        scr <- split.screen(figs=figs)
        screen(orig.dev)
        par(mar=c(2, 2, 0.2, 0.2))
        plotRGB(picture$picture, reverse=reverse)
        screen(thresh.dev)
        par(mar=c(2, 2, 0.2, 0.2))
        screen(menue.dev)
        par(mar=c(0.1, 0.1, 0.1, 0.1))  
      }
      screen(thresh.dev)
      image(1:dp[1], y.coord, 2 * result + old.picture, axes=!y1Gy2,
            col=col.thresh, zlim=c(-0.5, 3.5))
      if (y1Gy2) Axes()
      assign("old.param", fig[c("param", "lower", "upper")], envir=ENVIR)
      if (!min.parameter && !max.parameter)
        assign("old.picture", result, envir=ENVIR)
      return(fig)
    }
    screen(orig.dev)
    par(mar=c(2, 2, 0.2, 0.2))
    plotRGB(picture$picture, reverse=reverse)
    
    entry <-
      c(list(list(name=
                  if (missing(title)) rev(strsplit(picture$name, "/")[[1]])[1]
                  else title,
                  var=NULL, val=NULL, col=col.sep)),
        lapply(names(picture$param),
               function(x) list(names=x, var=paste("param$", x, sep=""),
                                delta=TRUE, 
                                val=function(d,v) quadratic(d=d, v=v, a=50))),
        list(list(name="refresh", var=NULL, val="simulate",  col=col.rect,
                  what="all")),
        lapply(names(picture$param),
               function(x) list(names=paste("lower bound",x),
                                var=paste("lower$", x, sep=""),
                                delta=TRUE, 
                                val=function(d,v) quadratic(d=d, v=v, a=50))),
        list(list(name="-----", var=NULL, val=NULL)),       
        lapply(names(picture$param),
               function(x) list(names=paste("upper bound", x),
                                var=paste("upper$", x, sep=""),
                                delta=TRUE, 
                                val=function(d,v) quadratic(d=d, v=v, a=50))),
        )
    
    screen(thresh.dev)
    par(mar=c(2, 2, 0.2, 0.2))
    screen(thresh.dev)
    image(1:dp[1], y.coord, old.picture, col = col.thresh[c(1, 4)],
          axes = !y1Gy2)
    if (y1Gy2) Axes()
    screen(menue.dev)
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    par(bg = "white")

    picture <-
      eval.parameters("fig", entry, update = TRUE, 
                      dev = menue.dev, create = FALSE, simulate = thresh, 
                      fig = picture, col.rect = col.rect, col.bg = col.bg, 
                      col.sep = col.sep, col.left = col.left, col.mid = col.mid, 
                      col.right = col.right, col.line = col.line,
                      col.txt = col.txt, what = "partial")
    ## choosing blue area
    ## `loc' for abbreviation only
    loc <- picture$loc
    if (is.null(loc)) {
      screen(menue.dev)
      plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
      text(0.5, 0.5, lab="choose rectangle by 2 points", col=col.txt)
    }     
    screen(thresh.dev, new=FALSE)
    par(new=TRUE)
    plot(Inf, Inf,  xlab="", ylab="", xaxs="i", yaxs="i", axes=FALSE,
         xlim=c(1, dp[1]), ylim=c(1,dp[2]))
    colrs <- c("red", "orange", "brown")
    if (is.null(loc)) loc[[1]] <- Locator(n=2, type="p", col=colrs[1])
    if (length(loc[[1]]$x)==2) {
      loc[[1]]$x <- range(loc[[1]]$x)
      loc[[1]]$y <- range(loc[[1]]$y)
      for (i in 1:3) {
        repeat {
          #str(loc)
         # points(300,200, pch=16, cex=2, col="green")
          #points(300,-200, pch=16, cex=2, col="magenta")
          
          screen(thresh.dev, new=FALSE)
          par(new=TRUE)
          plot(Inf, Inf,  xlab="", ylab="", xaxs="i", yaxs="i", axes=FALSE,
               xlim=c(1, dp[1]), ylim=c(1, dp[2]))
          for (j in 1:i)
            rect(loc[[j]]$x[1], loc[[j]]$y[1], loc[[j]]$x[2], loc[[j]]$y[2],
                 border=colrs[j])
          
          screen(menue.dev)
          par(mar=c(2, 1.2, 0.2, 0))
          if (!is.null(try(
            plot(apply(old.picture[loc[[i]]$x[1]:loc[[i]]$x[2],
                             loc[[i]]$y[1]:loc[[i]]$y[2]],
                     2, sum, na.rm=TRUE), loc[[i]]$y[1]:loc[[i]]$y[2],
               pch=16, cex=0.7, xaxs="i", yaxs="i", ylim=c(1,dp[2]))
                           ))) { ## fehler!!
            print("error !! -- mail the following output to schlather@cu.lu")
            plot(0,0)
            cat("analyse.profile : \n")
            str(list(old.picture=old.picture, loc=loc, dp))
          }
          yyy <- round(loc[[i]]$y)
          if (reverse) yyy <- 1 + dp[2] - rev(yyy)
          text(10, 5, adj=c(0, 0),
               paste(c("optimal", "lower", "upper")[i],
                     ": x=c(",round(loc[[i]]$x[1]), ",", round(loc[[i]]$x[2]),
                     "); y=c(", yyy[1], ",", yyy[2],
                     ")", sep=""), col=colrs[i])
         
          screen(thresh.dev, new=FALSE)
          par(new=TRUE)
          plot(Inf, Inf,  xlab="", ylab="", xaxs="i", yaxs="i", axes=FALSE,
               xlim=c(1, dp[1]), ylim=c(1, dp[2]))
          if (length(x <- Locator(n=1, type="p", col="green"))<1) break;

#          points(x$x, x$y, col="yellow", pch=20, cex=2)
          
          # algorithm to guess which corner should be changed into which position
          x$x <- min(dp[1], max(1, x$x))
          x$y <- min(dp[2], max(1, x$y))

          inside <- FALSE 
          if(x$x<loc[[i]]$x[1]) loc[[i]]$x[1] <- x$x else
          if(x$x>loc[[i]]$x[2]) loc[[i]]$x[2] <- x$x else
          inside <- TRUE
          
          if(x$y<loc[[i]]$y[1]) loc[[i]]$y[1] <- x$y else
          if(x$y>loc[[i]]$y[2]) loc[[i]]$y[2] <- x$y else 
          if (inside) {
            mn <- c(mean(loc[[i]]$x), mean(loc[[i]]$y))
            ru <- c(loc[[i]]$x[2], loc[[i]]$y[1]) - mn
            ru <- ru / sqrt(sum(ru^2))
            ro <- c(loc[[i]]$x[2], loc[[i]]$y[2]) - mn
            ro <- ro / sqrt(sum(ro^2))
            pnt <- c(x$x, x$y) - mn
            pnt <- pnt / sqrt(sum(pnt^2))
            Xru <- ru %*% pnt >= 0
            Xro <- ro %*% pnt >= 0
            alpha <- acos(ru %*% ro)
            betaU <- acos(min(1, abs(ru %*% pnt))) # min only for safety!
            betaO <- acos(min(1, abs(ro %*% pnt))) # min only for safety!
            if (Xru) {
              if (Xro) { # rechts
                loc[[i]]$x[2] <- x$x
                if (betaU < bDa * alpha) loc[[i]]$y[1] <- x$y else
                if (betaO < bDa * alpha) loc[[i]]$y[2] <- x$y
              } else {   # unten
                loc[[i]]$y[1] <- x$y
                if (betaU < bDa * (pi - alpha)) loc[[i]]$x[2] <- x$x else
                if (betaO < bDa * (pi - alpha)) loc[[i]]$x[1] <- x$x
              }
            } else {
              if (Xro) { # oben
                loc[[i]]$y[2] <- x$y
                if (betaU < bDa * (pi - alpha)) loc[[i]]$x[1] <- x$x else
                if (betaO < bDa * (pi - alpha)) loc[[i]]$x[2] <- x$x
              } else {   # links
                loc[[i]]$x[1] <- x$x
                if (betaU < bDa * alpha) loc[[i]]$y[2] <- x$y else
                if (betaO < bDa * alpha) loc[[i]]$y[1] <- x$y
              }
            }
          }
          
          screen(thresh.dev)
          image(1:dp[1], y.coord, 3 * old.picture, col=col.thresh,
                axes=!y1Gy2, zlim=c(-0.5,3.5))
          if (y1Gy2) Axes()
        }
        if (i==1 && length(loc)==1) loc[[3]] <- loc[[2]] <- loc[[1]]
      }
    
      ## to put a moved edge into the position of the ideal edge is
      ## quite hard, so if the line is closer than line.tolerance
      ## to an ideal one it is assumed that the user likes to have
      ## them identical
      
      for (i in 1:2) {
        for (j in 2:3) {
          if (abs(loc[[j]]$x[i] - loc[[1]]$x[i]) < line.tolerance)
            loc[[j]]$x[i] <-loc[[1]]$x[i]
          if (abs(loc[[j]]$y[i] - loc[[1]]$y[i]) < line.tolerance)
            loc[[j]]$y[i] <-loc[[1]]$y[i]
        }
      }
    }
    
   
    options(old.options)
    screen(menue.dev)
    par(mar=rep(0,4))
    plot(1, 1, col="white", axes=FALSE);
    text(1, 1, "calculating risk index ...", col="red")
    
    picture$loc <- loc
  } ## interactive

  if (debug) cat("interactive end\n")
  picture$r.i <- picture$risk.index <-
    picture$raw.risk.index <- NULL

  if (is.null(picture$loc))
    picture$loc <- rep(list(list(x=c(1, dp[1]), y=c(1,dp[2]))), 3)
  else stopifnot(is.list(picture$loc), length(picture$loc)==3,
                 is.list(picture$loc[[1]]))

  ## count the number of blue pixels
  old.picture <-
    old.picture[picture$loc[[1]]$x[1]:picture$loc[[1]]$x[2],
                picture$loc[[1]]$y[1]:picture$loc[[1]]$y[2]]
    
  absfreq <- apply(old.picture, 2, sum, na.rm=TRUE)
  weights <- apply(!is.na(old.picture), 2, sum)
  uw <- unique(weights)
  if (length(uw)!=1 && (length(uw)!=2 || all(uw!=0)))
    warning("number of NAs differs among the lines")
  
  ## correct/estimate number of blue points in case parts 
  ## are not visible
  freq <- absfreq # / weights * nrow(old.picture)
  weights <- sqrt(weights)
  absfreq[weights==0] <- NA
  picture$absfreq <- rev(absfreq)
  
  ## estimation of riskindex.
  if (!is.null(estimate.all)) {
    picture$r.i <-
      risk.index(cbind(length(freq):1, freq),
                         weights=weights, 
                         PrintLevel=PrintLevel,
                         selected.dist=selected.dist,
                         selected.rate = selected.rate,
                         measure=measure, method=method, 
                         endpoint.tolerance=endpoint.tolerance)
        
    if (PrintLevel>0)
      if (is.na(picture$r.i$risk.index)) cat("no shape parameter estimated.\n")
      else if (PrintLevel>1)
        cat("estimated risk index:", picture$r.i$risk.index, "\n")

    if (estimate.all) {
      ## mm[[1:4]] : ranges of x and y coordinates
      ## note that x- (or y-) coordinates are expanded simultaneously
      ##
      ## mm[[>4]]  : ranges of the parameter values for fct
      ## if upper and lower bound are equal mm[[]] contains only 1 element !
      mm <- list()
      loc <- picture$loc
      stopifnot(max(loc[[2]]$x[1], loc[[3]]$x[1]) <
                min(loc[[2]]$x[2], loc[[2]]$x[2]),
                max(loc[[2]]$y[1], loc[[3]]$y[1]) <
                min(loc[[2]]$y[2], loc[[2]]$y[2])
                )
      xy <- c("x", "y")
      ## the left and right edges (and the upper and lower one) are
      ## changed simulatenously
      for (nxy in 1:2) { ## xmin, xmax, ymin, ymax
        for (j in 1:2) {
          idx <- j + 2 * (nxy == 2)
          rg <- range(loc[[2]][xy[nxy]][[1]][j], loc[[3]][xy[nxy]][[1]][j])
          mm[[idx]] <-
            if (diff(rg)) seq(rg[j], rg[3-j],len=picture$ppi.xy[idx]) else rg[1] 
        }
        ## expand upper and lower set of interpolating points to get
        ## equal lengths of the pairs
        idx <- 2 * (nxy == 2) + 1
        rg <- max(length(mm[[idx]]), length(mm[[idx + 1]]))
        mm[[idx]] <- rep(mm[[idx]], len=rg)
        mm[[idx+1]] <- rep(mm[[idx+1]], len=rg)
      }

      nm <-  names(picture$param)
      for (i in 1:length(picture$param)) {     ## rgb parameters
        rg <- range(picture$lower[nm[i]][[1]], picture$upper[nm[i]][[1]])
        mm[[i+4]] <-
           if (diff(rg)) seq(rg[1], rg[2], len=picture$ppi.par[nm[i]][[1]]
                             ) else rg[1]
      }
     
      len.all.intv <- sapply(mm, function(x) length(x))[-2:-3] # due to the
      ## silmutaneous expansion of the x-coordinates (resp. y-coordinates)
      ## only one of the indices 1,2 and 3,4 are needed to calculate
      ## the number of parameters we run through
      picture$risk.index <- picture$raw.risk.index <-
        array(NaN, dim=len.all.intv) 
      prod.len.all.intv <- c(1, cumprod(len.all.intv[-length(len.all.intv)])) 
      len.par.intv <-  len.all.intv[-1:-2] ## delete remaining values for x & y
      total <- prod(len.par.intv)
      np <- length(picture$param)
      param <- numeric(np)
      names(param) <- names(picture$param)
      for (ix in 1:length(mm[[1]])) {
        for (iy in 1:length(mm[[3]])) {
          if (PrintLevel>2)
            cat(sep="","ix=", ix, "(",length(mm[[1]]),
              "), iy=", iy, "(",length(mm[[3]]), ")\n")
          figure <-
            picture$picture[mm[[1]][ix]:mm[[2]][ix], mm[[3]][iy]:mm[[4]][iy], ]

          idx <- rep(1, len=np)          
          for (i in 1:total) {
            for (p in 1:np) param[p] <- mm[[p+4]][idx[p]]
            result <- picture$fct(figure, as.list(param))
            freq <- apply(result, 2, sum, na.rm=TRUE)
            weights <- apply(!is.na(result), 2, sum)
            # freq <- freq / weights * nrow(figure)
            weights <- sqrt(weights)
            
            a <- risk.index(cbind(length(freq):1, freq),
                            weights=weights, 
                            PrintLevel=PrintLevel,
                            selected.rate = selected.rate,
                            measure=measure, method=method, 
                            endpoint.tolerance=endpoint.tolerance
                            )[c("risk.index", "raw.risk.index")]
            if (PrintLevel>1) cat(a[[1]],"\n")
            pos <- 1 + sum( (c(ix, iy, idx)-1) * prod.len.all.intv)
            picture$risk.index[pos] <- a[[1]]            
            picture$raw.risk.index[pos] <- a[[2]]
            
            if (i==total) break
            d <- 1
            while (idx[d] == len.par.intv[d]) {
              idx[d] <- 1
              d <- d + 1
            }
            stopifnot(d <= np)
            idx[d] <- idx[d] + 1
          }
        }
      }
    } # estimate.all
  }
  return(picture)
}




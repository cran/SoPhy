

## todo: modify sink2.f so that different plants type have different
##       P0,...
##       if roots of different plants are located in the same place
##       then take as lower/upper limits the min/max of the values

## atmopheric + root!
## pay attention that integral for all plants over potential uptake equals one
## 


create.stones <- function(h, trials=10){
  if (is.null(h$RF)) {
    warning("h$RF is NULL -- stones are not simulated")
    return(h)
  } 
  h$Stone.RF <- h$RF
  hxlim <- h$grid.x[c(1, length(h$grid.x))]
  hylim <- h$grid.y[c(1, length(h$grid.y))]
  min.size <- h$step / 100

  ## erase former simulation of stones
  idx <- (as.integer(h$idx.rf / h$max.horizons) %% 2) >= 1
  h$idx.rf[idx] <- h$idx.rf[idx] - h$max.horizons
  
#if (length(dev.list())==1) get(getOption("device"))() else dev.set(3)
#image(h$idx.rf,col=rainbow(30))
#readline("press return")
#dev.set(2)       
 
  for (i in 1:h$n) {    
    hor <- i-1
    lambda <- h[[i]]$stone$lambda
    if (!is.null(lambda) && (lambda>0)) {
      p <- rpois(1, lambda=lambda * sum(h$idx.rf == hor) * h$step^2)
      if (p>0) {
        xlim <- (h[[i]]$cut.x - 1 + c(-0.4999,0.4999)) * h$step + h$grid.x[1]
        ylim <- (h[[i]]$cut.y - 1 + c(-0.4999,0.4999)) * h$step + h$grid.y[1]
        n <- 1
        counter <- 0
        p10 <-  p * trials
        while (n<=p) {
          if ((counter <- counter + 1) > p10) {
            warning(paste("Too high required density of stones in horizon", i))
            break
          }
          count <-  0          
          repeat {
            x <- runif(1, min=xlim[1], max=xlim[2]) 
            y <- runif(1, min=ylim[1], max=ylim[2])   
            if (h$idx.rf[1 + round((x - h$grid.x[1]) / h$step),
                         1 + round((y - h$grid.y[1]) / h$step)] %% h$max.horizons
                == hor)
              break;
            if ((count <- count + 1) > trials) {
              warning(paste("Too small area of horizon",i))
              x <- NA
              break
            }
          }
          
          if (is.na(x)) break
          phi <- h[[i]]$stone$phi.distr(1,h[[i]]$stone$phi.mean,
                                        h[[i]]$stone$phi.s) %% pi
          ax <- pmax(min.size,
                     c(h[[i]]$stone$main.distr(1, h[[i]]$stone$main.mean,
                                               h[[i]]$stone$main.s),
                       h[[i]]$stone$sec.distr(1, h[[i]]$stone$sec.mean,
                                              h[[i]]$stone$sec.s)))
          v <- matrix(c(cos(phi), sin(phi), -sin(phi), cos(phi)), ncol=2)
          m <- v %*% diag(1/ax^2) %*% t(v)
          A <- m[1, 1]
          B <- m[1, 2]
          C <- m[2, 2]
          maxmin <- function(x,y) c(max(x[1], y[1]), min(x[2], y[2]))
          ellipse.xlim <- maxmin(hxlim, x + sqrt(1 / (A - B^2 / C)) * c(-1,1))
          ellipse.ylim <- maxmin(hylim, y + sqrt(1 / (C - B^2 / A)) * c(-1,1))
          ellipse.grid.x <-
            1 + as.integer(round((ellipse.xlim - h$grid.x[1]) / h$step))
          ellipse.grid.y <-
            1 + as.integer(round((ellipse.ylim - h$grid.y[1]) / h$step))
          ellipse.grid <-
            as.matrix(expand.grid(ellipse.grid.x[1]:ellipse.grid.x[2],
                                  ellipse.grid.y[1]:ellipse.grid.y[2]))
          e.grid <- t((t(ellipse.grid)-1) * h$step + c(h$grid.x[1], h$grid.y[1])
                      -c(x,y))
          ellipse <- ellipse.grid[apply(e.grid %*% m * e.grid, 1, sum) <= 1,,
                                  drop=FALSE] ## drop=FALSE, otherwise,
          ## if only 1 points is left, it is turned to a vector and the
          ## the access to h$idx.rf[ellipse] is interpreted an access to
          ## a vector!
          ellipse.value <- h$idx.rf[ellipse]
          already.stone <- as.integer(ellipse.value / h$max.horizons) %% 2 >= 1
          if (h[[i]]$stone$no.overlap && any(already.stone)) next
##27.12.03
# cat(h[[i]]$stone$no.lower, x, y);  print(ellipse.value)
# print(ellipse.xlim); print(ellipse.ylim)
# print(ellipse.grid.x); print(ellipse.grid.y)         
          ellipse.value <- ellipse.value %% h$max.horizons
          if (h[[i]]$stone$no.lower && any(ellipse.value < hor)) next
          if (h[[i]]$stone$no.upper && any(ellipse.value > hor)) next
          ellipse <- ellipse[!already.stone, , drop=FALSE] ## comment see above
          h$idx.rf[ellipse] <- h$idx.rf[ellipse] + h$max.horizons        
          n <- n + 1
          
 ##27.12.03
#  cat(i,x,y, 1 + round((x - h$grid.x[1]) / h$step),
#       1 + round((y - h$grid.y[1]) / h$step),":", dim(h$idx.rf), "\n")
#   if (length(dev.list())==1) get(getOption("device"))() else dev.set(3)
#   image(h$grid.x, h$grid.y, h$idx.rf,col=rainbow(30))
#   points(x + e.grid[,1] * h$step, y + h$step * e.grid[,2], pch=".")
#   points(x, y, pch=16)      
#    if (i>1) readline()
#    dev.set(2)

        }
      }
    }
  }
  storage.mode(h$idx.rf) <- "integer"

  stone.values <- unlist(sapply(h[1:h$max.horizons], function(x) x$stone$value))
  st.idx <- as.integer(h$idx.rf / h$max.horizons) %% 2 >= 1 
  h$Stone.RF[st.idx] <- stone.values[1 + h$idx.rf[st.idx] %% h$max.horizons]
  return(h) 
}

create.roots <- function(h, trials=10, PrintLevel=RFparameters()$PrintLevel,
                         message=NULL) {
  ## plants _ list of roots
  ## roots _ array([x, y, prev, next; children, level, dist, knotdist],
  ##               nrow=totallength of roots)
  ## if knot then sht>1 and nxt gives first shoot others following directly
  ## after within the table `root'

  if (is.null(h$Stone.RF)) {
    warning("h$Stone.RF is NULL -- roots are not simulated")
    return(h)
  } 
  h$Root.RF<- h$Stone.RF

  ### see also create.waterflow
  x <- 1     ## coordinates of the root segment
  y <- 2
  prv <- 7   ## index to previous root segment
  nxt <- 4   ## index to the next root segment(s)
  chld <- 5  ## number of subsequent roots segments in case of a knot
  lev <- 6   ## which knot level is it? (keeping track)
  #
  dst <- 3   ## which distance to first root segment (keeping track)
  kntdst <- 8## distance to the preceding knot (keeping track)
  h$plants <- list()
  h$plants.idx <- integer(1)
  lx <- length(h$grid.x)
  ly <- length(h$grid.y)
  xlim <- range(h$grid.x)
  ylim <- range(h$grid.y)
  pl <- 0 ## currently simulated, pl-th plant
  n.plants <- 0

  for (jj in 1:length(h$root)) { ## run through all plant types
    n.p <- rpois(1, h$root[[jj]]$plants.lambda * lx)
    n.plants <- n.plants + n.p
    if (n.p==0) next
    start <- numeric(n.p)
    for (plx in 1:n.p) {
      pl <- pl + 1
      h$plants.idx[pl] <- jj
      for (i in 1:(trials)) {
        ## if could not find a good place for trials times
        ## then put that plant anywhere (taking over last simulate x-coordinate)
        start[plx] <-
          1 + as.integer(round((runif(1, xlim[1] - h$step / 2,
                                      xlim[2] + h$step / 2) -
                                h$grid.x[1]) / h$step))
        stopifnot(start[plx] >= 1, start[plx]<=nrow(h$Stone.RF))
        ##                               to check whether programming is fine 
        top <- max(0, which(!is.nan(h$Stone.RF[start[plx], ])))
        h$plants[[pl]] <-
          t(c(x=start[plx],
              y=top+1,
              dst=rnorm(1, m=h$root[[jj]]$mean, sd=h$root[[jj]]$sd),
              ##                       aimed distance in first row!
              nxt=2,
              "#chldn"=1,
              lev=0, 
              prv=NA,
              kntdst=0
	    ))
        ## plant on stone may happen, but will not grow !!!
        if ((plx==1) || 
            all(abs(start[plx] - start[1:(plx-1)]) * h$step
                >= h$root[[jj]]$plants.mindist))
          break     
      }
      
      if (h$plants[[pl]][1, dst]<=0) next      ## if value of aimed distance
      ##                                          is negative
      if ((h$plants[[pl]][1, dst] < h$step) && 
	 (runif(1) > h$plants[[pl]][1, dst]/h$step)) next ## aimed distance
      ##                                             is less than one pixel

      ## add the next point vertically below if the plant has not started
      ## on the top of a
      
      if (top>0 &&  ## if only air in that  
          (!is.na(h$Stone.RF[h$plants[[pl]][1, x], top]) || 
           is.nan(h$Stone.RF[h$plants[[pl]][1, x], top])))
        h$plants[[pl]] <-                      
          rbind(h$plants[[pl]],
                c(h$plants[[pl]][1, x], top, h$step, NA, NA,
                  0,  1,
		  h$step + h$root[[jj]]$knot.mindist  ## 
		  ## oder doch lieber " + h$root[[jj]]$knot.mindist" weglassen?
                  ## dann haben alle Wurzelstoecke einen Hals
		  ))
    }
  }
  
  if (PrintLevel>5) print(h$plants)
  xtrim <- function(x) {
    x <- (x-1) %% (2 * lx)
    if (x>lx) 2 * lx - x + 1 else x + 1
  }
  ytrim <- function(y) { stopifnot(y <= ly); max(1, y) }

  if (n.plants > 0) for (pl in 1:n.plants) {

    if (!is.null(message))
      message(paste("generating roots:", pl, "out of", n.plants, "plants"))
    if (PrintLevel>3)
      cat("\n", pl, "out of", n.plants, "plants; step=", h$step,"\n")
    
    root <- h$plants[[pl]]
    len <- root[1,dst]      ## aimed distance
    if (PrintLevel>5)   print(c(pl,len))  
    k <- 2
    if ((nrow(root)==1) || (root[k, dst]>=len)) next
    jj <- h$plants.idx[pl]
    ## tab: 1st col: index to segment; 2nd col:value
    direction <- root[k,y] - root[root[k,prv], y]
    tab <-
      cbind(k, -h$root[[jj]]$depth.bonus * (direction < 0) 
            - h$root[[jj]]$side.bonus * (direction == 0) +
            h$root[[jj]]$rf.link(h$Stone.RF[root[k, x], root[k, y]],
                                 h$root[[jj]]$rf.m, h$root[[jj]]$rf.s)
            )
    j <- 3
    total <- h$step
    while (nrow(tab)>0) {
     if (PrintLevel>2) cat(j, nrow(root), "\b\b\b\b\b\b\b\b\b\b")      
      k <- tab[1,1]
      tab <- tab[-1,,drop=FALSE]
      if (runif(1) <
         h$root[[jj]]$stop.probab(root[k, lev], root[k, dst],
                                  h$root[[jj]]$stop.m, h$root[[jj]]$stop.s)) {
         next
      }
      if (nrow(tab)>0) tab[,2] <- tab[,2] * h$root[[jj]]$age.bonus


      ## neighbours of selected root segment
      
      p <- root[k,c(x,y)] +
        if (h$root[[jj]]$diagonal)
          t(as.matrix(expand.grid(c(-1, 0, 1), c(-1, 0, 1))))[, -5] # not (0,0)
        else rbind(c(-1, 1, 0, 0), c(0, 0, 1, -1))
      p <- p[, (p[1,]>=1) & (p[1,]<=lx) & (p[2,]>=1) & (p[2,]<=ly), drop=FALSE]
      

      ## no stone, not air (coded by NA)
      p <- p[, !is.na(h$Stone.RF[t(p)]), drop=FALSE]
      if (ncol(p)==0) next

      ## no own root
      if (h$root[[jj]]$no.own.root)
        p <- p[, apply(outer(root[1:k, x], p[1,], "!=") |
                       outer(root[1:k, y], p[2,], "!="), 2, all), drop=FALSE]
      if (ncol(p)==0)  next     
      
      ## direction change
      v <- - h$root[[jj]]$dir.ch * apply(abs(2 * root[k,c(x,y)] -
                                         root[root[k,prv], c(x,y)] - p), 2, sum)
      v <- rnorm(length(v), v, h$root[[jj]]$dir.ch.s * abs(v))

      ## depth.bonus
      v <- v + h$root[[jj]]$depth.bonus * (root[k, y] - p[2, ] == 1) +
        h$root[[jj]]$side.bonus * ((root[k, y] - p[2, ] == 0) ) +
            runif(length(v), max=1E-7) ## some small randomness

      ## random field
      v <- v + h$root[[jj]]$rf.link(h$Stone.RF[t(p)],
                                    h$root[[jj]]$rf.m, h$root[[jj]]$rf.s)
      ## the following three lines are necessary in case
      ## there are only NA or only Inf or only -Inf. Then `order'
      ## would leave the ordering unchanged!
      idx <- sample(1:length(v))
      v <- v[idx]
      p <- p[, idx, drop=FALSE]
      
      o <- order(v, decreasing=TRUE)
      p <- p[, o, drop=FALSE]
      v <- v[o]

      ## knot ?h$
      if (knot <- (h$root[[jj]]$knot.mindist<root[k, kntdst]) &&
          (runif(1) < h$root[[jj]]$knot.probab)) {
        root[k, chld] <-
          min(length(v),
              4 - sum(runif(1) >
                      cumsum(c(h$root[[jj]]$shoots.4, h$root[[jj]]$shoots.3))))
        level <-  root[k, lev] + 1
      } else {
        root[k, chld] <- 1
        level <- root[k, lev]       
      }
        
      root[k, nxt] <- j

      ## calculate the added total lengths
      ## restrict the number of children if the aimed total length
      ## is exceeded
      delta <- sqrt(apply((p[,1:root[k, chld], drop=FALSE] - root[k,c(x,y)])^2,
                          2, sum)) * h$step
      if ((sdelta <- sum(delta)) + total < len) {
        total <- total + sdelta
        n <- root[k, chld]
      } else {
        for (n in 1:root[k, chld]) {
          if (delta[n] + total > len) {
            if (runif(1) * delta[n] >= len - total) n <- n - 1	    
            break  
	  }
          total <- total + delta[n]
        }
        if (n==0) break
	delta <- delta[1:n]
      }
      p <- t(p[,1:n])
      v <- v[1:n]  
      root <- rbind(root,
                    cbind(p, root[k, dst] + delta, NA, NA, level, k,
                          (1 - knot) * root[k, kntdst] + delta))
      tab <- rbind(tab, cbind(j:(j+n-1), v))
      j <- j + n
      if (nrow(tab)>1) tab <-tab[order(tab[,2], decreasing=TRUE),,drop=FALSE]
    }
    h$plants[[pl]] <- root
  }

  if (is.function(h$root.zone)) {
    n.plants <- n.plants + 1
    pl <- pl + 1
    h$plants[[pl]] <- rep(NA, 3) # c(1, 2, 7)
    h$plants.idx[pl] <- 1
    root.zone <- matrix(FALSE, nrow=length(h$grid.x), ncol=length(h$grid.y))
    for (i in 1:length(h$grid.x)) {
      ylim <- h$root.zone(h$grid.x[i]);
      ylim <- which(h$grid.y>=ylim[1] & h$grid.y<=ylim[2])
      if (length(ylim)>0) { 
        h$plants[[pl]] <-
          rbind(h$plants[[pl]],
                cbind(i, ylim, (length(h$grid.y) - ylim) * h$step))
      }
    }
  }
  
  
  if (length(h$plants)>0) {
    ## 22.6.03, f muss doch immer 1 sein, oder?! Reduziert wird doch erst
    ## spaeter
    ## f <- max(1, h$water$red)
    f <- 1 ## 22.6.03
    for (i in 1:length(h$plants)) { ## lahme Programmierung!
      if (h$plants.idx[i]>0 && h$root[[h$plants.idx[i]]]$rf.Kf) {
        red.pl <- as.integer((h$plants[[i]][-1, x, drop=FALSE] - 1) / f + 1) +
          nrow(h$RF) * as.integer((h$plants[[i]][-1, y, drop=FALSE] - 1) / f) 
        h$Root.RF[red.pl] <- 
          h$Kf(plant=h$plants.idx[i],
               distances=h$plants[[i]][-1, dst],
               depth=(length(h$grid.y) - h$plants[[i]][-1, y] + 1) * h$dist,
	       rf=h$Stone.RF[h$plants[[i]][-1, c(x,y)]],
               param=h$root[[h$plants.idx[i]]])
      }
    }
  }
  return(h)
}


simulate.horizons <- function(h, first=1, what='all',
                              PrintLevel=RFparameters()$Print,
                              message=function(s) {print(s)}, stone.trials=50
                              ) {
  ## returns a Gaussian random field !
  what.options <- c("randomfield", "stone", "root", "all", "auto")
  RandomField <- 1 	
  Stone <- 2       
  Root <- 3
  All <- 4
  Auto <- 5

#  str(what); str(what.options)
  
  if (any(is.na(what <- pmatch(what, what.options)))) {
    if (PrintLevel>0) message("unmatched value for 'what'")
    return(h)
  }    
  if (All %in% what) what <- c(RandomField, Stone, Root, All) else 
  if (Auto %in% what) {
     what <- c(what, 
               if (is.null(h$RF)) c(RandomField, Stone, Root) else
               if (is.null(h$Stone.RF)) c(Stone, Root) else
	       if (any(c(is.null(h$Root.RF), is.null(h$plants.idx), 
	           is.null(h$plants)))) Root)
  }
  
  if (All %in% what || is.null(h$random.seed)) {
    h$random.seed <- list()
    for (i in 1:h$max.horizons) {
      runif(1)
      h$random.seed[[i]] <-
        get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
    }
    runif(1)
    h$stone.random.seed <-get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
    runif(1)
    h$root.random.seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
  }
  if (Root %in% what && !(Stone %in% what) && is.null(h$Stone.RF)) {
    if (PrintLevel>0) message("stones will also be generated")
    what <- c(what, Stone)
  }
  if (Stone %in% what && !(RandomField %in% what) && is.null(h$RF)) {
    if (PrintLevel>0) message("random field will also be generated")
    what <- c(what, RandomField)
  }
  
  if (is.null(h[[1]]$model)) {
    h$RF <- h$Stone.RF <- h$Root.RF <- NULL
    if (PrintLevel>0)
      message("no simulation: structure of main (first) horizon undefined")
    return(h)
  }

  if (RandomField %in% what) {
    rf <- matrix(NA, nrow=length(h$grid.x), ncol=length(h$grid.y))
    message("generating random field ...")
    if (is.null(h$RF)) first <- 1 else
    if (first>1) {
      for (i in 1:(first-1)) {
        if (is.null(h[[i]]$model)) {
          message(paste("structure of horizon #", i, "undefined!"))
          return(h)
        }
        cx <- h[[i]]$cut.x[1] : h[[i]]$cut.x[2]
        cy <- h[[i]]$cut.y[1] : h[[i]]$cut.y[2]
        if (i==1) rf[cx, cy] <- h$RF[cx, cy]
        else rf[cx, cy][h[[i]]$idx] <- h$RF[cx, cy][h[[i]]$idx]         
      }
    }


    h$rf.complete <- TRUE
    for (i in first:h$n) {
      assign(".Random.seed", h$random.seed[[i]], envir=.GlobalEnv)
      if (is.null(h[[i]]$model)) {
        mssge <- paste("structure of horizon #", i,
                       "undefined -- everything skipped afterwards")
        for (j in i:h$n) {
          cx <- h[[j]]$cut.x
          cy <- h[[j]]$cut.y
          rf[cx[1]:cx[2], cy[1]:cy[2]][h[[j]]$idx] <- NA
        }
        h$rf.complete <- FALSE
        break;
      }

      cx <- h[[i]]$cut.x
      cy <- h[[i]]$cut.y
      
      ##  "timedim=2" should be improved
      ## 21.6.03 is.nan changed to is.na
      if (any(is.na(unlist(PrepareModel(model=h[[i]]$model, timespacedim=2)
                           [c("param","mean")])))){
        if (h[[i]]$type!="H") { ## h[1]$type!="Start" !
          message("First horizon cannot be NaN (air)")
          return(h)
        }
        if (i < h$n && any(sapply(h[(i+1) : h$n], function(x) x$type) != "P")) {
          message("Only the last horizon can be NaN (air).")
          return(h)
        }
        
        rf[cx[1]:cx[2], cy[1]:cy[2]][h[[i]]$idx] <- NaN
        next
      }
      
      if (i==1) {
        rf[cx[1]:cx[2], cy[1]:cy[2]] <-
          GaussRF(x=h$grid.x[cx[1]:cx[2]], y=h$grid.y[cy[1]:cy[2]],
                              grid=TRUE, model=h[[i]]$model, register=i)
      } else {
        border <- h[[i]]$border
        if (any(select <- is.na(rf[border]))) {
          ## should only happen if horizons are crossing/cutting !          
          if (PrintLevel>1) {
            cat("horizon #", i,": some NA border points\n")
            if (PrintLevel>3) print(border[select,])
          }
          border <- border[!select, ]
        }
        
        ## standard case, not NaN
        mn <- h[[i]]$model$mean
        h[[i]]$model$mean <- 0
        
        rf[cx[1]:cx[2], cy[1]:cy[2]][h[[i]]$idx] <-
          mn + if (h[[i]]$materials$sharpness<1) {
            sqrt(1 - h[[i]]$materials$sharpness^2) *
              CondSimu("S", x=h$grid.x[cx[1]:cx[2]],
                       y=h$grid.y[cy[1]:cy[2]],
                       grid=TRUE, model=h[[i]]$model, register=i,
                       given=cbind(h$grid.x[border[,1]],
                         h$grid.y[border[,2]]),
                       data=rf[border] - mn, pch="", na.rm=TRUE
                       ) [h[[i]]$idx] +
                         h[[i]]$materials$sharpness *
                           GaussRF(x=h$grid.x[cx[1]:cx[2]],
                                   y=h$grid.y[cy[1]:cy[2]],
                                   grid=TRUE, model=h[[i]]$model,
                                   register=i)[h[[i]]$idx]
          } else {
            (zz <- GaussRF(x=h$grid.x[cx[1]:cx[2]], y=h$grid.y[cy[1]:cy[2]],
                           grid=TRUE, model=h[[i]]$model,
                           register=i))[h[[i]]$idx]
          }
        h[[i]]$model$mean <- mn
      }
    }
    h$RF <- rf
    h$Stone.RF <- h$Root.RF <- NULL
  }
  
  if (Stone %in% what) {
    assign(".Random.seed", h$stone.random.seed, envir=.GlobalEnv)
    message("generating stones ...")
    stone.values <-
      unlist(sapply(h[1:h$max.horizons], function(x) x$stone$value))
    if (any(is.nan(stone.values))) {
      message(paste("Forbidden stone$value for the",
                    which(is.nan(stone.values)), "th horizon/polygon", sep=""))
      h$Stone.RF <-  h$Root.RF <- NULL
      return(h)
    }
    h <- create.stones(h, trials=stone.trials)
    stopifnot(!is.null(h$Stone.RF))
    h$Root.RF <- NULL
  }
  
  if (Root %in% what) {
    assign(".Random.seed", h$root.random.seed, envir=.GlobalEnv)
    message("generating roots ...")
    h <- create.roots(h, message=message)
    stopifnot(!is.null(h$Root.RF))
  } 

  message("")
  return(h)
}

create.waterflow <-
  function(h,
           ## h$atm.data might be NULL
           #Units = c(LUnit = 'cm', TUnit = 'sec',  MUnit = '-',  BUnit = '-'),
           Kat = 2, # type of flow system; 0:horizontal, 1:axisymm., 2:vertical
           MaxIt = 20, # max # of iterations/time step; usually 20
           hTab = c(0.001,200),# interval within which table values used (ha,hb)
           dtMinMax = h$water$dtMinMax, # min and max permitted time increment
           DMul = c(1.1, 0.33),# iterations<=3 then time step multipl.d by dMul1
                               # iterations>=7 then time step multipl.d by dMul2
           atmadd=FALSE
           )
{
  ## no drain included yet -- keep it simple...
  ##     Nobs = 0
  root.x <- 1
  root.y <- 2
  root.dst <- 3 ## see create.roots

  seepage <- list()
  stopifnot(length(seepage)==0) ## huwe fragen unterschied zu free drainage und drains !

  NSeepD <- 2    #
  NumSPD <- 50
  NDrD <- 2
  NElDrD <- 8
  NMatD <- 20
  NumKD <- 6   ## always >= 5 !!
  NObsD <- 4

  Node <- NULL
   
  DrainF <-  FALSE ## !!! no drain included yet -- keep it simple...
  NDr <- DrCorr <- ND <- NED <- EfDim <- KElDr <- NA  
  NMat <- h$n

  lx <- length(h$grid.x)
  max.y <-
    ly <- length(h$grid.y)
  f <- max(1, h$water$red)
  water.x <- h$grid.x[water.idx.x <- seq(1, lx, f)]
  water.y <- h$grid.y[water.idx.y <- rev(seq(ly, 1, -f))]
  lx <- length(water.x)
  ly <- length(water.y)
  NumNP <- lx * ly # number of total (possible) nodes
  if (is.null(h$Root.RF)) return("error: is.null(h$Root.RF)")
  rf <- h$miller.link(h$Root.RF[water.idx.x, water.idx.y, drop=FALSE])
  if (any(dim(rf)<2)) return("error: grid too coarse")

  horizon <-
    as.vector((h$idx.rf[water.idx.x, water.idx.y] %% h$max.horizons) + 1)
  
  Par <- ## materials and Par are the same 
    t(sapply(h[1:h$n], function(x)
             unlist(x$materials[c("thr", "ths", "tha", "thm",
                                  "Alfa", "n", "Ks", "Kk", "thk")])))
  
  mat.other <-
    t(sapply(h[1:h$n],
             function(x) unlist(x$materials[c("Hseg", "Hslope", "POptm")])))
  stopifnot(is.matrix(Par), is.matrix(mat.other))
  HStart <- mat.other[, 1:2, drop=FALSE]
  POptm <- mat.other[, 3, drop=FALSE]

  NLay <- h$n
  NPar <- 9  ## only van Genuchten model used,

  if (AtmInf <- h$atmosphere$AtmInf){
    if (h$water$TPrint[1] < h$atmosphere$tInit + h$water$dt)
      return("error: 1st TPrint < atmosphere$tInit + water$dt")
 
    if (!is.matrix(atm.data <- h$atm.data) || (ncol(atm.data)!=10))
      return("error: atmosphere not a matrix of 10 columns")
    dummy <- diff(c(h$atmosphere$tInit, atm.data[,1]))
    if (!all(is.finite(dummy)) | any(dummy<=0))
      return("error: atmosph. time increments are not all positive")
    if (atm.data[nrow(atm.data), 1] < h$water$TPrint[length(h$water$TPrint)])
      return("error: last endpoint of atmospherical periods must be behind TPrint, see swms2d (water)")
    GWL0L <- h$atmosphere$GWL0L
    Aqh <- h$atmosphere$Aqh
    Bqh <- h$atmosphere$Bqh
    tInit <- h$atmosphere$tInit
    hCritS <- h$atmosphere$hCritS
  } else {
    atm.data <- NULL
    GWL0L <- Aqh <- Bqh <- hCritS <- NA
    tInit <- 0
    if (h$water$TPrint[1] < tInit + h$water$dt)
      return("error: 1st TPrint  <  water$dt")
  }
   
  stopifnot(ncol(Par)==NPar,
            nrow(Par)==h$n,
            (length(h$plants)==0) || !AtmInf ||  (nrow(Par)==length(POptm)))
  stopifnot(length(hTab)==2, length(dtMinMax)==2, length(DMul)==2,
            NSeepD>=length(seepage))
  NSP <- NP <- NSeep <- NA
  if (SeepF <- (NSeep <- length(seepage)) > 0) {
    NSP <- len(unlist(lapply(seepage, length)), NSeepD) 
    NP <- matrix(NA, nrow=NSeepD, ncol=NumSPD)
    for (i in 1:length(seepage)) {
      l <- length(seepage[[i]])
      stopifnot(l<=NumSPD) 
      NP[i, 1:l] <- seepage[[i]]
    }
  } 
  stopifnot(NSeep<=NSeepD)         

  hNew <-
    h$Hinit(HStart[horizon, 1],
            rep(h$grid.y[length(h$grid.y)] - water.y, each=length(water.x))) +
              HStart[horizon, 2] * h$millerH(rf)

  
  aK <- as.vector(h$millerK(rf))

  # determine the surface points
  warn.orig <- options()$warn
  options(warn=-1)
 
  surface.idx <-
    (1 : nrow(rf)) + nrow(rf) * 
      (apply(rf, 1, function(x) max(1, which(!is.nan(x)))) - 1)
  options(warn=warn.orig)
  
  ## 2.1.04 : checked: suface.idx always finite!
  ##if (!all(is.finite(surface.idx)))
  ##  return("columns of NaN in profile definition not allowed")


  ## boundary conditions
  FreeD <- FALSE
  qGWLF <-  FALSE
  bottomQ <- topQ <- 0
  topK <- bottomK <- 0
  bottom.list <- c("free drainage", "deep drainage", "impermeable",
                   "H constant (Dirichlet)", "Q constant (Neumann)")
  #, bottom.list, nomatch = length(bottom.list)+1),
  switch(h$water$bottom.bound,
         {bottomK <- -3; FreeD <- TRUE},
         {bottomK <- -3; qGWLF <- TRUE;
          if (!AtmInf) return("error: not available (atmospherical data missing)")
        },
         {bottomK <- 0},
         {bottomK <- 1;  hNew[,1] <- h$water$bottom.value},
         {bottomK <- -1; bottomQ <- h$water$bottom.value},
         {return("error: no match for bottom bound");}
         )
  if (!qGWLF) Aqh <- Bqh <- 0

  
  top.list <- c("impermeable", "H constant (Dirichlet)", "Q constant (Neumann)",
                "variable  H", "variable Q", "atmospheric")

   switch(h$water$top.bound,
         {topK <- 0},
         {topK <- 1;  hNew[surface.idx] <- h$water$top.value},
         {topK <- -1; topQ <- h$water$top.value},
         {topK <- 3
          if (!AtmInf)
            return("error: not available (atmospherical data missing)")},
         {topK <- -3
          if (!AtmInf)
            return("error: not available (atmospherical data missing)")},
         {topK <- -4
          if (!AtmInf)
            return("error: not available (atmospherical data missing)")},
         {return("error: no match for top bound");}
         )
  if (bottomK==topK) return("error: bottom and top conditions do not match")
 
  Kode <- c(rep(bottomK, lx), rep(0, lx * (ly - 1)))
  Q <- c(rep(bottomQ, lx), rep(0, lx * (ly - 1)))
  beta <- rep(0, NumNP)
  ## see after root uptake for top values of Kode and Q


################################################################################
################################################################################
  ## root uptake
  chem.root <- rep(NA, NumNP)
#  if(length(h$root)!=1) warning("more than 1 kind of root not programmed yet")

  SinkF <- FALSE
  rLen <- 0
  atm.idx <- NA

  if (length(h$plants)>0) {

    # get basic data out of h
    plts <- NULL
    plants.idx <- NULL
    for (i in 1:length(h$plants)) {    
      ## pl : coordinates of the roots
      ## 1 : x coordinates
      ## 2 : y coordinates
      ## 3 : distance to first root segment
      plts <- rbind(plts, h$plants[[i]][-1, c(root.x, root.y, root.dst),
                                        drop=FALSE])

      ## plants.idx : idx of the type of plant
      plants.idx <- c(plants.idx, rep(h$plants.idx[i], nrow(h$plants[[i]]) - 1))
    }
    ## reduce grid coordinates according to `reduction'
    idx <- ((((plts[,1] - 1) %% f) == 0) &
            (((length(h$grid.y) - plts[,2]) %% f) == 0)) 
    plts <- plts[idx, , drop=FALSE]
    plants.idx <- plants.idx[idx]

    ## note index start is now zero
    pl <- cbind((as.integer(plts[, 1] - 1) / f), ## as.int only if neumann
                ly - 1 - as.integer((length(h$grid.y) - plts[, 2]) / f),
                plts[, 3]) ## third column needed for .C call
        
    ## rows: variables, cols: materials
    root.mat <- sapply(h$root, function(x)
       unlist(x[c("P0", "P2H", "P2L", "P3", "r2H", "r2L")]))
    root.mat[1:4, ] <- -abs(root.mat[1:4, ]) 
    ## root.cond : 1:"dirichlet", 2:"neumann", 3:"atmospheric") for each type
    root.cond   <- sapply(h$root, function(x) x$root.condition)
    ## value for each type
    root.uptake <- sapply(h$root, function(x) x$root.uptake)
    ## unnormed beta, and for all root particles idependently of root.cond    

    pot.beta <- rep(NA, length(plants.idx))

    
    for (i in 1:length(h$root)) if (any(p.i <- plants.idx==i & pl[, 2]>0)) {
      ##                                                wachsen u.U. unten raus
      pot.beta <- h$beta(plant=i, 
                         distances=pl[p.i, 3], 
                         depth=(ly - pl[p.i, 2]) * h$step,
			 rf=h$Stone.RF[plts[,1:2, drop=FALSE]],
                         param=h$root[[i]]
                         ) # depth
    }

    if (any(is.na(pot.beta))) {
      str(h)
      print(pot.beta)
      if (f==1) # water$reduction
        return("error: NAs in potential water uptake")
      pot.beta[is.na(pot.beta)] <- 0
    }
    if (any(pot.beta<0)) return("error: negative values not allowed for beta")

    if (FALSE) {
      print("OK")
      print(plants.idx)
      print(cbind(pl[,1] + pl[,2] * lx, plants.idx) )
      print("pl")
      print(pl)
      print(idx)
      print("(nrow(pl)")
      print(nrow(pl))
      print(dim(pl))
      print(lx)
      print(cbind(pl[,1] + pl[,2] * lx, plants.idx))
    }
    
    prep <- .C("preproot", 
               as.double(root.mat), as.integer(root.cond), 
               as.double(root.uptake), as.integer(length(root.cond)),
	       as.double(pot.beta),
               ## indices pl[...] start from zero !
               as.integer(cbind(pl[,1] + pl[,2] * lx, plants.idx)), 
	       as.integer(nrow(pl)),
               as.integer(NumNP),
               as.integer(c(2,0,1)), ## defines the priority of the
               ## the conditions: atm (2) has hightest priority,
               ## i.e., e.g., if there are 3 roots of atm condition
               ## and 3 roots of neumann condition and 2 roots of dirichlet
               ## condition then the node get the atm condition.
               ## it would get the dirichlet condition if there were 4 roots
               ## of dirichlet conditions in our example.
               ## output:
               dirich.idx = integer(NumNP),
               neumann.idx = integer(NumNP),
               atm.idx = integer(NumNP),
               dirich.val = double(NumNP),
               neumann.val = double(NumNP),
               ## atm.val is the sum of the raw potential root uptakes
               atm.val = double(NumNP),
               root = double(NumNP * 6), ## liegend
               atmadd = as.integer(atmadd), 
               error = integer(1),
               NAOK=TRUE, DUP=FALSE, PACKAGE="SoPhy")
    if (prep$error!=0) {
      return(paste("error nr.", prep$error," in 'preproot' occured"));
    }
    root <- t(matrix(prep$root, nrow=6))
    
    # dirichlet
    idx <- as.logical(prep$dirich.idx)
    Kode[idx] <- 1
    hNew[idx] <- prep$dirich.val[idx] ## == min(root.uptake)
    chem.root[idx] <- 6

    # neumann
    idx <- as.logical(prep$neumann.idx)
    Kode[idx] <- -1
    Q[idx] <- prep$neumann.val[idx] * aK[idx] ## == multiple of root.uptake
    chem.root[idx] <- -6
    
    # atmosphere
    atm.idx <- as.logical(prep$atm.idx) & (prep$atm.val>0)
    if (!any(atm.idx)) atm.idx <- NA
    else if ((s <- sum(prep$atm.val)) > 0) {
      if (!AtmInf) return("error: !AtmInf and condition=atmosphere")
      beta[atm.idx] <- prep$atm.val[atm.idx]  ## / s
      ## devision by not necessary since done by the Fortran-Code
      rLen <- h$atmosphere$rLen
      SinkF <- TRUE
      ##      chem.root[atm.idx] <- 0 ## leave it as it is...
    }
    
    ## check
    stopifnot(all(prep$dirich.idx + prep$neumann.idx + prep$atm.idx <= 1))

  #  red.pl <- as.integer((pl[, 1, drop=FALSE] - 1) / f + 1) +
  #    lx * as.integer((pl[, 2, drop=FALSE] - 1) / f) ## works!

    ## should be full length, independently whether root is active or not
    POptm <-  POptm[horizon]

    
    if (FALSE)
     { 
       print(POptm)
       print("prep$dirich.idx")
       print(prep$dirich.idx)
       print(prep$dirich.val)
       print("prep$neumann.idx")
       print(prep$neumann.idx)
       print(prep$neumann.val)
       print("prep$atm.idx")
       print(prep$atm.idx)
       print(prep$atm.val)
       print(prep$dirich.idx | prep$neumann.idx | prep$atm.idx)
       readline()
     }
  } else root <- NULL

  ##########################  top values -- overwrite root conditions
  Kode[surface.idx] <- topK
  Q[surface.idx] <- topQ
  if (length(atm.idx)>1 || !is.na(atm.idx))
    Kode[atm.idx] <- pmin(0, Kode[atm.idx])
  ## Do not change values of Kode after this line -- atm.idx must be the last
  
  ##########################  create finite elements
  x <- rep(1 : (lx - 1), ly - 1) + as.vector(outer(rep(lx, lx-1), 0:(ly-2)))
  x <- cbind(x, x + 1, x + lx + 1, x + lx)
  ## stones are coded by NA !!
  horizonZ <- matrix(is.na(rf[as.vector(x)]), ncol=4)
  idx <- apply(horizonZ, 1, sum)
  idx1 <- idx == 1
  
  ## eliminate the na in each row (but we do not know in which column it is)
  x1 <- as.vector(t(x[idx1, ])) ## 4 row, in each column there is 1 NA
  x1 <- t(matrix(x1[as.vector(!t(horizonZ[idx1,]))], nrow=3))

  x[idx1, ] <- cbind(x1, x1[,3]) ## only triples, so duplicate last point
  x <- x[idx <= 1, ] ## meshes with two or more stones in the corner are
  ##                    eliminated
  
  horizonZ <- matrix(horizon[as.vector(x)], ncol=4)
  minhor <- apply(horizonZ, 1, min)
  ## triangular meshes may not be treated further (idx1 reduced by idx<=1, cp x)
  slct <- (apply(horizonZ==minhor, 1, sum) == 1) && !idx1[idx <= 1]
  if (any(slct)) {
    horizonZ <- horizonZ[slct,,drop=FALSE]    
    min2hor <- apply(horizonZ, 1, function(x) min(x[x!=min(x)]))
    x2 <- x[slct, , drop=FALSE]
    idx <-  which(minhor[slct] == horizonZ, arr.ind=TRUE)[,2]
    l2 <-  nrow(horizonZ)
    xP1 <- x2[cbind(1:l2, (idx %% 4) + 1)]
    xP3 <- x2[cbind(1:l2, ((idx+2) %% 4) + 1)]
    x <- rbind(x[!slct,],
                cbind(x2[cbind(1:l2,idx)], xP1, xP3, xP3),
                cbind(xP1, x2[cbind(1:l2, ((idx+1) %% 4) + 1)], xP3, xP3)
                )
    minhor <- c(minhor[!slct], minhor[slct], min2hor)
  }

 
  ##########################  isolated grid points (not part of an element)
  flux <- rep(FALSE, length(rf)) ## note rf is a matrix.
  flux[x] <- TRUE
  newidx <- rep(NA, NumNP)
  NumNP <- sum(flux)
  newidx[flux] <- 1:NumNP
  Kode <- Kode[flux]

  if (!is.null(root)) {
    root <- root[flux, ]
  }

  geometry <- t(sapply(h[1:h$n], function(x)
                       unlist(x$materials[c("angle", "first", "second")])))
                     
  KXR <- cbind(matrix(newidx[x], ncol=4), geometry[minhor, ], minhor)
  KXB <- (1:NumNP)[(Kode!=0)]  
  KodCB <- (chem.root[(Kode!=0)])[!is.na(KXB)]
  KXB <-  KXB[!is.na(KXB)]

  Width <- 0.5 * h$step * h$water$red ## is this a good choice when top boundary
  ##  see manual of swms2d
  ##  is not a horizontal line --- unclear!
  ## NOTIZEN :
  ## Effective Laenge der
  ## elemente neu berechnen, abhaengig davon ob horizontal (1) order
  ## schraeg (sqrt(2)) oder senkrecht (0) oder ueberhaengend(0)
  ##
  ## Idee: nur der alleroberste Punkt einer Spalte hat 
  ## Atmosphaeren-Kontakt, diesen Punkt suchen. Somit existieren genau
  ## length(x) OF-Punkte bis auf Steine.
  ## Achtung OF-finite Elemente koennen u.U. nur einen  OF-Punkt
  ## enthalten, z.b. falls Horizongrenze senkrecht verlaeuft:
  ## oberster Punkt OF-Punkt, aber unterster Punkt der senkr. Strecke
  ## gehoert zum finiten Element
  ##
  ##       O ++++   O:Oberflaechenpkt, +:atmosph. element
  ##       |
  ##    -O+*
  ##     | |
  ## ABER
  ##       O ++++   O:Oberflaechenpkt, +:atmosph. element
  ##      /|
  ##    -O-*
  ##     | |

  
  ## solute transport  
  if (lChem <- h$chem$lChem) {
    Epsi <- c(0, 0.5, 1)[h$chem$Epsi]
    lUpW <- c(TRUE, FALSE)[h$chem$lUpW]
    PeCr <- if (lUpW) 0 else h$chem$PeCr
    ChPar <-
       t(sapply(h[1:h$n], function(x)
                unlist(x$materials[c("Bulk.d", "Diffus", "longiDisper",
                                     "transvDisper", "Adsorp", "SinkL1",
                                     "SinkS1", "SinkL0", "SinkS0")])))
    switch(h$chem$bottom.bound,
           {bottomC <- -3},
           {bottomC <- -3;
            if (!AtmInf)
              return("error: not available (atmospherical data missing)")},
           {bottomC <- 0},
           {return("error: no match for bottom bound")}
         )
     switch(h$chem$top.bound,
         {topC <- 0},
         {topC <- 1},
         {topC <- -1},
         {topC <- 3},
         {topC <- -3},
         {topC <- -4},
         {return("error: no match for top bound")}
         )
    if ((topC==bottomC) && (h$chem$top.value!=h$chem$bottom.value))
      return("error: concurrent definition on top and bottom")

    ### TODO: include h$chem$root.uptake and h$chem$intern.source
    cBound <- rep(0,6)
    if (bottomC!=0) cBound[abs(bottomC)] <- h$chem$bottom.value
    if (topC!=0) cBound[abs(topC)] <- h$chem$top.value

    KodCB <- c(rep(bottomC,lx), rep(0, lx * (ly - 1)))
    KodCB[surface.idx] <- topC
    KodCB <- KodCB[flux]
    KodCB <- KodCB[Kode!=0]
    stopifnot(length(KodCB)==length(KXB))
  } else {
    Epsi <- lUpW <- PeCr <- ChPar <- KodCB  <- NA
    cBound <- rep(NA, 6)
  }
  
  z <- rep(water.y, each=lx)[flux]
   
 ###################
 if (FALSE)
  {
  print(1:NumNP)                     # n
  print(Kode)                   # Kode
  print(rep(water.x, ly)[flux])      # x
  print("z")  
  print(z)                             # z
  print(as.vector(hNew)[flux])         # initial pressure head
  print(Q[flux])                       # discharge rate Q (water)
  print(horizon[flux])
  print("beta")
  print(beta[flux])                    # beta, water uptake distribution
  print(as.vector(h$millerH(rf[flux])))
  print(aK[flux])
  print(as.vector(h$millerT(rf[flux])))
  print(flux)
  print(c(length((1:NumNP)), length( Kode), length( rep(water.x, ly)[flux]),
          length(z), length(as.vector(hNew)[flux]), length(0), length(Q[flux]),
          length(horizon[flux]), length(beta[flux]),
          length(as.vector(h$millerH(rf[flux]))), length(aK[flux]),
          length(as.vector(h$millerT(rf[flux])))
          ))
  }
###################

  
  return(list(#Units=Units,
              Kat=Kat,
              MaxIt=MaxIt,
              TolTh=h$water$TolTh,
              TolH=h$water$TolH,
              lWat=h$water$lWat,     
              FreeD=FreeD,
              
              ## material
              NLay=NLay, 
              hTab=hTab, # c(hTab1, hTabN)
              Par=Par, # material information
              
              dt=h$water$dt,
              dtMinMax=dtMinMax,  # dtmin, dtmax
              DMul=DMul, # 
              TPrint=h$water$TPrint,
              
              ## root(sink infotmation), see after atmoshere
              
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
              lArtD=h$chem$lArtD,
              PeCr=PeCr,
              ChPar=ChPar,
              KodCB=KodCB,
              cBound=cBound, 
              tPulse=h$chem$tPulse,
              
              ## nodes
              nCodeM=cbind((1:NumNP),       # n
                Kode=Kode,                  # Kode
                x=rep(water.x, ly)[flux],   # x
                z=z,                        # z
                H=as.vector(hNew)[flux],    # initial pressure head
                C=0,                        # initial concentration Conc
                Q=Q[flux],                  # discharge rate Q (water)
                hor=horizon[flux], 
                b=beta[flux],               # beta, water uptake distribution
                aH=as.vector(h$millerH(rf[flux])),
                aK=aK[flux], 
                aT=as.vector(h$millerT(rf[flux]))
                ),
              
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
              atmosphere=atm.data,
              
              ## Block D table 8.4, logically after atmosphere
              ## (sink information)
              root=root, # 6 variables
              POptm=POptm[flux],
              
              ## additional information
              flux = flux,
              water.x = water.x,
              water.y = water.y
              )
         )
}
  
  
  
draw.horizons <- function(h, areas=TRUE,
                          col.hor=c('#000000', '#996600', '#660000', '#CC9933',
                            '#666600', '#CCCC99', '#CCCCCC', '#990000',
                            '#FFCC00', '#FFFFCC',rep('white',h$max.horizons)),
                          border.col=NULL, picture=NULL,
                          lwd = 2, quadratic=TRUE, all=TRUE,
                          cex=1, cex.leg=1
                          ) {
  stopifnot(length(col.hor) == h$max.horizons * 2)

  xlim <- range(h$grid.x) + c(-1,1) * h$step / 2 # 31.12.03: c(-1,1) * h$step / 2
  ylim <- range(h$grid.y) + c(-1,1) * h$step / 2
  
  if (quadratic) {
    m <- max(xd <- diff(xlim), yd <- diff(ylim))
    xlim[1] <- xlim[1] - (m - xd) / 2
    xlim[2] <- xlim[1] + m
    ylim[1] <- ylim[2] - m
  }

  if (all) 
    plot(Inf, Inf, xlim=xlim, ylim=ylim[2]-ylim, xlab="x", ylab="z",
         xaxs="i", yaxs="i", cex.axis=cex, cex.lab=cex)
  par(new=TRUE)
  plot(Inf, Inf, xlim=xlim, ylim=ylim, xlab="x", ylab="z",
         xaxs="i", yaxs="i", axes=FALSE)
  if (areas) {
    image(h$grid.x, h$grid.y, h$idx.rf %% h$max.horizons,
          col=col.hor, axes=FALSE,
          xlab="", ylab="", zlim=c(0, 2 * h$max.horizons - 1),
          xlim=xlim, ylim=ylim, xaxs="i", yaxs="i", add=TRUE)
    if (!is.null(border.col)) {
      border.col <- rep(border.col, length=h$max.horizons)
      for (i in 1:h$n)
        if (!is.null(h[[i]]$border)) {
          points((h[[i]]$border[,1] - 1) * h$step + h$grid.x[1],
                 (h[[i]]$border[,2] - 1) * h$step + h$grid.y[1],
                 pch=as.character(i), col=border.col[i], cex=0.4)
      }
    }
      
    legend(xlim[1], ylim[1], xj=0, yj=0, bg="white",
           pch=15, col=rep(col.hor, length=h$n), legend=paste(1:h$n),
           cex=cex.leg)
  } else {
    if (!is.null(picture) && all) {
      plotRGB(picture, x=h$grid.x, y=h$grid.y, 
                                    add=TRUE, xaxs="i", yaxs="i")
    }
    for (i in 1:h$n) 
      if (h[[i]]$type!="Start") {
        lines(h[[i]]$points$x, h[[i]]$points$y, col=col.hor[i], lwd=lwd)
      }
  }
  invisible(list(xlim=xlim, ylim=ylim))
}


calculate.horizons <- function(h) {
  linear <- function(x1, y1, x2, y2, x) 
    ((y2-y1) * x + y1 * x2 - y2 * x1) / (x2 - x1)
  xlim <- h$grid.x[c(1, length(h$grid.x))]
  if (h$n != as.integer(h$n)) stop("h$n is not an integer")
  if (is.null(h$max.horizons)) h$max.horizons <- RFparameters()$maxmodels
  h$n <- as.integer(h$n)
  stopifnot(h$n <= h$max.horizons)
  remove.history <- function(h) {
    for (i in 1:length(h)) {
      if (is.list(h[[i]]) && length(h[[i]])>0)
        h[[i]] <- remove.history(h[[i]])
    }
    h$.history <- NULL
    return(h)
  }
  h <- remove.history(h)
  h$hQThFlC <- NULL
  if (h$n>1) for (i in 2:h$n) if (h[[i]]$type=="H") {
    pts <- h[[i]]$points
    n.pts <- length(pts$x)
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
      pts$y <- c(linear(pts$x[1], pts$y[1], pts$x[2],
                        pts$y[2], xlim[1] - h$step / 2), pts$y)
      pts$x <- c(xlim[1] - h$step / 2, pts$x)
    }
    h[[i]]$points <- pts
  }
  .Call("GetHorizons", h, n=as.integer(c(1,h$n)), PACKAGE="SoPhy")
  h$plants <- h$plants.idx <- h$RF <- h$Stone.RF <- h$Root.RF <- NULL
  return(h)
}


modify.horizons <- function(h, percent=5, level.percent=5, rdistr=rnorm) {
  if (h$n==1) return(h)
  
  linear <- function(x1, y1, x2, y2, x) 
       ((y2-y1) * x + y1 * x2 - y2 * x1) / (x2 - x1)
  
  xlim <- range(h$grid.x)
  ylim <- range(h$grid.y)
  genuine.hor <- sum(unlist(lapply(h, function(x) x$type=="H"))) + 1
  for (i in 2:h$n) {
    s.x <- sqrt(var(h[[i]]$points$x)) * percent / 100
    s.y <- sqrt(var(h[[i]]$points$y)) * percent / 100
    l <- length(h[[i]]$points$x)
    if (h[[i]]$type=="H") {
      rho.y <- diff(range(h$grid.y)) /  genuine.hor * level.percent / 100
      h[[i]]$points$x <- h[[i]]$points$x + rdistr(l, s=s.x)
      if (h[[i]]$points$x[1] >= h[[i]]$points$x[2])
        h[[i]]$points$x[1] <- h[[i]]$points$x[2] - s.x + runif(1,max=s.x)
      h[[i]]$points$y <- h[[i]]$points$y + rdistr(1, s=rho.y) + rdistr(l, s=s.y)
      if (h[[i]]$points$x[l] <= h[[i]]$points$x[l-1])
        h[[i]]$points$x[l] <- h[[i]]$points$x[l-1] + s.x - runif(1,max=s.x)
    } else { ## type=="P"      
      rho.x <- diff(range(h$grid.x)) / (h$n - genuine.hor) * level.percent / 100
      rho.y <- diff(range(h$grid.y)) / (h$n - genuine.hor) * level.percent / 100
      h[[i]]$points$x <- h[[i]]$points$x + rdistr(1, s=rho.x) + rdistr(l, s=s.x)
      h[[i]]$points$y <- h[[i]]$points$y + rdistr(1, s=rho.y) + rdistr(l, s=s.y)
      h[[i]]$points$x[l] <- h[[i]]$points$x[1]
      h[[i]]$points$y[l] <- h[[i]]$points$y[1]
    }
  }
  return(calculate.horizons(h))
}

plotRF <- function(h, col.txt='black', col.stones='white',
                   col.rf=if (is.null(h$col.rf)) rainbow(100) else h$col.rf,
                   titl=TRUE, line=0.7,
                   lim=1, quadratic=TRUE, cex=1, cex.leg=0.8,
                   pch=5,
                   cex.pch = max(0.05, 3/min(length(h$grid.x),
                     length(h$grid.y))),
                   root.col=grey(1 - (i-1) / (1.6 * length(h$plants))),
                   what=c("Root.RF", "Stone.RF", "RF"),
                   transf=c("K", "H", "theta", "none"),
                   legend=TRUE, ylim, zlim) {
  root.x <- 1
  root.y <- 2

  what <- match.arg(what)
  transf <- match.arg(transf)
  HEX <- c(0:9, LETTERS[1:6])
  xlim <- range(h$grid.x)
  if (missing(ylim)) ylim <- range(h$grid.y)
  if (quadratic) {
    m <- max(xd <- diff(xlim), yd <- diff(ylim))
    xlim[1] <- xlim[1] - (m - xd) / 2
    xlim[2] <- xlim[1] + m
    ylim[1] <- ylim[2] - m
  }  
  
  txt <-
    paste(if (transf=="none") {if (what=="RF") "Gauss. " else
                                '"Gaussian field "' } else
      paste('alpha["', transf, ' "]', sep=""),
      switch(what, RF='("without stones & roots")',
             Stone.RF='("without roots")', Root.RF=''), sep="")

  rf <- switch(what, RF=h$RF, Stone.RF=h$Stone.RF, Root.RF=h$Root.RF)
  if (is.null(rf)) {return(paste("error:", txt, "has not been generated yet"))}
  if (transf!="none") rf <- h$miller.link(rf)
  rf <- switch(transf, K=h$millerK(rf), H=h$millerH(rf), theta= h$millerT(rf),
               none=rf)

  if (missing(zlim))
    zlim <- quantile(rf[is.finite(rf)], probs=c(0, lim), na.rm=TRUE)
  plot(Inf, Inf, xlim=xlim, ylim=ylim[2]-ylim, xlab="x", ylab="z",
       xaxs="i", yaxs="i", cex.axis=cex, cex.lab=cex)
  if (!is.logical(titl) || titl)
    title(if (is.logical(titl)) parse(text=txt) else titl,
          cex.main=cex * 1.2, col.main=col.txt, line=line)
  if (any(!is.finite(zlim))) return("error: too many non finite values")
   
  image(h$grid.x, h$grid.y, rf>zlim[2], axes=FALSE, add=TRUE,
        col=c(col.stones, "black"), xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
  par(new=TRUE) ## otherwise my.legend does not work!
  image(h$grid.x, h$grid.y, rf, ann=FALSE, axes=FALSE, zlim=zlim, # add=TRUE,
        col=col.rf, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
  if (what=="Root.RF" && length(h$plants)>0) {
    lc <- seq(100,200,l=length(h$plants))
    col <-  paste("#",
                  HEX[lc / 16], HEX[lc %% 16],
                  HEX[lc / 16], HEX[lc %% 16],
                  HEX[(lc-100) / 16], HEX[(lc-100) %% 16],sep="")
    for (i in 1:length(h$plants)) {
      ##if (!h$root[[h$plants.idx[i]]]$rf.Kf) ## ?! undecided
      {
        if (nrow(h$plants[[i]])>0) {
          points( t(t(h$plants[[i]][, c(root.x, root.y), drop=FALSE] - 1) *
                    h$step + c(h$grid.x[1], h$grid.x[1])),
                 # cex=max(0.178, 10/min(length(h$grid.x), length(h$grid.y))),
                 # cex=max(0.15, 8.5/min(length(h$grid.x), length(h$grid.y))),
                 cex=cex.pch, pch=pch, col=root.col)
        }
      }
    }
  }
  if (legend)
    my.legend(zlim=zlim, col=col.rf, lb.y=ylim[1], lb.x=xlim[1], cex=cex.leg)
  ## the following is necessary since otherwise the frame of the image plot
  ## is broken (due to error?! in R 1.8.0)
  par(new=TRUE)
  plot(Inf, Inf, xlim=xlim, ylim=ylim[2]-ylim, xlab="x", ylab="z",
       xaxs="i", yaxs="i", cex.axis=cex, cex.lab=cex)
  dimnames(rf) <- list(h$grid.x, h$grid.y)
  invisible(rf)
}


plotWater <- function(h, instance,
                      what=c("H", "Q", "theta", "vx", "vz", "Conc", "logH"),
                      col.txt="black",
                      col.simu=if (is.null(h$col.simu)) rainbow(100)
                      else if (what.nr==7) rev(h$col.simu) else h$col.simu,
                      col.exception = c("yellow", "darkred"), lim=1,
                      titl=TRUE, line=0.2, quadratic=TRUE, cex=1, cex.leg=0.8,
                      legend=TRUE, ylim, zlim){
  what.list <- eval(formals()$what)
  what <- if (!is.null(swp <- h$water$print) && length(what)==7) swp else what
  what.nr <- if (is.numeric(what)) what else pmatch(what, what.list)
  what <- what.list[what.nr]
  
  stopifnot(is.numeric(what.nr), length(what.nr)==1) 
  result <- matrix(nrow=length(h$water.x),ncol=length(h$water.y))
  if (missing(instance)) instance <- NULL
  
  result[h$flux] <-
    if (is.null(h$hQThFlC)) NA
    else h$hQThFlC[ 1 + ((what.nr - 1) %% 6), ,
                        if(is.null(instance)) dim(h$hQThFlC)[3]
                        else instance]
  if (!any(is.finite(result)))
    return("error: no finite value in the swms2d output")
  if (what.nr > 6) { # logH,what.nr=7
    idx <- result < 0  & !is.na(result)
    result[idx] <- log(-result[idx]) / log(10)
    result[!idx] <- -Inf
  }
  xlim <- range(h$water.x)
  if (missing(ylim)) ylim <- range(h$water.y)
  if (quadratic) {
    m <- max(xd <- diff(xlim), yd <- diff(ylim))
    xlim[1] <- xlim[1] - (m - xd) / 2
    xlim[2] <- xlim[1] + m
    ylim[1] <- ylim[2] - m
  }
  
  if (missing(zlim) || is.null(zlim)) {
    zlim <- result[is.finite(result)]
    zlim <- switch(what, 
                   H=c(min(0, quantile(zlim, probs=1-lim)), 0), 
                   Q=quantile(zlim, probs=c(1-lim, lim)), 
                   theta=range(0, unlist(lapply(h[1:10], 
                     function(x) x$materials$ths)), max(zlim)), 
                   vx=quantile(zlim, probs=c(1-lim, lim)), 
                   vz=quantile(zlim, probs=c(1-lim, lim)), 
                   Conc=c(0, max(0, quantile(zlim, probs=lim))),
                   logH=quantile(zlim, probs=c(1-lim, 1))
                   )
  }
  plot(Inf, Inf, xlim=xlim, ylim=ylim[2]-ylim, xlab="x", ylab="z",
       xaxs="i", yaxs="i", cex.axis=cex, cex.lab=cex)
  if (!is.logical(titl) || titl)
    title(if (is.logical(titl)) what else titl,
          cex.main=cex * 1.2, col.main=col.txt, line=line)
  par(new=TRUE)
  
  outside <- (result<zlim[1]) + 2 * (result>zlim[2])
  image(h$water.x, h$water.y, outside,
        ann=FALSE, axes=FALSE, col=c("white", col.exception),
        xlim=xlim, ylim=ylim, xaxs="i", yaxs="i", zlim=c(0,2))
  image(h$water.x, h$water.y, result, ann=FALSE, axes=FALSE,
        col=col.simu, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i", zlim=zlim,
        add=TRUE)
  if (legend)
    my.legend(zlim=zlim, col=col.simu, lb.y=ylim[1], lb.x=xlim[1], cex=cex.leg,
              y.i=2/length(col.simu))
  dimnames(result) <- list(h$water.x, h$water.y)
  invisible(result)
}


my.legend <- function(lb.x=0, lb.y=0, zlim, col, cex=0.8, y.intersp=0.02,
                      bg='white') {
  ## uses already the legend code of R-1.3.0
  cn <- length(col)
  filler <- vector("character", length=(cn-3)/2)
  bg.orig <- par()$bg
  par(bg=bg)
  m <- mean(zlim)
  if (abs(m) < zlim[2] / 1000) m <- 0
  leg <- c(zlim[1], m, zlim[2])

  ## neu
  leg <- round(leg, floor(2 - sort(log(abs(leg))/log(10))[2]))
  
  legend(lb.x, lb.y, y.i=y.intersp, x.i=0.1, yj=0,
         legend=   c(formatC(leg[3], dig=2), filler,
           formatC(leg[2], dig=2), filler,
           formatC(leg[1], dig=2)),
          lty=1, lwd=0.1, col=rev(col), cex=cex)
  par(bg=bg.orig)
  invisible(NULL)
}










tortuosity <- function(depth, freq, range.depth, range.xi, len.xi=100,
                       lower.bound.m = 1.005, tolerance.m=0.01, PrintLevel=0) {
  
    H.asterix <- function(x, xi, s=1, lower.tail=TRUE, a, b) {
      stopifnot(xi!=0)
      res <- (pmax(0, 1 + xi * a * x / s)^(1 - 1/xi) -
              pmax(0, 1 + xi * b * x / s)^(1 - 1/xi)) /
                ((1-xi) * x / s * (b - a))
      res[x==0] <- 1
      if (lower.tail) 1 - res else res
    }

    target <- function(x) {
      sum(abs(mx.fr * H.asterix(dpth, xi, s=x[1], FALSE, 1, b=x[2]) -
              frq)^2)
    }

    lower.bound.b <- 2 * lower.bound.m - 1
    max.frq <-  max(freq, na.rm=TRUE)
    stopifnot(all(range.xi < 0), tolerance.m >= 0)
    xi.s <- seq(min(range.xi), max(range.xi), len=len.xi)
    d.s <- which(depth >= min(range.depth) & depth <= max(range.depth))    
    s <- x <- b <- v <- matrix(ncol=length(d.s), nrow=length(xi.s))
    Dpth0 <- numeric(length(d.s))
    for (ni in 1:length(d.s)) {
      i <- d.s[ni]
      if (PrintLevel>1) cat(i, "")
      frq <- freq[-1:-i]
      Dpth0[ni] <- depth[i+1]
      dpth <- depth[-1:-i] - Dpth0[ni]
      idx <- is.finite(frq)
      frq <- frq[idx]
      dpth <- dpth[idx]
      mx.fr <- max(frq)
      for (nxi in 1:length(xi.s)) {
        xi <- xi.s[nxi]
        init <- c(40, 2)
   #     print(target(init)) 
        opt <- optim(init, # 1. scale; 2. b
                     target, lower=c(0.1, lower.bound.b),
                     upper=c(1000, 100), method="L-BFGS-B",
                     control=list(parscale=c(50,2),
                       fnscale=target(init) / 10)
                     )
  #      print(opt$val / target(init) * 10)
        if (opt$val<Inf && all(is.finite(opt$par))) {
          s[nxi, ni] <- opt$par[1]
          b[nxi, ni] <- opt$par[2]
          x[nxi, ni] <- xi
          v[nxi, ni] <- opt$val # currently not used
        }
      }
    }

    threshold <- !apply(is.na(x), 1, any)
    X <- x[threshold, 1]
    B <- diff(apply(b, 1, range, na.rm=TRUE))[threshold]
    min.i <- B == min(B[apply(b > lower.bound.b, 1, all)])
    Xmin <- X[min.i]
    diffB0 <- B[min.i]
    if (PrintLevel>2) print(min.i)
    min.i <- B <= diffB0 * (1 + tolerance.m) & apply(b > lower.bound.b, 1, all)
    if (PrintLevel>2) print(min.i)
    
    b.sel <- as.vector(b[min.i, ])
    s.sel <- as.vector(s[min.i, ])
    opt.idx <- order(b.sel)[round(length(b.sel) / 2 + 0.5 +
                                  (tolerance.m==0) * c(-0.1, 0.1))]
    opt.b <- mean(b.sel[opt.idx])
    opt.s <- mean(s.sel[opt.idx])
    opt.dpth0 <- mean(rep(Dpth0, each=sum(min.i))[opt.idx])
    dpth <- depth[depth >= min(range.depth)]

    max.frq <- max(freq[depth >=  opt.dpth0], na.rm=TRUE)
 
    return(list(xi=X, raw.m=(1+b)/2, span.m=B/2,
                opt=list(xi=Xmin, m = (1+opt.b)/2, s = opt.s,
                  mspan=diffB0 / 2, D=opt.dpth0), 
                input=list(freq=freq, depth=depth), 
                fitted=list(depth=dpth,
                  p=if (tolerance.m ==0) max.frq * H.asterix(dpth-opt.dpth0,
                          Xmin, s=opt.s, FALSE, 1, opt.b)
                  else rep(NA, length(dpth))
                )
              )
         )
}


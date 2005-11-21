## this function behaves quite different from tortuosity
## not investigated way; but guess, for numerical reasons.

tortuosity2 <- function( depth, freq, range.depth, range.xi, len.xi=100,
                       lower.bound.m = 1.005, PrintLevel=0) {
  
    H.asterix <- function(x, xi, s=1, lower.tail=TRUE, m) {
      stopifnot(xi!=0)
      res <- (pmax(0, 1 + xi * x / s)^(1 - 1/xi) -
              pmax(0, 1 + xi * (2 * m - 1) * x / s)^(1 - 1/xi)) /
                ((1-xi) * x / s * 2 * (m - 1))
      res[x==0] <- 1
      if (lower.tail) 1 - res else res
    }

    target <- function(x) {
      sum(abs(mx.fr * H.asterix(dpth, xi, s=x[1], FALSE, m=x[2]) -
              frq)^2)
    }

    max.frq <-  max(freq, na.rm=TRUE)
    xi.s <- seq(min(range.xi), max(range.xi), len=len.xi)
    d.s <- which(depth >= min(range.depth) & depth <= max(range.depth))    
    s <- x <- m <- v <- matrix(ncol=length(d.s), nrow=length(xi.s))
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
        init <- c(40, 1.5)
        print(target(init))
        opt <- optim(c(40, 2), # 1. scale; 2. m
                     target, lower=c(0.1, lower.bound.m),
                     upper=c(1000, 50.5), meth="L-BFGS-B",
                     control=list(parscale=c(50, 1.5),
                       fnscale=target(init) / 10))
        print(opt$val / target(init) * 10)
        if (opt$val<Inf && all(is.finite(opt$par))) {
          s[nxi, ni] <- opt$par[1]
          m[nxi, ni] <- opt$par[2]
          x[nxi, ni] <- xi
          v[nxi, ni] <- opt$val # currently not used
        }
      }
    }

    threshold <- !apply(is.na(x), 1, any)
    X <- x[threshold, 1]
    M <- diff(apply(m, 1, range, na.rm=TRUE))[threshold]
    min.i <- which(M == min(M[apply(m > lower.bound.m, 1, all)]))
    print(min.i)
    Xmin <- X[min.i]
    diffM0 <- M[min.i]
    
    opt.idx <- order(m[min.i, ])[round(ncol(m) / 2 + +0.5 + c(-0.1, 0.1))]
    opt.m <- mean(m[min.i, opt.idx])
    opt.s <- mean(s[min.i, opt.idx])
    opt.dpth0 <- mean(Dpth0[opt.idx])
    dpth <- depth[which(depth >= min(range.depth))]

    return(list(raw.x=x, raw.m=m, span.x=X, span.m=M,
                opt=list(x=Xmin, m = opt.m, s = opt.s,
                  mspan=diffM0), # da f linear!!!
                input=list(freq=freq, depth=depth), 
                fitted=list(depth=dpth, p=max.frq * H.asterix(dpth-opt.dpth0,
                                        Xmin, s=opt.s, FALSE, opt.m)
                )
              )
         )

}


if (FALSE) {

   mtort2 <- tortuosity2(truedepth, freq, 
                    range.depth=c(85, 95-7), range.xi=c(-2, -0.6),
                    Print=2
                    )
    str(mtort2)
    
    mtort <- mtort1 <-tortuosity(truedepth, freq, 
                    range.depth=c(85, 95-7), range.xi=c(-2, -0.6),
                    Print=2
                    )
    str(mtort1)
    
    for (n in names(mtort1))
      print(c(n, sum(abs(unlist(mtort1[n][[1]]) -
                         unlist(mtort2[n][[1]])), na.rm=TRUE)))

    print(mtort1$raw.m)
    print(mtort2$raw.m)
    print(mtort2$raw.m-mtort1$raw.m)

 }

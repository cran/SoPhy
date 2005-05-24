# nice -20 R --no-save < runJCH.R > out
# tail -f out
# source("runJCH.R")

if (EXTENDED.TESTING <- file.exists("source.R")) source("source.R")

if (EXTENDED.TESTING) {

dev <- TRUE
excep.dev <- 2
final <- TRUE
final <- FALSE

path <- "sh.dir/"
if (!file.exists(path)) dir.create(path)
ps <- paste(path, "eps", sep="")
txt <- paste(path, "txt", sep="")
simu <- paste(path, "simu", sep="")

name.jch <- paste(path, "jconthydr.data.tex", sep="")
if (file.exists(name.jch)) load(name.jch) else {
  runif(1)
  seed.jch.simu5 <- seed.jch.start <- seed.jch.sketches <- seed.jch.simu1 <-
    seed.jch.simu3 <- seed.jch.simu2 <- seed.jch.simu4 <- .Random.seed
}

Pr <- 2 + (!is.logical(dev)) 

## runif(1); seed.jch.simu5 <- .Random.seed;
cat("\nrunJCH.R data\n")
.Random.seed <- seed.jch.simu5
sh.jch(input=c(9, 10,
         11, 0), dev=dev, Pr=Pr, final=final, ps=ps, txt=txt, simu=simu)

## runif(1); seed.jch.start <- .Random.seed;
cat("\nrunJCH.R start\n")
.Random.seed <- seed.jch.start
sh.jch(input=0, final=final, Pr=Pr, ps=ps, txt=txt, simu=simu)


## runif(1); seed.jch.sketches <- .Random.seed;
cat("\nrunJCH.R sketches\n")
.Random.seed <- seed.jch.sketches
sh.jch(input=c(1, 0), dev=dev, Pr=Pr, final=final, ps=ps, txt=txt, simu=simu)

## runif(1); seed.jch.simu1 <- .Random.seed;
cat("\nrunJCH.R simu1\n")
.Random.seed <- seed.jch.simu1
sh.jch(input=c(14, 2, 0), dev=dev, Pr=Pr, final=final, ps=ps, txt=txt, simu=simu)


## runif(1); seed.jch.simu3 <- .Random.seed;
cat("\nrunJCH.R simu3\n")
.Random.seed <- seed.jch.simu3
sh.jch(input=c(12, 13,
         0), dev=excep.dev, Pr=Pr, final=final, ps=ps, txt=txt, simu=simu)

## runif(1); seed.jch.simu2 <- .Random.seed;
cat("\nrunJCH.R simu2\n")
.Random.seed <- seed.jch.simu2
sh.jch(input=c(3, 4, 0), dev=dev, Pr=Pr, final=final, ps=ps, txt=txt, simu=simu)

## runif(1); seed.jch.simu4 <- .Random.seed;
cat("\nrunJCH.R simu4\n")
.Random.seed <- seed.jch.simu4
sh.jch(input=c(6,
         8, 7, 0), dev=dev, Pr=Pr, final=final, ps=ps, txt=txt, simu=simu)

save(file=name.jch,
     seed.jch.simu1, seed.jch.simu2, seed.jch.simu3, seed.jch.simu4,
     seed.jch.simu5,
     seed.jch.sketches, seed.jch.start)
}

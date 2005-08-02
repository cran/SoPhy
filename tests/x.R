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

## runif(1); seed.jch.sketches <- .Random.seed;
cat("\nrunJCH.R sketches\n")
.Random.seed <- seed.jch.sketches
sh.jch(input=c(1, 2, 5, 0), dev=dev, Pr=Pr, final=final, ps=ps, txt=txt,
       simu=simu)
}

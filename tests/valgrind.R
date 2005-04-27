## R -d "valgrind --tool=memcheck --leak-check=yes" < valgrind.R 
## R -d "valgrind --tool=memcheck --leak-check=yes"
## source("valgrind.R")

if (FALSE)
  {

#detach(package:SoPhy)
library(SoPhy)
source("/home/schlather/article/R/SOPHY/SoPhy/R/swms2d.R")
if (TRUE) {
  runif(1)
  seed <- .Random.seed
  save(file="seed", seed)
} else {
  load("seed")
   .Random.seed <- seed
}

h <- xswms2d(xlim=c(1, 100), ylim=c(1, 100), step=1, new=NULL)
h$H2$model$model[[1]]$model <- "nugget" ## * * * * * * * * * * * *
h <- simulate.horizons(h)        ## simulate stochastic components
swms2d.out <- swms2d(h, iter.print=1) ## numerical simulation


example(calculate.horizons) ## nicht ok, dauert
example(CDE)
# example(analyse.profile)  ## dauert viel zu lange
example(create.roots)
example(create.stones)
example(create.waterflow)
example(draw.horizons) # ok
example(flowpattern) #ok
example(modify.horizons)
example(my.legend) #ok
example(plotFlow) #ok
example(plotRF)#ok
example(quader)#ok
example(read.picture) #ok
example(read.swms2d.table) # ok
example(risk.index) #ok
example(sh.jh) # ok
example(simulate.horizons) #ok
example(swms2d)  #  *** not ok!!
example(tracer) #ok
example(xswms2d)#ok


}

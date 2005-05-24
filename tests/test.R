# source("test.R")

if (EXTENDED.TESTING <- file.exists("source.R")) source("source.R")

## redo example 1 of SWMS_2d

path <- paste(system.file(package='SoPhy'), 'swms2d', sep="/")
x <- read.swms2d.table(path)
res <- swms2d(x,
              iter.print = 1,
              max.iteration = 2, message=function(t) TRUE,
              breakpoint=1, intermediate.result=function(...) {cat("#")},
              )
print(paste("cat ", path, "/BALANCE.OUT", sep=""))
system(paste("cat ", path, "/BALANCE.OUT", sep=""))
print(res$balance)


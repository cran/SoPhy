
library(SoPhy, lib=if (file.exists("~/TMP/SoPhy")) "~/TMP")

if (EXTENDED.TESTING <- FALSE) {
s <- "/home/schlather/R/old.randomfields/RandomFields/tests/source.R"
if (EXTENDED.TESTING <- file.exists(s)) source(s) else {
  s <- "/home/schlather/R/RF/RandomFields/tests/source.R"
  if (EXTENDED.TESTING <- file.exists(s)) source(s)
  ## EXTENDED.TESTING <- FALSE
}
}

.path <- "~/R/SOPHY/SoPhy/R/"
if (file.exists(paste(.path, "swms2d.R", sep=""))) {
  Source <- function(x) {
    cat(x, "...\n")
    x <- paste(.path, x, sep="")
    source(x)
  }
  Source("3dplot.R")
  Source("Sophy.R")
  Source("analytic.R")
  Source("extremal.coeff.R")
  Source("rgb.R")
  Source("sh.R")
  Source("simu.R")
  Source("swms2d.R")
  Source("tortuosity.R")
}

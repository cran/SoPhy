
library(SoPhy, lib=if (file.exists("~/TMP/SoPhy")) "~/TMP")

s <- "/home/schlather/article/R/NEW.RF/RandomFields/tests/source.R"
if (file.exists(s)) source(s)

.path <- "/home/schlather/article/R/SOPHY/SoPhy/R/"
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
  Source("simu.R")
  Source("swms2d.R")
}

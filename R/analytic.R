ADE <- function(z, time, C0, dispersion, velocity) {
  nt <- as.integer(length(time))
  nz <- as.integer(length(z))
  conc <- .Fortran('ade', as.double(C0), as.double(dispersion),
                   as.double(velocity), as.double(time), nt,
                   as.double(z), nz, conc=double(nt * nz),
                   PACKAGE="SoPhy")$conc
  return(matrix(conc, nrow=nz))
}


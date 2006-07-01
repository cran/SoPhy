
read.picture <- function(picture, extensions=c("tif", "tiff",
                               if (.Platform$OS.type=="unix") c("gif", "jpeg")),
                         PrintLevel=RFparameters()$Print, tmp.dir=".") {
#  if (.Platform$OS.type!="unix") {
#    warning("`read.picture' currently only available on unix systems")
#    return(NULL)
#  }
  ## search for file -- either the file should have one of the `extensions'
  ## or the file name expanded by one of the `extension', in the given sequence
  ##
  ## if it is not tif it converts to tif first
  extensions <- unique(tolower(extensions))
  path <- dirname(picture)
  pic <- basename(picture)
  if (is.null(tmp.dir)) tmp.dir <- path
  if (substr(tmp.dir, nchar(tmp.dir), nchar(tmp.dir))!="/")
    tmp.dir <- paste(tmp.dir, "/", sep="")
  pic <- strsplit(pic,"\\.")[[1]]
  if (length(pic)>1) {
    ext <- tolower(ext.save <- pic[length(pic)])
    pic <- paste(pic[-length(pic)], collapse=".")
  } else {
    ext.save <- NULL
    ext <- ""
  }
  s <- paste(path, pic, sep="/")
  sfound <- picture
  if (length(i <- which(ext==extensions)) == 0) {
    ## was ".ext" not an extension, but part of the name ?
    if (!is.null(ext.save)) {
      s <- paste(s, ext.save, sep=".")
      pic <- paste(pic, ext.save, sep=".")
    }
    i <- 1
    while (i<=length(extensions)) {
      if (PrintLevel>2)
        cat(paste(s, extensions[i], sep="."),"\n")
      if (file.exists(paste(s, exten <- extensions[i], sep=".")))
        break
      if (file.exists(paste(s, exten <- toupper(extensions[i]), sep=".")))
        break
      i <- i + 1
      }
    if (i>length(extensions)) {
      stop("file ", picture, " not found")
    } else {
      sfound <- paste(s, exten, sep=".")
    }
    if (PrintLevel>1) cat("found", sfound, "\n")
  }
  if (is.na(pmatch(extensions[i], c("tif", "tiff")))) {## not tiff
    if (.Platform$OS.type!="unix")
      stop("cannot handle", extensions[i], "files within ", .Platform$OS.type)
    s <- paste(tmp.dir, pic, ".tif", sep="")
    if (file.exists(s)) warning(paste(s, "used in instead of", sfound))
    else {
      if (PrintLevel>2) cat("converting", sfound, "to", s, "\n")
      system(paste("convert ", sfound, s))
      if (!file.exists(s)) stop("conversion to tif did not work out.")
    }
  } else s <- sfound
  .Call("readtiff", s, PACKAGE="SoPhy")
}

write.picture <- function(picture, tif) {
#  if (.Platform$OS.type!="unix")
#    stop("`write.picture' currently only available on unix systems")
  stopifnot(is.array(picture), is.numeric(picture))
  storage.mode(picture) <- "integer"
  if (!(rev(strsplit(tif, "\\.")[[1]])[1] %in% c("TIFF", "TIF", "tiff", "tif")))
    tif <- paste(tif, ".tif", sep="")
#  print(tif)
  ppm <- paste(tif, "...ppm", sep="")
  while (file.exists(ppm)) ppm <- paste(ppm, "...ppm", sep="")
  .Call("writetiff", picture, tif, PACKAGE="SoPhy")
  txt <- paste("rawtoppm", dim(picture)[1], dim(picture)[2],
               tif, ">", ppm, "; convert -flip ", ppm, tif, "; rm", ppm)
  invisible(system(txt))
}

plotRGB <- function(picture, x=1:dp[1], y=if (reverse) dp[2]:1 else 1:dp[2],
                     reverse=TRUE, cex.axis=1, ...) {  
  max.pixels <- 2000000
  dp <- dim(picture)
  stopifnot(length(dp)==3, dp[3] %in% 3:4)
  
  if (dp[1] * dp[2] > max.pixels) {
    f <- ceiling(sqrt(dp[1] * dp[2] / max.pixels)) - 1
    xi <- rep(c(TRUE, rep(FALSE, f)), len=length(x))
    yi <- rep(c(TRUE, rep(FALSE, f)), len=length(y))
    picture <- picture[xi, yi, 1:3]
    x <- x[xi]
    y <- y[yi]
    dp <- dim(picture)
  }
  HEX <- c(0:9, LETTERS[1:6]) 
  HEX <- paste(HEX[1 + (0:255) / 16], HEX[1 + (0:255) %% 16], sep="")
  m <- matrix(1:prod(dp[-3]), nrow=dp[1], ncol=dp[2])
  picture[is.na(picture)] <- 255
  col <- paste("#",HEX[picture[,,1] + 1], HEX[picture[,,2] + 1],
                  HEX[picture[,,3] + 1], sep="")
  if (y1Gy2 <- y[1]>y[2]) y <- -y
  image(x, y, m, col=col, axes=!y1Gy2, frame.plot=TRUE, cex.axis=cex.axis, ...)
  if (y1Gy2) {
    axis(1, cex.axis=cex.axis)
    pry <- pretty(y)
    axis(2, at=pry, labels=-pry, cex.axis=cex.axis)
  }
  par(new=TRUE)
  plot(Inf, Inf, axes=FALSE, frame.plot=TRUE, xlim=0:1, ylim=0:1,
       xlab="", ylab="")
  invisible()
}

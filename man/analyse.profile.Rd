\name{analyse.profile}
\alias{analyse.profile}
\title{Analysis of a flow pattern in soil}
\description{
  The function fits the Pareto distribution to
  a flow pattern obtained from a dye experiment
}
\usage{
 analyse.profile(picture, fct, param, lower, upper, loc,
           estimate.all=NULL, selected.dist=0.95, selected.rate=c(0.5, 0.8),
           measure=function(x) x^2, method=c("fix.m","optim.m", "ml"),
	   endpoint.tolerance=-10,
           ppi.par = 5,
           ppi.xy = c(xlow=2, xup=2, ylow=2, yup=2),

           interactive=TRUE, PrintLevel=RFparameters()$Print,
           extensions=c("tif", "tiff", "gif", "jpeg"),

           X11.width=13, new=TRUE,
           col.thresh=c("white", "yellow", "black", "blue"),
           col.rect="red", col.bg="yellow", col.sep="gray",
           col.left="red", col.mid="darkgreen", col.right="darkgreen",
           col.line = "red", col.txt="black",
           reverse=TRUE, title
          )
}

\arguments{
  \item{picture}{list or array or character.
    If is a list, then the programme expects that the structure is
    of the output format
    of \command{analyse.profile}.  Then the parameters
    \code{fct}, \code{param},
    \code{upper}, \code{lower} and \code{loc} must be missing, and they
    are taken from the list.
    The list element \code{picture} might be missing and is reread from
    the orginal file.
    \cr

    If it is an array, a three dimensional array is expected where the third
    dimension has three or four components (the forth component is
    ignored if any). The array is intepreted as RGB-A coded picture.\cr

    If \code{picture} is a character string, the function reads in a
    file with name or name base 
    \code{picture}, see \command{\link{read.picture}}.
  }
  \item{fct}{threshold function that extracts the stained area from an
    RGB figure. \code{fct} has two input variables, the first is a three
    dimensional array, the second is an arbitrary list of \emph{named}
    parameters. Further, the function must return a list of standard
    values of the parameters, if the second parameter is missing,
    independently of value of the first variable.

    If \code{fct} is not given, then

    \preformatted{fct <- function(i, p) {
        if (missing(p)) list(minR=160, maxGB=200)
        else {
	  gb <- i[,,2] + i[,,3]
          (gb >= i[,,1] * p$minR / 100 & 80<=gb &
	   p$maxGB>=gb)
        }
      }
    }
  }
  \item{param}{list.  Initial values for the parameter list
    of the function \code{fct}.
    If  \code{param} is not named the names will be \sQuote{parameter1},
    \sQuote{parameter2}, etc.  Be sure that the parameter names match the
    names in \code{fct}.
    \cr
    
    The values of \code{param} are changed in the interactive
    plot.  If \code{param} is not given \code{param} takes the value
    returned by \code{fct()}.
  }
  \item{lower}{inital values for the lower bounds of the
    range of the parameters in \code{fct}.  The values can be changed in
    the interactive plot.  If not given \code{lower} is set to \code{param}.}
  \item{upper}{inital values for the upper bounds of the
    range of the parameters in \code{fct}.  The values can be changed in
    the interactive plot. 
    If not given \code{upper} is set to \code{param}.}
  \item{loc}{missing or list or list of lists.
    If it is a list then it must have the
    the components \code{x} and \code{y}.  Each of the
    components are vectors of two components containing the x and y
    coordinates of two opposite corners of a rectangle, respectively.
    The rectangle gives the clipping area for which the number of
    stained pixels is calculated.  The values can be changed interactively.
    \cr
    
    If \code{loc} is missing \code{loc} is determined in the interactive plot.
    \cr
    
    If \code{loc} is a
    list of lists then it has three components, each of
    which is a list as described above.
    The first component gives the clipping reactangle (that might be
    thought to be best), the other two the
    minimal and maximal tolerable cutting areas.
    If an edge of the last two rectangle agrees with an edge of the
    best rectangle up to 2 pixels, it is assumed that the user
    wishes that the two edges agree.

    Note that the ideal rectangle is first given by two points which
    are the opposite corners of the rectangle.
    This rectangle and the following ones are changed by giving a further
    point; the programme guesses suitably well which corner should
    be moved into which position.
    \cr
    
    if \code{loc} is not given and \code{interactive=FALSE} then
    the whole picture is taken.
  }
  \item{estimate.all}{\code{NULL} or logical. If \code{NULL} then
    nothing is estimated. If \code{FALSE} only a single risk index
    according to \code{param} is
    estimated; if \code{TRUE} the risk index is calculated for a grid
    given by \code{lower}, \code{upper} and the
    second and third component of \code{loc}.
    The marginal number of grid points is \code{ppi.par} for the
    components of \code{lower} and \code{upper} and \code{ppi.xy}
    for the definition of the clipping rectangle;
    see also \command{\link{risk.index}}}
  \item{selected.dist}{number in \eqn{(0,1)} or a vector of integers. 
    Distances for which the form parameter of the pareto distribution
    is estimated; see \command{\link{risk.index}}. 
    If \code{selected.dist} is a number in \eqn{(0,1)} the distances are
    \code{1:(1 + round(max(data[,1]-1) * selected.dist))}.
 }
  \item{selected.rate}{additionally to the indices given by
    \code{selected.dist}, the form parameter is estimated also for
    those distances where the corresponding number of observed relative stained
    pixels (w.r.t. to the maximum number of observed pixels)
    is within the interval given by \code{selected.rate}. 
    The risk index is calculated as the median of the estimated form parameters.}
  \item{measure}{instead of the default least squares another distance
    function can be given.}
  \item{method}{the number of observed paths is a free parameter when
    fitting the Pareto distribution. It can either be set as the maximum
    number of stained pixels for the currently considered distances or
    depths (\code{'fix.m'}) or fitted within the optimisation algorithm
    (\code{'optim.m'}). Usually, it is not worth using the slower
    \code{'optim.m'} option. As third new option the maximum likelihood
    method "ml" is available.}
  \item{endpoint.tolerance}{optimisation parameter.  If the shape
    parameter is negative then the distribution has a finite upper
    endpoint.  Hence, mathematically, the lowest upper end point of the
    Pareto distribution is given as the largest distance for which at
    least one stained pixel is observed. 
    For stability reasons and because the observed data might be a
    scale mixure of Pareto distribution it is advantageous to allow for
    some tolerance of the minimal upper end point.
    
    If \code{endpoint.tolerance} is positive then the lower threshold
    for the upper end point is the largest distance for which the
    number of observed stained pixels is larger than
    \code{endpoint.tolerance}.

    If \code{endpoint.tolerance} is negative then the lower threshold
    equals the largest distance for which at least one stained pixel is
    observed minus the modulus of \code{endpoint.tolerance}.
  }
  \item{ppi.par}{list or scale. Number of points by which
    a given parameter interval is divided
    (if the interval is not trivial). The points of all intervals (and those for
    the clipping area) build a grid. For each point of the grid, the
    risk index is estimated, if \code{estimate.all=TRUE}. \code{ppi.par}
    must be given by name (or a single integer is given then the number
    of grid points are the same for each parameter);
    see also the variable \code{param}.

    If \code{ppi.par} is \code{NULL} and \code{picture} is a list,
    the value in \code{picture} is used.
  }
  \item{ppi.xy}{vector of 4 components; number of points by which the interval
    for the positions of the lower, upper, left and right edge of the
    clipping area is divided.

    If \code{ppi.xy} is \code{NULL} and \code{picture} is a list,
    the value in \code{picture} is used.
  }
  \item{interactive}{logical.  If \code{TRUE} the user can
    modify the
    parameters interacitvely.}
  \item{PrintLevel}{numeric.
    \itemize{
      \item{\eqn{<}=0}{no messages are printed}
      \item{\eqn{>}=1}{the higher the number the more information is given}
    }
  }
  \item{extensions}{accepted extensions for the graphics.  Used if
    \code{picture} is a character}
  \item{X11.width}{width of the interactive plot. The height is
    determined automatically.  Only used if \code{new=TRUE}}
  \item{new}{if \code{TRUE} a new window for the interactive plot is opened}
  \item{col.thresh}{vector of 4 colour specifications. The components
    are the colours for
    \enumerate{
      \item the pixels that are  not considered as stained
      by the current and the preceding parameter choice
      \item  the pixels that are considered as stained by the preceding
      parameter choice, but not by the current one.
      \item  the pixels that are considered as stained by the current parameter
      choice, but not by the preceding one
      \item  the pixels that are considered as stained by the current
      and the preceding parameter choice
    }
    If the values of an interval for a parameter
    is changed than the parameters in \code{param}
    are considered as the the
    preceding parameters.
  }
  \item{col.rect}{colour for the interactive menue: colour of a button for free input.}
  \item{col.bg}{colour for the interactive menue: colour of an interactive bar}
  \item{col.sep}{colour for the interactive menue: colour of a separating line}
  \item{col.left}{colour for the interactive menue: colour of a preceding element}
  \item{col.mid}{colour for the interactive menue: colour for a message}
  \item{col.right}{colour for the interactive menue: colour of a subsequent element}
  \item{col.line}{colour for the interactive menue: colour of a marking
    line in interactive bars of direct choice.}
  \item{col.txt}{colour for the interactive menue: colour of headers}
  \item{reverse}{logical; see details}
  \item{title}{title for the interactive plot; if not given the picture
    name is used (usually the file name)}
}
\details{
  If the picture is read in from the file the following particular case
  may appear.  In soil profiles often a metal grid is placed of size 1m
  x 1m.
  This metal grid must be first removed by usingGIMP (www.gimp.org),
  for example, and,
  in general, also the figure has to be rotated.
  The pixels of the metal bars are removed by coding them 255, 255, 255
  in RGB code.  To get a reliable estimation of the risk index, the
  number of removed pixels must be the same  for each horizontal line
  in the chosen clipping area, except the horizontal line consists of
  white pixels only that should be removed.
  Now the programme sets each pixel of a horizontal line to NA
  if at least 70\% of
  the pixels of that line are coded 255, 255, 255.  Otherwise only the
  white pixels are set to NA.
  \cr
  
  Plotting aspects:
  
  \code{reverse=TRUE} is for the plot of soil profiles, where the
  numbering of the y-axis start at the top and is counted positively
  downwards.  

  \code{reverse=FALSE}
  The usual presentation of the vertical axis.
 
  The functions checks whether \code{y} is in decreasing or increasing
  order, and plots the y axis accordingly. 
}
\value{
  returns a list of the following arguments. Since the arguments do not
  have a fixed order they should be called by name, only.
  \item{picture}{a three dimensional array, containing the RGB-A coding}
  \item{fct}{see the variable \code{fct} above}
  \item{param}{see the variable \code{param} above}
  \item{lower}{see the variable \code{lower} above}
  \item{upper}{see the variable \code{upper} above}
  \item{loc}{see the variable \code{loc} above}
  \item{ppi.par}{see the variable \code{ppi.par} above}
  \item{ppi.xy}{see the variable \code{ppi.xy} above}
  \item{r.i}{output of \command{\link{risk.index}} for the ideal parameters}
  \item{risk.index}{array of \code{2 + length(param)} dimensions
    containing the estimated risk indeces obtained by
    the sensitivity analysis. The latter is
    based on intervals of values for the edges of the clipping area, see
    \code{loc} and \code{ppi.xy} and the intervals for the parameters, see
    \code{lower}, \code{upper} and \code{ppi.par}.  The first two
    dimensions of \code{risk.index}
    are  for the horizontal and vertical variation of the
    rectangle, the following dimensions are for the parameters, given in
    the order of \code{param}. The number of components in each
    dimension is given by \code{ppi.xy} and \code{ppi.par} except an
    interval is trivial.  Then the number of components is 1.
  }
  \item{raw.risk.index}{
    same as \code{risk.index} except that the raw risk index is stored,
    see the function \command{\link{risk.index}}.

    The raw risk index is usually not used.
  }
  \item{absfreq}{vector. Number of blue pixels per line (from top)
  }
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/~schlather}}

\references{
  Schlather, M. and Huwe, B. (2005)
  A risk index for characterising flow pattern in soils
  using dye tracer distributions
  \emph{Submitted to J. Contam. Hydrol.} % \bold{}, .
}


\seealso{%\command{\link{define.horizons}},
  \command{\link{risk.index}}, \code{\link{SoPhy}}
}


\examples{
%% library(SoPhy, "~/TMP")
%% load('~/soil.contamination/R/data/K06.G.tif.l')
%% K06 <- l; 
%% plotRGB(K06$picture[100:500, 200:400, 1:3])
%% K06$picture <- K06$picture[100:500, 200:400 , 1:3];
%% K06$loc <- rep(list(list(x=c(50,350), y=c(50,150))), 3)
%% environment(K06$fct) <- NULL
%% save(file=paste('~/R/SOPHY/SoPhy/data/K06.RData'), K06, compress=TRUE)

% library(SoPhy, lib="~/TMP")
% source("/home/schlather/R/SOPHY/SoPhy/R/extremal.coeff.R")
% close.screen(close.screen()); dev.off(); dev.off(); dev.off();

data("K06")
path <- paste(system.file(package='SoPhy'),"tracer", sep="/")
K06$name <- paste(path, "K06.G.tif", sep="/")
rate <- if (interactive()) c(0.5, 0.8) else c(0.6, 0.67)
ri <- analyse.profile(K06, estimate.all=FALSE, method="fix.m",
                      selected.rate=rate, selected.dist=NULL,
                      Print=2)
cat("estimated risk index =", ri$r.i$risk.index, "\n")
}

\keyword{spatial}

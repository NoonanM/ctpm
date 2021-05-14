\name{plot.variogram}
\alias{plot.variogram}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting methods for variogram objects.}
\description{
Produces plots of varigram objects (semi-variance vs. time lag) and model semi-variance functions, with approximate confidence intervals around the semi-variance estimates.
}
\usage{
plot(x, CTPM=NULL, col="black", col.ctpm = "red", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A variogram object calculated using \code{\link{variogram}}.}
  \item{CTPM}{A ctpm model object in the same format as the output of ctpm.fit. Not currently supported.}
  \item{col}{Color for the empirical variogram.}
  \item{col.ctpm}{Color for the model. Not currently supported.}
  \item{\dots}{ Additional \code{plot} function parameters.}
  }

\details{
For highly irregularly phylogenetic distances with few species, it may be useful to set \code{complete = FALSE} to coarsen the variogram. When this is the cases, species are binned across lags, with the number of lags estimated using \code{nclass.scott}.
}

\value{
Returns a plot of semi-variance vs. time lag, with the empirical variogram in black and the ctmm semi-variance function in red if specified. zoom includes a log-scale zoom slider to manipulate fraction.
}
\references{
 M. J. Noonan, W. F. Fagan, and C. H. Fleming (in prep.) ``A semi-variance approach to visualising phylogenetic autocorrelation''.
}

\author{ M. J. Noonan, C. H. Fleming.}

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{vignette("variogram")}, \code{\link{variogram}}, \code{\link{nclass.scott}}. }

\examples{
#Load package and data
library(ctpm)
data("moid_traits")
data("musteloids")

#Extract the trait of interest from the full dataset
SSD <- moid_traits$SSD

#Calculate variogram
SVF <- variogram(SSD, musteloids)

#Plot the variogram
plot(SVF)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{variogram}
\alias{variogram}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate an empirical variogram from phylogenetic data}
\description{
This function calculates the empirical variogram of phylogenetic data for visualizing stationary (time-averaged) autocorrelation structure.
One of three weighting algorithms can be used, based upon Noonan, Fagan, and Fleming (in prep).
}
\usage{
variogram(data, phylo, weights = "BM", complete = TRUE, progress = TRUE, level = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A vector of continuous species trait data. This vector needs to be of the same length and same order as phy$tip.label. }
  \item{phylo}{An object of class 'phylo'.}
  \item{weights}{The weights to apply. Can be one of \code{IID}, \code{BM}.}
  \item{complete}{A logical value indicating whether the semi-variance is to be calculated across all possible lags. Defalts to \code{TRUE}}
  \item{progress}{A logical value indicating whether to include a progress bar.}
  \item{level}{Confidence level of confidence bands (95\% default CIs)}.
  }

\details{
For highly irregularly phylogenetic distances with few species, it may be useful to set \code{complete = FALSE} to coarsen the variogram. When this is the cases, species are binned across lags, with the number of lags estimated using \code{nclass.scott}.
}

\value{
Returns a variogram object (class variogram) which is a dataframe containing the time-lag, \code{Distance}, the semi-variance estimate at that lag, \code{Gamma}, and the upper and lower confidence intervals.
}
\references{
 M. J. Noonan, W. F. Fagan, and C. H. Fleming (in prep.) ``A semi-variance approach to visualising phylogenetic autocorrelation''.
}

\author{ M. J. Noonan. }

\note{ Can be slow on very large phylogenies. }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{vignette("variogram")}, \code{\link{plot.variogram}}, \code{\link{nclass.scott}}. }

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

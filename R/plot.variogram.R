
##############
# Functions for plotting variograms

plot.variogram <- function(x, CTPM = NULL, col="black", col.CTPM = "red", fraction = 1, y.units = NULL, ...){
  
  SVF <- ctmm:::new.variogram(x)
  
  SVF@info$axes <- x@info$axes
  
  suppressWarnings(ctmm:::plot.variogram(SVF,
                        CTMM = CTPM,
                        col.CTMM = col.CTPM,
                        fraction = fraction,
                        ...))
}

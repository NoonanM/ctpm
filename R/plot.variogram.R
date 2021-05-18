
##############
# Functions for plotting variograms

plot.variogram <- function(x, CTPM = NULL, col="black", col.CTPM = "red", fraction = 1, y.units = NULL, ...){
  
  
  suppressWarnings(ctmm:::plot.variogram(x,
                        CTMM = CTPM,
                        col.CTMM = col.CTPM,
                        fraction = fraction,
                        ...))
}

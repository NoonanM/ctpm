
##############
# Functions for plotting variograms

# plot.variogram <- function(x, CTPM = NULL, col="black", col.CTPM = "red", fraction = 1, ...){
#   
#   SVF <- ctmm.variogram(x)
#   
#   SVF@info$axes <- x@info$axes
#   
#   attr(SVF,"info")$lags <- attr(x,"info")$lags
#   
#   suppressWarnings(plot.variogram(SVF,
#                                   CTMM = CTPM,
#                                   col.CTMM = col.CTPM,
#                                   fraction = fraction,
#                                   ...))
# }

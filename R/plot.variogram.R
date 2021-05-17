
##############
# Functions for plotting variograms

# svf.func <- function(CTPM)
# {
#   
#   # FIRST CONSTRUCT STANDARD ACF AND ITS PARAMTER GRADIENTS
#   if(any(class(CTPM) == "lm")) # IID
#   {
#     sigma <- summary(CTPM)$sigma^2
#     acf <- function(t){ if(t==0) {1} else {0} }
#     acf.grad <- function(t){ NULL }
#     SVF <- function(t){if(t==0) {0} else {sigma} }
#     COV <- unname(vcov(FIT))
#     # variance of SVF
#     VAR <- function(t)
#     {
#       { 1^2 * COV }
#     }
#     
#   }
#   else if(any(class(CTPM) == "slouch") && length(CTPM$evolpar) == 1) # BM
#   {
#     sigma <- CTPM$evolpar$sigma2_y
#     acf <- function(t){ 1-t }
#     acf.grad <- function(t){ NULL }
#     SVF <- function(t){ t * sigma}
#     COV <- unname(CTPM$beta_primary$vcov)
#     # variance of SVF
#     VAR <- function(t)
#     { t^2 * COV }
#     
#   }
#   else if(any(class(CTPM) == "slouch") && length(CTPM$evolpar) == 2) # OU
#   {
#     sigma <- CTPM$evolpar$vy
#     tau <- CTPM$evolpar$hl
#     acf <- function(t){ exp(-t/tau) }
#     acf.grad <- function(t){ t/tau^2*acf(t) }
#     SVF <- function(t){sigma * (1 - exp(-(t/tau)))}
#     COV <- unname(CTPM$beta_primary$vcov)
#     # variance of SVF
#     VAR <- function(t)
#     { (1 - exp(-(t/tau)))^2 * COV }
#   }
#   
#   # finish off a few pieces
#   ACF <- function(t) { acf(t) }
#   grad <- function(t,...) { c(svf(t)/sigma, -sigma*acf.grad(t)) }
# 
#   
#   # chi-square effective degrees of freedom
#   DOF <- function(t) { return( 2*SVF(t)^2/VAR(t) ) }
#   
#   return(list(svf=SVF,VAR=VAR,DOF=DOF,ACF=ACF))
# }


plot.variogram <- function(x, CTPM = NULL, col="black", col.CTPM = "red", fraction = 1, y.units = NULL, ...){
  
  #Convert to ctmm variogram object
  # SVF <- data.frame(SVF=x$SVF,
  #                   DOF=x$DOF,
  #                   lag=x$lag)
  # 
  # SVF <- ctmm:::new.variogram(SVF)
  # 
  # SVF@info$axes <- x@info$axes
  
  
  ctmm:::plot.variogram(SVF,
                        CTMM = CTPM,
                        col.CTMM = col.CTPM,
                        fraction = fraction,
                        ...)
  
  
  #Old version of the function
  # if(!is.null(ext))
  # {
  #   xlim <- ext$x
  #   ylim <- ext$y
  # }
  # 
  # # empirical variograms
  # x <- ctmm:::listify(x)
  # n <- length(x)
  # 
  # # theoretical models
  # CTPM <- ctmm:::listify(CTPM)
  # m <- length(CTPM)
  # 
  # 
  # 
  # # subset the data if xlim or fraction provided
  # if(!is.null(xlim))
  # {
  #   fraction <- 1 # xlim overrides fraction
  #   x <- lapply(x,function(y){ y[xlim[1]<=y$lag & y$lag<=xlim[2],] })
  # }
  # else
  # {
  #   # maximum lag in data
  #   max.lag <- sapply(x, function(v){ ctmm:::last(v$Distance) } )
  #   max.lag <- max(max.lag,xlim)
  #   # subset fraction of data
  #   max.lag <- fraction*max.lag
  #   
  #   # subset all data to fraction of total period
  #   #if(fraction<1) { x <- lapply(x, function(y) { y[y$Distance<=max.lag,] }) }
  #   
  #   xlim <- c(0,max.lag)
  # }
  # 
  # 
  # # calculate ylimits from all variograms !!!
  # if(is.null(ylim)) {
  #   max.y <- sapply(x, function(v){ max(v$CI_max) } )
  #   ylim <- c(0,max.y)
  # }
  # 
  # 
  # #Axis labels
  # ylab <- expression(paste(gamma, "( ", tau, " )"))
  # 
  # xlab <- paste("Phylogentic Distance (", units, ")", sep = "")
  # 
  # 
  # # fix base plot layer
  # plot(xlim,ylim, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, col=grDevices::rgb(1,1,1,0), ...)
  # 
  # 
  # #ADD CIs
  # polygon(c(x[[1]]$Distance,rev(x[[1]]$Distance)),
  #         c(x[[1]]$CI_max,rev(x[[1]]$CI_min)),
  #         col="grey85",
  #         border = NA)
  # 
  # #ADD POINT ESTIMATE
  # lines(y = x[[1]]$Gamma,
  #       x = x[[1]]$Distance,
  #       col = col)
  # 
  # 
  # 
  # #If there's a model, add it to the plot
  # if(!is.null(CTPM)){
  #   
  #   lag <- x[[1]]$Distance
  #   
  #   # color array for plots
  #   col.ctpm <- array(col.ctpm,m)
  #   
  #   for(i in 1:m){
  #   plot.svf(lag, CTPM[[i]], col=col.ctpm[[i]])
  #   }
  # }
}

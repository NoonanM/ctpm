
##############
# Functions for plotting variograms

#########
# Generate SVF etc

svf.func <- function(CTPM)
{
  
  # FIRST CONSTRUCT STANDARD ACF AND ITS PARAMTER GRADIENTS
  if(any(class(CTPM) == "lm")) # IID
  {
    sigma <- summary(CTPM)$sigma^2
    acf <- function(t){ if(t==0) {1} else {0} }
    acf.grad <- function(t){ NULL }
    SVF <- function(t){if(t==0) {0} else {sigma} }
    COV <- summary(CTPM)$cov.unscaled * summary(CTPM)$sigma^2
    # variance of SVF
    VAR <- function(t)
    {
      { 1^2 * COV }
    }
    
  }
  else if(any(class(CTPM) == "slouch") && length(CTPM$evolpar) == 1) # BM
  {
    sigma <- CTPM$evolpar$sigma2_y
    acf <- function(t){ 1-t }
    acf.grad <- function(t){ NULL }
    SVF <- function(t){ t * sigma}
    COV <- as.numeric(CTPM$beta_primary$vcov * sigma)
    # variance of SVF
    VAR <- function(t)
    { t^2 * COV }
    
  }
  else if(any(class(CTPM) == "slouch") && length(CTPM$evolpar) == 2) # OU
  {
    sigma <- CTPM$evolpar$vy
    tau <- CTPM$evolpar$hl
    acf <- function(t){ exp(-t/tau) }
    acf.grad <- function(t){ t/tau^2*acf(t) }
    SVF <- function(t){sigma * (1 - exp(-(t/tau)))}
    COV <- as.numeric(CTPM$beta_primary$vcov * sigma)
    # variance of SVF
    VAR <- function(t)
    { (1 - exp(-(t/tau)))^2 * COV }
  }
  
  # finish off a few pieces
  ACF <- function(t) { acf(t) }
  grad <- function(t,...) { c(svf(t)/sigma, -sigma*acf.grad(t)) }
  
  
  # variance of SVF
  # VAR <- function(t)
  # {
  #   g <- grad(t)
  #   return( c(g %*% COV %*% g) ) 
  # }
  
  # chi-square effective degrees of freedom
  # DOF <- function(t,error=0) { return( 2*SVF(t,error=error)^2/VAR(t,error=error) ) }
  DOF <- function(t) { return( 2*SVF(t)^2/VAR(t) ) }
  
  return(list(svf=SVF,VAR=VAR,DOF=DOF,ACF=ACF))
}


##########
plot.svf <- function(lag,CTPM,alpha=0.05,col="red",type="l",...)
{
  # changed from max lag to all lags
  # changed from error=0 or number/logical to error=NULL or array
  
  # number of pixels across diagonal of display
  PX <- ceiling(sqrt(sum(grDevices::dev.size("px")^2)))
  
  {
    lag <- seq(0,ctmm:::last(lag),length.out=PX)
    # error <- 0 -> e0
  }
  # else if(all(diff(error[-1])==0))
  # { error <- error[2] -> e0 } # can still plot curve because all errors it the same
  
  SVF <- svf.func(CTPM)
  svf <- SVF$svf
  DOF <- SVF$DOF
  
  # point estimate plot
  # SVF <- Vectorize(function(t,error=e0) { svf(t,error=error) })
  SVF <- Vectorize(function(t) { svf(t) })
  
  lag[1] <- lag[2]/1000 # almost go to origin, but not quite to avoid nugget
  
  # if(length(error)==1) # can plot curve
  { graphics::curve(SVF,from=0,to=ctmm:::last(lag),n=PX,add=TRUE,col=col,...) }
  # else
  # { graphics::points(lag,SVF(lag,error),col=col,type=type,...) }
  
  # confidence intervals if COV provided
    SVF <- Vectorize(function(t){ svf(t) })(lag)

    for(j in 1:length(alpha))
    {
      # dof <- Vectorize(function(t,error=e0) { DOF(t,error=error) })(lag,error)
      dof <- Vectorize(function(t) { DOF(t) })(lag)
      svf.lower <- Vectorize(function(df){ ctmm:::CI.lower(df,alpha[j]) })(dof)
      svf.upper <- Vectorize(function(df){ ctmm:::CI.upper(df,alpha[j]) })(dof)

      graphics::polygon(c(lag,rev(lag)),c(SVF*svf.lower,rev(SVF*svf.upper)),col=ctmm:::malpha(col,0.1/length(alpha)),border=NA)
    }

}


plot.variogram <- function(x, CTPM = NULL, col="black", col.ctpm = "red", xlim=NULL, ylim=NULL, ext=NULL, fraction = 1, units = "Ma",...){
  
  if(!is.null(ext))
  {
    xlim <- ext$x
    ylim <- ext$y
  }
  
  # empirical variograms
  x <- ctmm:::listify(x)
  n <- length(x)
  
  # theoretical models
  CTPM <- ctmm:::listify(CTPM)
  m <- length(CTPM)
  
  
  
  # subset the data if xlim or fraction provided
  if(!is.null(xlim))
  {
    fraction <- 1 # xlim overrides fraction
    x <- lapply(x,function(y){ y[xlim[1]<=y$lag & y$lag<=xlim[2],] })
  }
  else
  {
    # maximum lag in data
    max.lag <- sapply(x, function(v){ ctmm:::last(v$Distance) } )
    max.lag <- max(max.lag,xlim)
    # subset fraction of data
    max.lag <- fraction*max.lag
    
    # subset all data to fraction of total period
    if(fraction<1) { x <- lapply(x, function(y) { y[y$lag<=max.lag,] }) }
    
    xlim <- c(0,max.lag)
  }
  
  
  # calculate ylimits from all variograms !!!
  if(is.null(ylim)) {
    max.y <- sapply(x, function(v){ max(v$CI_max) } )
    ylim <- c(0,max.y)
  }
  
  
  #Axis labels
  ylab <- expression(paste(gamma, "( ", tau, " )"))
  
  xlab <- paste("Phylogentic Distance (", units, ")", sep = "")
  
  
  # fix base plot layer
  plot(xlim,ylim, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, col=grDevices::rgb(1,1,1,0), ...)
  

  #ADD CIs
  polygon(c(x[[1]]$Distance,rev(x[[1]]$Distance)),
          c(x[[1]]$CI_max,rev(x[[1]]$CI_min)),
          col="grey85",
          border = NA)
  
  #ADD POINT ESTIMATE
  lines(y = x[[1]]$Gamma,
        x = x[[1]]$Distance,
        col = col)
  

  
  #If there's a model, add it to the plot
  if(!is.null(CTPM)){
    
    lag <- x[[1]]$Distance
    
    # color array for plots
    col.ctpm <- array(col.ctpm,m)
    
    for(i in 1:m){
    plot.svf(lag, CTPM[[i]], col=col.ctpm[[i]])
    }
  }
}

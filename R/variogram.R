
variogram <- function(data, phylo, weights = "BM", complete = TRUE, progress = TRUE, level = 0.95){
  
  #For testing
  data("moid_traits")
  data("musteloids")
  data <- moid_traits$SSD
  phylo <- musteloids
  
  SPECIES <- phylo$tip.label
  names(data) <- SPECIES
  
  #Get all pairwise phylogenetic distances
  DISTS <- as.data.frame(as.table(cophenetic(phylo)))
  DISTS <- DISTS[order(DISTS$Freq),]
  
  if(complete){
    TAU <- unique(DISTS$Freq)
  } else {
    TAU <- DISTS$Freq
    TAU <- sort(kmeans(TAU,
                       centers = sqrt(length(TAU))+1,
                       iter.max = 100)$centers[,1])
    
    
    # TAU <- seq(min(DISTS$Freq),
    #            max(DISTS$Freq),
    #            length.out = nclass.scott(DISTS$Freq))
  }
  
  
  # Calculate gamma 
  GAMMA <- vector("numeric", length(TAU))
  CI_min <- vector("numeric", length(TAU))
  CI_max <- vector("numeric", length(TAU))
  
  if(progress){
    pb <- txtProgressBar(min = 0, max = length(TAU)-1, style = 3)
  }
  
  for(k in 2:length(TAU)){
    
    #Species with a common ancestor separated by lag tau
    if(complete){
      SUB <- DISTS[DISTS$Freq == TAU[k],]
    } else{
      SUB <- DISTS[DISTS$Freq <= TAU[k] & DISTS$Freq > TAU[k-1],]
    }
    
    #Check for lags without data. Happens when complete = FALSE
    if(nrow(SUB) >0){
      
      
      #######################################################
      ##########        Calculate the weights       #########
      #######################################################
      
      #uniform weights (not statistically appropriate, but I've included them in case they're needed for testing)
      if(weights == "uniform"){
        PAIRS <- expand.grid(SUB$Var1,SUB$Var2)
        COV <- diag(nrow(PAIRS))
        ONE <- rep(1, nrow(COV))
        W <- ctmm:::PDsolve(COV) %*% ONE
        W <- W/sum(W)
      }
      
      
      #IID weights
      if(weights == "IID"){
        PAIRS <- combn(unique(as.character(SUB[,1])),2)
        
        #IID weights
        N <- length(unique(SUB[,1]))
        COV <- array(0,c(N,N,N,N)) # where N <- length(SUBSET)
        # then set the off-diagonals in a loop (which also incorrectly sets the diagonals to 1/2, but that's fine if we do it first)
        for(i in 1:N){
          COV[i,,i,] <- 1/4
          COV[i,,,i] <- 1/4
          COV[,i,i,] <- 1/4
          COV[,i,,i] <- 1/4
        }
        
        # then flatten 
        dim(COV) <- c(N^2,N^2)
        
        # Drop any that aren't needed
        KEEPERS <- match(do.call("paste", SUB[,1:2]), do.call("paste", expand.grid(unique(SUB[,1]),unique(SUB[,1]))))
        COV <- COV[KEEPERS,KEEPERS]
        
        # then set the diagonal to 1
        diag(COV) <- 1
        
        # Vector of 1s
        ONE <- rep(1, nrow(COV))
        
        W <- solve(COV) %*% ONE # can switch to ctmm:::PDsolve if the behaviour of solve is not ideal
        W <- W/sum(W)#} else {
      }
      
      
      
      #BM weights
      if(weights == "BM"){
        
        PAIRS <- combn(unique(as.character(SUB[,1])),2)
        
        #Species at time lag tau
        phylo_sub <- ape::keep.tip(phylo, as.character(SUB$Var1))
        
        VCV <- ape::vcv(phylo_sub)
        
        N <- nrow(unique(SUB))
        COV <- array(0,c(N,N)) # where N <- length(SUBSET)
        
        for(i in 1:nrow(SUB)){
          for(j in 1:nrow(SUB)){
            IJ <- VCV[as.character(SUB[,1][i]),as.character(SUB[,2][i])]
            KL <- VCV[as.character(SUB[,1][j]),as.character(SUB[,2][j])]
            COV[i,j] <- (IJ + KL)/TAU[k]
          } #closes loop over j
        } #closes loop over i
        
        # then set the diagonal
        diag(COV) <- 1
        
        #Squared correlation matrix
        COV <- COV^2
        
        #Vector of 1s
        ONE <- rep(1, nrow(COV))
        
        #Calculate the weights
        W <- ctmm:::PDsolve(COV) %*% ONE # can switch to ctmm:::PDsolve if the behaviour of solve is not ideal
        W <- W/sum(W)
        
      }
      
      
      # Weighted semi-variance
      DIFFS <- vector("numeric", nrow(SUB))
      for(l in 1:length(DIFFS)){
        DIFFS[l] <- W[l]*0.5*abs((data[SUB$Var1[l]] - data[SUB$Var2[l]]))^2
      }
      GAMMA[k] <- sum(DIFFS)
      
      
      #Variance
      VAR <- 2 * (t(W) %*% COV %*% W) * (GAMMA[k])^2
      
      #Chi square CI
      CIs <- ctmm:::chisq.ci(MLE = GAMMA[k], VAR = VAR, level = level)
      
      CI_min[k] <- CIs[1]
      CI_max[k] <- CIs[3]
      
      
      #Update progress bar is turned on
      if(progress){setTxtProgressBar(pb, k)}
    } #Closes the check for any data in the lag
  }
  
  SVF <- data.frame(Distance = TAU,
                    Gamma = GAMMA,
                    CI_min = CI_min,
                    CI_max = CI_max)
  
  SVF <- new.variogram(SVF)
  
  return(SVF)
}


##############
# Function for plotting the variograms

plot.variogram <- function(x, CTPM = NULL, col="black", col.ctpm = "red", units = "Ma",...){
  
  plot(y = x$Gamma,
       x = x$Distance,
       ylab = expression(paste(gamma, "( ", tau, " )")),
       xlab = paste("Phylogentic Distance (", units, ")", sep = ""),
       type = "l",
       col = col,
       ...)
  polygon(c(x$Distance,rev(x$Distance)),
          c(x$CI_max,rev(x$CI_min)),
          col="grey85",
          border = NA)
  lines(y = x$Gamma,
        x = x$Distance,
        col = col)
  
  #If there's a model, add it to the plot
  if(!is.null(CTPM)){
    
    #Range of the variogram
    RANGE <- c(0, max(x$Dist))
    
    #Sequence of time lags to plot over for lines()
    TAU <- seq(RANGE[1], RANGE[2], 0.001)
    
    # Check which model we're working with and plot accordingly
    
    #IID
    if(any(class(CTPM) == "gls")){
      lines(y = rep(CTPM$sigma^2,
                    length(TAU+1)),
            x =  TAU,
            col = col.ctpm)
      lines(y = seq(0, CTPM$sigma^2,
                    length.out = 5),
            x =  rep(0,5), col = col.ctpm)
    }
    
    #BM
    if(any(class(CTPM) == "slouch") && length(CTPM$evolpar) == 1){
      lines(y = (TAU*CTPM$evolpar$sigma2_y), x = TAU, col = col.ctpm)
    }
    
    #OU
    if(any(class(CTPM) == "slouch") && length(CTPM$evolpar) == 2){
      lines(y = CTPM$evolpar$vy * (1 - exp(-(TAU/CTPM$evolpar$hl))),
            x =  TAU, col = col.ctpm)
    }
  }
}

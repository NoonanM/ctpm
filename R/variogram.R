
variogram <- function(data, phylo, weights = "IID", complete = FALSE, units = "Ma", progress = TRUE,  ...){
  
  #For testing
  # data("moid_traits")
  # data("musteloids")
  # data <- moid_traits$SSD
  # phylo <- musteloids
  
  SPECIES <- phylo$tip.label
  names(data) <- SPECIES
  
  #Get all pairwise phylogenetic distances
  DISTS <- as.data.frame(as.table(ape::cophenetic.phylo(phylo)))
  DISTS <- DISTS[order(DISTS$Freq),]
  
  #Lag BINS
  if(complete){
    TAU <- unique(DISTS$Freq)
  } else {
    TAU <- DISTS$Freq
    TAU <- sort(kmeans(TAU,
                       centers = sqrt(length(TAU))+1, ...)$centers[,1])
  }
  
  
  # Calculate gamma 
  GAMMA <- vector("numeric", length(TAU))
  CI_min <- vector("numeric", length(TAU))
  CI_max <- vector("numeric", length(TAU))
  VAR <- vector("numeric", length(TAU))
  DOF <- vector("numeric", length(TAU))
  
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
      VAR[k] <- 2 * (t(W) %*% COV %*% W) * (GAMMA[k])^2
      
      DOF[k] <- 2*GAMMA[k]^2/VAR[k]
      
      #Chi square CI
      # CIs <- ctmm:::chisq.ci(MLE = GAMMA[k], VAR = VAR, level = level)
      # 
      # CI_min[k] <- CIs[1]
      # CI_max[k] <- CIs[3]
      
      
      #Update progress bar is turned on
      if(progress){setTxtProgressBar(pb, k)}
    } #Closes the check for any data in the lag
  }
  
  #Convert TAU for correct plotting
  LAG <- TAU %#% units
  
  SVF <- data.frame(SVF=GAMMA,DOF=DOF,lag=LAG)
  
  # SVF <- data.frame(Distance = TAU,
  #                   Gamma = GAMMA,
  #                   CI_min = CI_min,
  #                   CI_max = CI_max)
  
  SVF <- new.variogram(SVF)
  
  SVF@info$axes <- "x"
  
  return(SVF)
}

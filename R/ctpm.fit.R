ctpm.fit <- function(data, phylo, model = NULL, units = "Ma"){
  
  #Fit an IID model
  if(model == "IID"){
    fit <- lm(data ~ 1)
    COV <- matrix(vcov(fit))
    row.names(COV) <- "major"
    colnames(COV) <- "major"
    
    FIT <- ctmm::ctmm(mean = "stationary",
                      sigma =  var(data),
                      errors = FALSE,
                      mu = fit$coefficients[1],
                      tau = NULL,
                      axes = c("x"),
                      COV.mu = vcov(fit),
                      COV = COV,
                      circle = 0,
                      range = TRUE,
                      AIC = AIC(fit),
                      BIC = BIC(fit))
  }
  
  
  #Fit the BM model
  if(model == "BM"){
    
    fit <- slouch::brown.fit(phy = phylo,
                            species = phylo$tip.label,
                            response = data,
                            hessian = T)
    
    tau <- Inf; names(tau) <- "position"
    
    COV <- -solve(fit$hessian) %#% units
    row.names(COV) <- "major"
    colnames(COV) <- "major"
    
    FIT <- ctmm::ctmm(mean = "stationary",
                      sigma =  var(data),
                      errors = FALSE,
                      mu = fit$beta_primary$coefficients[1],
                      tau = tau,
                      axes = c("x"),
                      COV.mu = matrix(Inf),
                      COV = COV,
                      circle = 0,
                      range = FALSE,
                      AIC = fit$modfit$AIC,
                      BIC = NULL)
    
  }
  
  #Fit the OU model
  if(model == "OU"){
    
    phylo$edge.length <- phylo$edge.length %#% units
    
    fit <- slouch::slouch.fit(phy = phylo,
                             species = phylo$tip.label,
                             response = data,
                             hessian = T)
    
    tau <- fit$evolpar$hl; names(tau) <- "position"
    
    tau <- tau %#% units
    
    COV <- -solve(fit$hessian) %#% units
    
    #re-arrange
    COV_2 <- COV
    COV[1,1] <- COV_2[2,2]; COV[2,2] <- COV_2[1,1]
    
    row.names(COV) <- c("major", "tau position")
    colnames(COV) <- c("major", "tau position")
    
    FIT <- ctmm::ctmm(mean = "stationary",
                      sigma = fit$evolpar$vy,
                      errors = FALSE,
                      mu = fit$beta_primary$coefficients[1],
                      tau = tau,
                      axes = c("x"),
                      COV.mu = matrix(fit$beta_primary$vcov),
                      COV = COV,
                      circle = 0,
                      range = TRUE,
                      AIC = fit$modfit$AIC,
                      BIC = NULL)
    
  }
  
  #Return the fitted model
  return(FIT)
}

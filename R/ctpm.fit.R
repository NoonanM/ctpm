ctpm.fit <- function(data, phylo, model = NULL){
  
  #Fit an IID model
  if(model == "IID"){
    FIT <- nlme::gls(data ~ 1)
  }
  
  
  #Fit the BM model
  if(model == "BM"){
    FIT <- slouch::brown.fit(phy = phylo,
                            species = phylo$tip.label,
                            response = data)
  }
  
  #Fit the OU model
  if(model == "OU"){
    FIT <- slouch::slouch.fit(phy = phylo,
                             species = phylo$tip.label,
                             response = data)
  }
  
  #Return the fitted model
  return(FIT)
}

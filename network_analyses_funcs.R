#### Functions needed to run the network analyses
#### Project: gallers
#### First edit: 2020.04.06
#### Latest edit: 2021.10.25
#### Author: Kate Wootton
#### R version 3.6.3

library(vegan)

# Get the fomula for traits in RandomForest
get_RF_traits_formula <- function(alienData){
  # Base of formula
  fromTraits <- paste(colnames(alienData$traitFrom), collapse = "+")
  toTraits <- paste(colnames(alienData$traitTo), collapse = "+")
  
  # Formula
  FormulaFromTo <- as.formula(paste("~", fromTraits, "+", toTraits))
  
  return(FormulaFromTo)
}


# Get the fomula for traits in GLM or fourth corner
# traits2use is a vector of names of traits to use in the analysis OR the string "all" to use all traits
get_GLM_traits_formula <- function(alienData,
                                   traits2use = "all",
                                   which_model = c("GLM","4C")){

  if(traits2use == "all"){
    traits2use = c(colnames(alienData$traitFrom),colnames(alienData$traitTo))
  }
  fromTraits <- colnames(alienData$traitFrom)[which(colnames(alienData$traitFrom) %in% traits2use)]
  
  fromTraitsSq <- sapply(fromTraits, 
                         function(i) 
                           ifelse(is.numeric(unlist(alienData$traitFrom[i])), 
                                  return(paste0(i,"_Sq")), 
                                  return(NULL))) %>%
    unlist()
  
  toTraits <- colnames(alienData$traitTo)[which(colnames(alienData$traitTo) %in% traits2use)]
  
  toTraitsSq <- sapply(toTraits, 
                         function(i) 
                           ifelse(is.numeric(unlist(alienData$traitTo[i])), 
                                  return(paste0(i,"_Sq")), 
                                  return(NULL))) %>%
    unlist()
  # Independent term
  indep <- paste(c(fromTraits, toTraits), collapse = "+")
  
  if(which_model == "GLM"){
    combTr <- combn(c(fromTraits, toTraits), 2)  
  }
  if(which_model == "4C"){
    combTr <- t(expand.grid(fromTraits, toTraits))
  }
  
  inter <- paste(combTr[1,], ":", combTr[2,], collapse = "+")
  
  # Build formula
  FormulaFromTo <- as.formula(paste("~", indep, "+", inter))
  
  return(FormulaFromTo)
}
  
# Get phylogeny "traits" for RandomForest
get_RF_phylo_traits <- function(alienData){
  fromPhylo <- wcmdscale(alienData$phyloDistFrom, add = TRUE)
  colnames(fromPhylo) <- paste0("fromPhylo", 1:ncol(fromPhylo))
  
  toPhylo <- wcmdscale(alienData$phyloDistTo, add = TRUE)
  colnames(toPhylo) <- paste0("toPhylo", 1:ncol(toPhylo))
  
  fromToPhylo <- alienData(adjMat = alienData$adjMat,
                             traitFrom = as.data.frame(fromPhylo),
                             traitTo = as.data.frame(toPhylo))
  
  return(fromToPhylo)
}

# Get the phylogeny formula for RandomForest
get_RF_phylo_formula <- function(fromToPhylo){
  # Base of formula
  fromPhyloTraits <- paste(colnames(fromToPhylo$traitFrom), collapse = "+")
  toPhyloTraits <- paste(colnames(fromToPhylo$traitTo), collapse = "+")
  
  # Formula
  FormulaPhylo <- as.formula(paste("~", fromPhyloTraits, "+", toPhyloTraits))
  
  
  return(FormulaPhylo)
}
  


  #### Function to plot interaction, phylogenetic, and trait distances against each other and calculate correlations
  #### Project: gallers
  #### First edit: 2020.04.04
  #### Latest edit: 2021.10.25
  #### Author: Kate Wootton
  #### R version 3.6.3


# Load libraries
library(psych)
library(reshape2)
library(dplyr)

# Calculate distances between species based on shared interactions
calc_dist_ints <- function(alienData,
                           troph = c("from", "to"),
                           vegMethod = "jaccard", 
                           toMelt = FALSE){
  if(troph == "to"){
    distInts <- vegan::vegdist(t(alienData$adjMat), 
                               method = vegMethod,
                               na.rm = TRUE)
    distInts <- as.matrix(distInts)   
  }
  
  if(troph == "from"){
    distInts <- vegan::vegdist(alienData$adjMat, 
                               method = vegMethod,
                               na.rm = TRUE)
    distInts <- as.matrix(distInts)   
  }
  
  if(toMelt){
    distInts <- melt(as.matrix(distInts), value.name = "ints")
  } 
  
  return(distInts)
}


# Calculate distances between species based on traits
calc_dist_traits <- function(formulaTrait = ~ -1 +., 
                             alienData, 
                             troph = c("from", "to"),
                             vegMethod = "euclidean", 
                             toMelt = FALSE){
  
  if(troph == "from"){
    trophTraits <- stats::model.matrix(formulaTrait, data = alienData$traitFrom)
  }

  if(troph == "to"){
    trophTraits <- stats::model.matrix(formulaTrait, data = alienData$traitTo)
  }
  
  distTraits <- vegan::vegdist(trophTraits, method = "euclidean", na.rm = TRUE)
  distTraits <- as.matrix(distTraits)

  if(toMelt){
    distTraits <- melt(distTraits, value.name = "traits")
  } 
  
  return(distTraits)
}

# Calculate correlations in distances based on traits, interactions, and phylogenies
# Option to plot correlations
compare_distances <- function(formulaTrait = ~ -1 +.,
                              alienData, 
                              troph = c("from", "to"),
                              vegMethodSp = "jaccard",
                              vegMethodTr = "euclidean", 
                              toMelt = TRUE,
                              toPlot = FALSE, 
                              maintitle=NA){
  
  # Calculate distances based on interactions
  distInts <- calc_dist_ints(alienData, troph = troph, 
                             vegMethod = vegMethodSp, toMelt = toMelt)
  
  # Calculate trait distances
  distTraits <- calc_dist_traits(formulaTrait = formulaTrait,
                                 alienData, troph = troph,
                                 vegMethod = vegMethodTr, toMelt = toMelt)
  
  # Include phylogenetic distances
  if(troph == "from"){
    distPhylo <- melt(as.matrix(alienData$phyloDistFrom),
                      value.name = "phylo")
    dists <- left_join(distPhylo,distInts) %>% left_join(distTraits)
  }
 
  if(troph == "to"){
    distPhylo <- melt(as.matrix(alienData$phyloDistTo),
                      value.name = "phylo")
    dists <- left_join(distPhylo,distInts) %>% left_join(distTraits)
  }
  
  #Plot correlations
  if(toPlot){
    pairs.panels(dists[,3:5], main = maintitle)
  }
  
  # calculate correlations
  distCors <- cor(dists[,3:5], use = "complete.obs")
  
  return(distCors)
}


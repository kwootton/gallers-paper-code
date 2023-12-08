#### Functions to edit data
#### "Differential imprints of trait-matching in a network: complementary insights from complementary methods"
#### Wootton et al 2024
#### Project: gallers
#### First edit: 2023.02.06
#### Latest edit: 2023.12.08
#### Author: Kate Wootton
#### R version 4.2.2


#####################################
## Set up 
#####################################

# load packages

trimTraits <- function(data){
  # Set traits we want to use
  # We are removing TREE.HEIGHT because it is a component of tree volume, and correlated with leaf toughness
  salTraits <- c("LEAF.SIZE", "LEAF.HAIRINESS", "LEAF.THICKNESS",
                 "GLUCOSIDES", "LEAF.TOUGHN", "LEAF.SLA", 
                 "TREE.VOLUME", "PHENOLOGY.SAL")
  # We are removing GALL.TOUGHN because it is correlated with GALL.WALL (and GALL.WALL is easier to measure),
  # OVERWINTERING.SITE because only shoot galls and a few bud gall species overwinter in the gall, the rest are in the soil
  # REPRODUCTION.GAL because only 5 species reproduce parthenogenetically and they are all leaf blade bean galls
  galTraits <- c("GALLTYPE", "BODYLENGTH.GAL", "PHENOLOGY.GAL",
                 "OVIPOS.STRATEGY", "DEVELOPMENT.TIME", "GALL.VOLUME",
                 "GAL.POS", "GALL.WALL")
  # We are removing "P.I" (Parasitoid/inquiline) because inquiline is also a category of endo/ecto
  # KOINO.IDIO because it's almost perfectly correlated with ENDO.ECTO
  # and ATTACK.STRATEGY because inquilines are the only chewers
  parTraits <- c("ENDO.ECTO", "ATTACK.STAGE", "PHENOLOGY.PAR",
                 "BODYLENGTH.PAR", "OVIPOS.LNTH")
  
 if("GALLTYPE" %in% colnames(data$traitTo)){
   traits2keep <- c(salTraits, galTraits)
  }else{
    traits2keep <- c(galTraits, parTraits)
  }
  
  data <- removeTraits(data, traits2keep)
  
  return(data)
}

removeTraits <- function(data, traits2keep){
  data$traitFrom <- data$traitFrom[,which(colnames(data$traitFrom) %in% traits2keep)]
  data$traitTo <- data$traitTo[,which(colnames(data$traitTo) %in% traits2keep)]
  
  return(data)
}

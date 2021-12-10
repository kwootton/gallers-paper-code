#### Function to run each model with all combination of n traits
#### Project: gallers
#### First edit: 2020.04.04
#### Latest edit: 2021.10.25
#### Author: Kate Wootton
#### R version 3.6.3


# Run each model with ntraits
run_trait_combinations <- function(alienData, ntraits = 4, to_save = FALSE){
  
  # Get all names of all possible traits
  from_traits <- colnames(alienData$traitFrom)
  to_traits <- colnames(alienData$traitTo)
  
  # Get the number of unique trait values for each trait (e.g. "P.I." has 2 - P and I)
  from_traits_nvals <- alienData$traitFrom %>% apply(2, function(x) length(unique(x)))
  to_traits_nvals   <- alienData$traitTo %>% apply(2, function(x) length(unique(x)))
  traits_nvals <- c(from_traits_nvals, to_traits_nvals)
  
  # Get all combinations of ntraits, where there's at least one trait from each trophic group
  trait_combos <- t(combn(c(from_traits, to_traits), ntraits)) %>%
    tibble %>%
    mutate(both_levels = apply(X = ., MARGIN = 1, function(x) any(x %in% from_traits) & any(x %in% to_traits))) %>%
    filter(both_levels == TRUE) %>%
    select(!last_col())
  
  
  trait_results <- list()
  iter <- 0
  
  # Run through all trait combinations and run the models
  for(i in 1:nrow(trait_combos)){
    
    trait_subs <- trait_combos[i,] %>% unlist                                  # Which traits for this combination
    from_traits_subs <- trait_subs[which(trait_subs %in% from_traits)]         # Which are from traits
    to_traits_subs   <- trait_subs[which(trait_subs %in% to_traits)]           # and which are to traits
    
    iter <- iter + 1
    
    alienData2 <- alienData                                                    # Set up a new aliendata object 
    alienData2$traitFrom <-alienData$traitFrom[,from_traits_subs,drop=FALSE]   # and only keep those traits in our subset
    alienData2$traitTo <-alienData$traitTo[,to_traits_subs,drop=FALSE]
    
    # run the models that use real traits
    KNN_traits <- fitKNN(alienData2,
                         distFrom = "jaccard",
                         distTo = "jaccard",
                         distTraitFrom = "euclidean",
                         distTraitTo = "euclidean", 
                         weight = 1, 
                         nNeig = 3)
    
    
    fourC <- fitGLM(polyTrait(alienData2),
                    formula = get_GLM_traits_formula(alienData2, traits2use = "all", which_model = "4C"),
                    family = binomial(link = "logit"))
    
    
    GLM <- fitGLM(polyTrait(alienData2),
                  formula = get_GLM_traits_formula(alienData2,traits2use = "all",which_model = "GLM"),
                  family = binomial(link = "logit"))
    
    RF <- fitRF(alienData2,
                formula = get_RF_traits_formula(alienData2),
                ntree = 2000,
                nodesize = 1)
    
    # Get the fits
    fits <- data.frame()
    fits["ll","KNN_traits"] <- logLik(KNN_traits, error = 0.001)
    fits["ll","fourC"] <- logLik(fourC, error = 0.001)
    fits["ll","GLM"] <- logLik(GLM, error = 0.001)
    fits["ll","RF"] <- logLik(RF, error = 0.001)
    
    fits["tjurD","KNN_traits"] <- tjur(KNN_traits)
    fits["tjurD","fourC"] <- tjur(fourC)
    fits["tjurD","GLM"] <- tjur(GLM)
    fits["tjurD","RF"] <- tjur(RF)
    
    trait_results[[i]] <- list()
    trait_results[[i]][["traits"]] <- trait_subs
    trait_results[[i]][["fits"]] <- fits
  }
  
  # Rbind all the traits into a dataframe
  traits_fits <- lapply(trait_results, function(x) x[["traits"]]) %>% 
    do.call("rbind",.) 
  colnames(traits_fits) <- paste0("t",seq(1:ntraits))
  
  # Get the fits associated with each model
  KNN_tjurD <- unlist(lapply(trait_results, function(x) x[["fits"]]$KNN_traits[2]))
  fourC_tjurD <- unlist(lapply(trait_results, function(x) x[["fits"]]$fourC[2]))
  GLM_tjurD <- unlist(lapply(trait_results, function(x) x[["fits"]]$GLM[2]))
  RF_tjurD <- unlist(lapply(trait_results, function(x) x[["fits"]]$RF[2]))
  
  # And put together in one big dataframe
  fits <- data.frame(traits_fits,KNN_tjurD = KNN_tjurD, fourC_tjurD = fourC_tjurD, GLM_tjurD = GLM_tjurD, RF_tjurD = RF_tjurD) 
  fits$total.ntraits <- apply(fits, 1, function(x) traits_nvals[x["t1"]]+traits_nvals[x["t2"]]+traits_nvals[x["t3"]]+traits_nvals[x["t4"]])
  colnames(fits)[5:8] <- c("KNN","4-corner","GLM","RF")
  
  
  # Is Salix or galler the from species?
  from_sp <- ifelse("GALLTYPE" %in% from_traits, "gp", "sg")
  
  if(to_save){
    save(fits, file=paste0("submission_202105/data/traits",ntraits,"_modelfits_",from_sp,".Rdata"))
  }
  
  return(fits)
}

# scale to 2 dp
scaleFUN <- function(x) sprintf("%.2f", x)

# Plot the fits of each model with different traits as a box plot faceted by model
plot_ntrait_fits <- function(alienData,fits,ntraits){
  
  from_traits <- colnames(alienData$traitFrom)
  to_traits <- colnames(alienData$traitTo)
  
  # Is Salix or galler the from species?
  from_sp <- ifelse("GALLTYPE" %in% from_traits, "G", "S")
  
  # Make a list of all the fits for each combination that includes each trait "t"
  trait_scores <- list()
  for(t in c(from_traits,to_traits)){
    trait_scores[[t]] <- fits[which(fits[,1:ntraits] == t, arr.ind=TRUE)[,1],] %>% mutate(trait = t)
  }
  trait_scores_df <- do.call(rbind, trait_scores)
  

  trait_scores_long <- trait_scores_df %>%
    pivot_longer(KNN:RF, names_to = "Model", values_to = "tjurD") %>%
    mutate(troph = ifelse(trait %in% from_traits, "from","too")) %>%
    mutate(trait = factor(trait, levels = c(from_traits[order(from_traits,decreasing=TRUE)],to_traits[order(to_traits,decreasing=TRUE)])))
  
  # Change the names of traits to make them friendlier to read
  if(from_sp == "G"){
    trait_scores_long <- trait_scores_long %>%
      mutate(trait = dplyr::recode(trait, GALL.WALL = "Gall wall", GALL.POS = "Gall position", GALL.TOUGHN = "Gall toughness", GALL.VOLUME = "Gall volume", KOINO.IDIO = "Koino/idiobiont", GALLTYPE = "Gall type", BODYLENGTH.GAL = "Bodylength (G)", PHENOLOGY.GAL = "Phenology (G)", OVERWINTERING.SITE = "Overwintering site", OVIPOS.STRATEGY = "Oviposition strategy", REPRODUCTION.GAL = "Reproduction (G)",DEVELOPMENT.TIME = "Development time", ATTACK.STRATEGY = "Attack strategy", BODYLENGTH.PAR = "Bodylength (P)", OVIPOS.LNTH = "Ovipositor length", ENDO.ECTO = "Endo/Ectoparasitoid", ATTACK.STAGE = "Attack stage", P.I = "Parasitoid/Inquiline", PHENOLOGY.PAR = "Phenology (P)"))
  }else{
    trait_scores_long <- trait_scores_long %>%
      mutate(trait = dplyr::recode(trait, GALL.WALL = "Gall wall", GALL.POS = "Gall position", GALL.TOUGHN = "Gall toughness", GALL.VOLUME = "Gall volume", GALLTYPE = "Gall type", BODYLENGTH.GAL = "Bodylength", PHENOLOGY.GAL = "Phenology (G)", OVERWINTERING.SITE = "Overwintering site", OVIPOS.STRATEGY = "Oviposition strategy", REPRODUCTION.GAL = "Reproduction",DEVELOPMENT.TIME = "Development time", LEAF.SIZE = "Leaf size", LEAF.THICKNESS = "Leaf thickness", LEAF.HAIRINESS = "Leaf hairiness", GLUCOSIDES = "Glucosides", LEAF.TOUGHN = "Leaf toughness", LEAF.SLA = "Leaf SLA", TREE.HEIGHT = "Tree height", TREE.VOLUME = "Tree volume", PHENOLOGY.SAL = "Phenology (S)")) 
  }
  
  # Set colors based on whether it's salix-galler or gal-para we're plotting
  if(from_sp == "G"){
    cols = c("darkorange","darkblue")
  }else{
    cols = c("darkgreen","darkorange")
  }
  
  # plot
  p4_box <- ggplot(trait_scores_long) + geom_boxplot(aes(x = trait, y = tjurD, group = trait, fill = troph), alpha = 0.5,outlier.alpha = 0.2, outlier.size=0.5)+
    scale_y_continuous(labels=scaleFUN) + coord_flip()+scale_fill_manual(values = cols) + ylab("Model fit (Tjur's D)")+ xlab("")  +
    theme_pubr(legend = "none") +
    facet_grid(.~Model, scales = "free_x")+ theme(panel.spacing.x = unit(8,"mm"))
  
  return(p4_box)
}

#### Functions needed in the analysis of the galler data
#### Project: gallers
#### First edit: 2020.04.06
#### Latest edit: 2021.10.25
#### Author: Kate Wootton
#### R version 3.6.3

## Load libraries
library(magrittr)
library(dplyr)
library(ape) # for function GetSalixPhylo

# reads in and cleans full interaction data
GetInteractionData <- function(datadir){
  # This is the adjacency matrix of gallers and parasitoids plus info on sites and salix species
  interaction_data <- read.csv2(paste0(datadir,"interaction_data.csv"), stringsAsFactors = FALSE, row.names = NULL)
  # Remove the species "alpina?" because it's uncertain and we don't have enough info on it
  interaction_data <- interaction_data[which(interaction_data$Salix != "alpina?"),]
  # let's be consistent with name of column RGALLER
  colnames(interaction_data)[which(colnames(interaction_data) == "Rgaller")] <- "RGALLER"
  
  return(interaction_data)
}


# Reads in and cleans full interaction data, but only returns species that were found in the field when collecting traits 2018
GetFieldInteractionData <- function(datadir){
  # This is the adjacency matrix of gallers and parasitoids plus info on sites and salix species
  interaction_data <- read.csv2(paste0(datadir,"interaction_data.csv"), stringsAsFactors = FALSE, row.names = NULL)
  # Remove the species "alpina?" because it's uncertain and we don't have enough info on it
  interaction_data <- interaction_data[which(interaction_data$Salix != "alpina?"),]
  # let's be consistent with name of column RGALLER
  colnames(interaction_data)[which(colnames(interaction_data) == "Rgaller")] <- "RGALLER"
  
  # If we only care about species found in the field, cut out the rest
  fieldspp <- read.csv(paste0(datadir,"list_of_field_species.csv"))
  interaction_datafield <- interaction_data %>% filter(RGALLER %in% fieldspp$Species) %>% filter(RSAL %in% fieldspp$Species)# %>% select(c(1:5), which(colSums(.[6:ncol(.)], na.rm = TRUE) >0))
  interaction_data <- cbind(interaction_datafield[,c(1:5)],interaction_datafield[,(which(colSums(interaction_datafield[6:ncol(interaction_datafield)], na.rm = TRUE) > 0)+5)])
  
  return(interaction_data)
}

# Find which species we have fewer than "rarespp_threshold" records in the data
# (where a "record" is a line in the df_interact.csv file, which corresponds to a unique combination
# of site visit, salix, galler, and parasitoid species)
# Note there are many other ways to calculate rarespp and this is the most conservative
FindRareSpecies <- function(rarespp_threshold,taxa){
  
  df_interact <- read.csv2(paste0(datadir,"df_interact.csv"))
  species     <- unique(df_interact[,taxa])
  nrecords    <- sapply(species, function(spp) length(which(df_interact[,taxa]==spp)))
  rarespp     <- species[nrecords <= rarespp_threshold] %>% as.character()
  
  return(rarespp)
}

# RemoveRareSpecies <- function(interaction_data, rarespp_threshold){
#   
#   # Remove species we don't see frequently
#   rare_salix   <- FindRareSpecies(rarespp_threshold, "RSAL")
#   rare_gallers <- FindRareSpecies(rarespp_threshold, "RGALLER")
#   rare_paras   <- FindRareSpecies(rarespp_threshold, "RPAR")
#   
#   interaction_data2 <- interaction_data %>%
#                       filter(!RSAL %in% rare_salix) %>%
#                       filter(!RGALLER %in% rare_gallers) %>%
#                       select(!all_of(rare_paras))
#   print(paste0("WARNING!!! REMOVING: ",length(rare_salix)," Salix spp, ", length(rare_gallers)," gallers and ", length(rare_paras), " parasitoids"))
#   print(paste0("Because we have ",rarespp_threshold, " or fewer records for them. Set with rarespp_threshold"))
#   
#   return(interaction_data2)   
# }

RemoveRareSpecies <- function(interaction_data, rarespp_threshold){

  # Remove species we don't see frequently
  rare_salix   <- FindRareSpecies(rarespp_threshold, "RSAL")
  rare_gallers <- FindRareSpecies(rarespp_threshold, "RGALLER")
  rare_paras   <- FindRareSpecies(rarespp_threshold, "RPAR")

  interaction_data2 <- interaction_data %>%
                      filter(!RSAL %in% rare_salix) %>%
                      filter(!RGALLER %in% rare_gallers) %>%
                      select(!all_of(rare_paras))
  print(paste0("WARNING!!! REMOVING: ",length(rare_salix)," Salix spp, ", length(rare_gallers)," gallers and ", length(rare_paras), " parasitoids"))
  print(paste0("Because we have ",rarespp_threshold, " or fewer records for them. Set with rarespp_threshold"))

  return(interaction_data2)
}

# Read in trait data, scale or impute if required
GetTraits <- function(taxa, trait_names, datadir, to_scale=FALSE, to_impute=FALSE){
  if(!taxa %in% c("RSAL","RGALLER","RPAR")){
    print("ERROR! taxa must be RSAL, RGALLER or RPAR")
    break
  }
  # read in all our data
  if(taxa == "RSAL"){
    # read in all our data
    trait_df <- read.csv(paste0(datadir,"salix_traits.csv"),stringsAsFactors = FALSE) %>%
      filter(SPECIES != "alpina")  # Remove alpina because it's not in our main dataset
    colnames(trait_df)[which(colnames(trait_df) == "SPECIES")] <- "Salix"
  }
  if(taxa == "RGALLER"){
    trait_df  <- read.csv(paste0(datadir,"gall_traits.csv"),stringsAsFactors = FALSE) 
  }
  if(taxa == "RPAR"){
    trait_df  <- read.csv(paste0(datadir,"para_traits.csv"),stringsAsFactors = FALSE)
    # Make the trait OVIPOS.LNTH include rostrum length as a proxy for weevils
    trait_df$OVIPOS.LNTH <- as.numeric(as.character(trait_df$OVIPOS.LNTH))
    print("Getting NAs here is normal!")
    trait_df[which(!is.na(trait_df$ROSTRUM)),"OVIPOS.LNTH"] <- trait_df[which(!is.na(trait_df$ROSTRUM)),"ROSTRUM"]
  }
  if(to_impute){
    # fill in missing traits
    trait_df <- ImputeMissing(trait_df, trait_names, taxa)
  }
  if(to_scale){
    # scale all numeric columns
    trait_df <- trait_df %>% mutate_at(all_of(trait_names), ScaleNumeric)
  }
  
  #Make non-numeric columns to factors
  trait_df[sapply(trait_df, is.character)] <- lapply(trait_df[sapply(trait_df, is.character)], 
                                         as.factor)
  
  row.names(trait_df) <- trait_df[,taxa]
  # select just the columns with the trait values and the species names
  trait_df <- trait_df %>% select(all_of(trait_names))
  if(ncol(trait_df) < length(trait_names)){
    print("Warning: requesting traits that aren't in the data")
  }
  return(trait_df)
}


GetAdjacencyMatrix <- function(interaction_data, resource, consumer, edgelist = FALSE, weighted = FALSE){
  resource_spp <- unique(interaction_data[,resource])
  if(consumer == "RGALLER")  consumer_spp <- unique(interaction_data[,consumer])
  if(consumer == "RPAR") consumer_spp <- colnames(interaction_data)[6:ncol(interaction_data)]
  
  matr <- matrix(data = NA, nrow = length(resource_spp), ncol = length(consumer_spp))
  coexist <- as.data.frame(matr, row.names = resource_spp)
  interact <- as.data.frame(matr, row.names = resource_spp)
  colnames(coexist) <- consumer_spp
  colnames(interact) <- consumer_spp
  
  for(resSp in resource_spp){
    resource_subs <- interaction_data[which(interaction_data[,resource] == resSp),]
    
    for(conSp in consumer_spp){
      if(consumer == "RGALLER"){
        consumer_subs <- interaction_data[which(interaction_data[,consumer] == conSp),]
        ints <- interaction_data[which(interaction_data[,resource] == resSp & interaction_data[,consumer] == conSp),]
        coexist[resSp,conSp] <- length(intersect(unique(resource_subs$site), unique(consumer_subs$site)))
        interact[resSp,conSp] <- length(unique(ints$site))
      }
      
      if(consumer == "RPAR"){
        coexist[resSp,conSp] <- length(intersect(unique(resource_subs$site), unique(interaction_data[which(!is.na(interaction_data[,conSp])),"site"])))
        interact[resSp,conSp] <- length(unique(resource_subs[which(!is.na(resource_subs[,conSp])),"site"]))
      }
      
      
      
    }
  }
  matr <- interact
  matr[coexist == 0] <- NA
  
  if(!weighted){
    matr[matr>0] <- 1
  }
  
  if(edgelist){
    el <- im2el(matr, weighted = FALSE) 
    return(el)
  } else{
    return(as.matrix(matr))
  }
  
}



ImputeMissing <- function(trait_df, trait_names, taxa, verbose = FALSE){
  # For missing traits, apply the mean for that group (starting with genus and working up)
  # For strings, assign the most frequent trait
  
  # Start with the lowest taxa level we have for this group, work up
  if(taxa == "RSAL")   levels <- "all"
  if(taxa == "RGALLER") levels <- c("SPECIES.GROUP.L","GENUS.K","all")
  if(taxa == "RPAR")   levels <- c("GENUS","FAMILY","SUPERFAMILY","ORDER", "all")
  
  # run through each column of traits
  for(trait in trait_names){
    #print(trait)
    # Start by trying to group by the lowest taxonomic level
    for(level in levels){
      #print(level)
      # if we have to take the mean of everything:
      if(level == "all"){
        # if the trait is numeric, take the mean
        if(is.numeric(trait_df[,trait])){
          meanval <- mean(trait_df[,trait], na.rm = TRUE)
          trait_df[which(is.na(trait_df[,trait])),trait] <- meanval
        } else{
          # otherwise take the most common string
          tbl <- table(trait_df[,trait])
          modeval <- names(tbl[which(tbl == max(tbl))])
          trait_df[which(is.na(trait_df[,trait])),trait] <- modeval
        }  
        
      } else {
        # group by the relevant taxa level
        if(is.numeric(trait_df[,trait])){
          # if numeric, take the mean
          #trait.means <- trait_df %>% dplyr::group_by_(.dots = level) %>%
          #                dplyr::summarize(meantrait = mean(trait, na.rm = TRUE))
          trait.means <- vector(mode = "numeric", length = length(unique(trait_df[,level])))
          names(trait.means) =unique(trait_df[,level])
          for(grp in unique(trait_df[,level])){
            trait.means[grp] <- mean(trait_df[which(trait_df[,level] == grp),trait])
          }
          for(r in 1:nrow(trait_df)){
            # run through each row in the data frame
            if(is.na(trait_df[r,trait])){
              # and if it's na, find the mean of that taxa level
              #if(!is.na(trait.means[which(trait.means[level] == trait_df[r, level]),"meantrait"])){
              if(!is.na(trait.means[trait_df[r, level]])){
                trait_df[r,trait] <- trait.means[trait_df[r, level]]
              } else {
                # otherwise leave as NA
                trait_df[r,trait] <- NA
              }
            }
          }
        } else {
          # if it's not numeric, take the most common string
          trait.mode <- trait_df %>% 
            group_by_at(level) %>%
            dplyr::count(across(all_of(trait))) %>% 
            slice_max(1)
          
          for(r in 1:nrow(trait_df)){
            # run through by row and find the NAs
            if(is.na(trait_df[r,trait])){
              # if there are others of that taxa group, assign the most common trait they have
              if(!is.na(trait.mode[which(trait.mode[,level] == as.character(trait_df[r, level])),trait])){
                trait_df[r,trait] <- trait.mode[which(trait.mode[,level] == as.character(trait_df[r, level])),trait]
              } else {
                # otherwise leave as NA
                trait_df[r,trait] <- NA
              }
            }
          }
          
        }
        
      }
      
    }
    
    
  }
  
  return(trait_df)
  
}

# Check if column is numeric : if so, scale it, if not, don't
ScaleNumeric <- function(column){
  if(is.numeric(column)) 
    scaled_col <- scale(column)
  else 
    scaled_col <- column
  
  return(scaled_col)
}



GetSalixPhylo <- function(datadir,interaction_data){
  sal_phylo <- read.nexus(paste0(datadir, "salix_taxon.nex"))#"salix_phylogeny.nex"))
  
  # Get RSAL names for tip labels on sal_phylo
  sal_phylo_sp <- as.data.frame(sal_phylo$tip.label, stringsAsFactors = FALSE)
  names(sal_phylo_sp) <- "Salix"
  new_labels <- left_join(sal_phylo_sp,interaction_data[,c("Salix","RSAL")]) %>% unique
  new_labels[which(is.na(new_labels$RSAL)),"RSAL"] <- "SalExtra"
  sal_phylo$tip.label <- new_labels[,"RSAL"]
  
  return(sal_phylo)
}

# Convert an incidence (adjacency) matrix to an edge list
im2el <- function(im, extra.cols = NA, weighted = TRUE){
  M <- gather_(data = as.data.frame(im), 
               key_col = "consumer", 
               value_col = "IS", 
               gather_cols = colnames(im)[-which(colnames(im) == extra.cols)])
  M$IS <- as.numeric(M$IS)
  M[is.na(M$IS),"IS"] <- NA
  if(is.na(extra.cols)){
    M$resource <- rep(rownames(im), ncol(im))
  }else{
    M$resource <- rep(rownames(im), ncol(im)-length(extra.cols))  
  }
  if(!weighted){
    M$IS[which(M$IS > 0)] <- 1
  }
  
  #M <- M[M$IS >0,]
  return(M[,c(1,3,2)])
}


# Compute sum of log-likelihood for a binomial process
ll_fn <- function(L, p, rmvNaN = TRUE) {
  #p[which(p < 0.000000001)] <- 0.000000001
  ll = p*0
  ll[which(L==1)] = log(p[which(L==1)])
  ll[which(L==0)] = log(1-p[which(L==0)])
  if(rmvNaN) ll[which(is.infinite(ll) | is.nan(ll))] <- 0
  return(sum(ll))
}


tjur_D <- function(L,p) {
  D = mean(p[which(L==1)]) - mean(p[which(L==0)])
  return(D)
}

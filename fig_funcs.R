#### Functions to plot interaction matrices
#### Project: gallers
#### First edit: 2020.04.06
#### Latest edit: 2021.12.10
#### Author: Kate Wootton
#### R version 3.6.3

## Load libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# melt the probability matrix
melt_prb <- function(prb){
  prb_melt <- melt(matrix(prb, nrow = nrow(prb), ncol = ncol(prb), dimnames = list(rownames(prb),colnames(prb))), value.name = "Probability")
  colnames(prb_melt) <- c("resource","consumer", "Probability")
  return(prb_melt[,c(2,1,3)])
}

# Make a figure: Matrix of gallers (rows) and parasitoids or salix as columns.
# prb is the model prediction
# data is the alienData object
# res is the resource species (i.e. "Salix" or "Galler")
# con is the consumer species (i.e. "Galler" or "Parasitoid")
# group_res_by is a dataframe of the resource species names and the grouping variable 
compare_predictions_fig <- function(prb, data, res, cons, showpredict = TRUE, legpos = "right", plottitle = "", group_res_by=NA, group_con_by=NA, showcoexist = FALSE, flip = FALSE){
  
  # Convert the adjacency matrix into an edge list
  adjMat <- data$adjMat
  adjMat_melt <- im2el(adjMat, weighted = FALSE)
  
  # Do the same with the model output
  prb_melt <- melt_prb(prb)
  
  # Set the size of the points in the plot depending on number of data points
  size = case_when(
    sum(dim(adjMat)) > 135 ~ 1,
    sum(dim(adjMat)) < 100 ~ 3,
    TRUE                   ~ 2
  )
   
  # Set colors, gallers are always rows
  if(res == "Galler"){
    col <- "darkblue" 
    flip = FALSE
  }
  if(res == "Salix"){
    col <- "darkgreen"
    flip = TRUE
  }
  
  #who's on the x vs y axis? (This may get flipped later)
  yname = res
  xname = cons
  
  # Set up the grouping
  if(!is.na(group_res_by)){
    if(all(unique(adjMat_melt$resource) %in% group_res_by[,1])){
      colnames(group_res_by) <- c("resource","grouping_factor_res")
    }else{
      if(all(unique(adjMat_melt$resource) %in% group_res_by[,2])){
        colnames(group_res_by) <- c("grouping_factor_res","resource")
      } else{
        print("ERROR! Missing at least some resource species from the grouping factor")
        break
      }
    }
    prb_melt <- left_join(prb_melt, group_res_by, by = "resource") 
    adjMat_melt <- left_join(adjMat_melt, group_res_by, by = "resource")
    
    prb_melt$resource <- factor(prb_melt$resource, levels = unique(group_res_by$resource)[order(group_res_by$grouping_factor_res)])
    adjMat_melt$resource <- factor(adjMat_melt$resource, levels = unique(group_res_by$resource)[order(group_res_by$grouping_factor_res)])
  }else{
    prb_melt$resource <- factor(prb_melt$resource, levels = unique(adjMat_melt$resource)[order(adjMat_melt$resource)])
    adjMat_melt$resource <- factor(adjMat_melt$resource, levels = unique(adjMat_melt$resource)[order(adjMat_melt$resource)])
  }
  
  if(!is.na(group_con_by)){
    if(all(unique(adjMat_melt$consumer) %in% group_con_by[,1])){
      colnames(group_con_by) <- c("consumer","grouping_factor_con")
      }else{
      if(all(unique(adjMat_melt$consumer) %in% group_con_by[,2])){
        colnames(group_con_by) <- c("grouping_factor_con","consumer")
        prb_melt <- left_join(prb_melt, group_con_by, by = "consumer") 
        adjMat_melt <- left_join(adjMat_melt, group_con_by, by = "consumer") 
      } else{
        print("ERROR! Missing at least some consumer species from the grouping factor")
        break
      }
    }
    prb_melt <- left_join(prb_melt, group_con_by, by = "consumer") 
    adjMat_melt <- left_join(adjMat_melt, group_con_by, by = "consumer")
    prb_melt$consumer <- factor(prb_melt$consumer, levels = unique(group_con_by$consumer)[order(group_con_by$grouping_factor_con)])
    adjMat_melt$consumer <- factor(adjMat_melt$consumer, levels = unique(group_con_by$consumer)[order(group_con_by$grouping_factor_con)])
  }else{
    prb_melt$consumer <- factor(prb_melt$consumer, levels = unique(adjMat_melt$consumer[order(adjMat_melt$consumer)]))
    adjMat_melt$consumer <- factor(adjMat_melt$consumer, levels = unique(adjMat_melt$consumer)[order(adjMat_melt$consumer)])
    
  }
  
  if(showcoexist){
    prb_melt[which(is.na(adjMat_melt$IS)), "Probability"] <- 0
  }
  
  # If we want to show the model prediction as well as the interactions
  if(showpredict){ 
    # set fill as the value of the model prediction (Probability)
    p <- ggplot() + geom_tile(data = prb_melt, aes(prb_melt[,1],prb_melt[,2],fill = Probability)) #+
    size <- size - 0.4
  } else {
    # Otherwise set fill to zero
    p <- ggplot() + geom_tile(data = prb_melt, aes(prb_melt[,1],prb_melt[,2],fill = 0)) #+ 
  }
    
  # Set up the right grouping and x vs y axes
  if(!is.na(group_res_by) & !is.na(group_con_by) & flip == TRUE)
    fg <- facet_grid(rows = vars(grouping_factor_con), cols = vars(grouping_factor_res), scales = "free", space = "free") 
  if(!is.na(group_res_by) & !is.na(group_con_by) & flip == FALSE)
    fg <- facet_grid(rows = vars(grouping_factor_res), cols = vars(grouping_factor_con), scales = "free", space = "free")
  if(!is.na(group_res_by) & is.na(group_con_by) & flip == TRUE)
    fg <- facet_grid(cols = vars(grouping_factor_res), scales = "free", space = "free") 
  if(!is.na(group_res_by) & is.na(group_con_by) & flip == FALSE)
    fg <- facet_grid(rows = vars(grouping_factor_res), scales = "free", space = "free")
  if(is.na(group_res_by) & !is.na(group_con_by) & flip == TRUE)
    fg <- facet_grid(rows = vars(grouping_factor_con), scales = "free", space = "free") 
  if(is.na(group_res_by) & !is.na(group_con_by) & flip == FALSE)
    fg <- facet_grid(cols = vars(grouping_factor_con), scales = "free", space = "free")
  if(is.na(group_res_by) & is.na(group_con_by))
    fg <- NA
  

  pfig <-  p +
    scale_fill_gradient(low = "white", high = col, limits = c(0,1)) + 
    geom_point(data = adjMat_melt[which(adjMat_melt$IS == 1),], aes(adjMat_melt[which(adjMat_melt$IS == 1),1],adjMat_melt[which(adjMat_melt$IS == 1),2]), color = "orange", size = size ) +#0.8
    theme_classic() +
    ylab(yname) +
    xlab(xname)  +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title = element_text(size = 16),
          axis.line = element_blank(),
          legend.position = legpos,
          # Put a border around each facet panel
          panel.border=element_rect(color="black",size=0.5, fill=NA),
          # Change the size of the facet labels
          strip.text.y = element_text(size = 12, angle = 0),
          strip.background = element_rect(color = "grey", fill = "grey")) +
    ggtitle(plottitle)
  
  # Option to show combinations that never coexist as grey squares
  if(showcoexist){
    pfig <- pfig + 
      #geom_tile(data = adjMat_melt[which(is.na(adjMat_melt$IS)),], aes(adjMat_melt[which(is.na(adjMat_melt$IS)),1],adjMat_melt[which(is.na(adjMat_melt$IS)),2]), fill = "white") +
      geom_point(data = adjMat_melt[which(adjMat_melt$IS==0),], aes(adjMat_melt[which(adjMat_melt$IS==0),1],adjMat_melt[which(adjMat_melt$IS==0),2]), color = "grey", size = size, shape = 1) 
      #geom_point(data = adjMat_melt[which(is.na(adjMat_melt$IS)),], aes(adjMat_melt[which(is.na(adjMat_melt$IS)),1],adjMat_melt[which(is.na(adjMat_melt$IS)),2]), color = "darkgrey", size = size, shape = 3) 
    #pfig <- pfig + geom_tile(data = adjMat_melt[which(is.na(adjMat_melt$IS)),], aes(adjMat_melt[which(is.na(adjMat_melt$IS)),1],adjMat_melt[which(is.na(adjMat_melt$IS)),2]), fill = "grey") 
  }
  
  if(flip){
    pfig <- pfig +
      coord_flip() 
  }
  
  pfig <- pfig + fg 
  
  return(pfig)
  
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



#### Functions to compare model predictions
#### Project: gallers
#### First edit: 2020.04.04
#### Latest edit: 2021.12.10
#### Author: Kate Wootton
#### R version 3.6.2


## Load libraries
library(psych)
library(reshape2)

## Melt each matrix of predictions then join all data frames
make_prediction_df <- function(modelList){
  
  meltedList <- lapply(seq_along(modelList), function(x) melt(as.data.frame(as.table(modelList[[x]])), value.name = names(modelList)[[x]]))
  modelPreds <- meltedList %>% purrr::reduce(left_join)
  
  return(modelPreds)
}


# Get correlation between models, with option to plot
get_prediction_correlations <- function(modelList, plot = TRUE){
  
  # Get the predictions for each model into one dataframe
  modelPreds <- make_prediction_df(modelList)
  
  # Rename the columns to be prettier for printing
  colnames(modelPreds) <- recode(colnames(modelPreds), 
                                 KNN_phylo = "KNN (phylo)", 
                                 KNN_traits = "KNN (traits)", 
                                 GLM = "GLM", 
                                 RF_traits = "RF (traits)", 
                                 fourC = "Fourth Corner", 
                                 RF_phylo = "RF (phylo)")
  # Reorder columns
  modelPreds <- modelPreds %>% select(c(Var1:KNN, starts_with("KNN"),starts_with("RF"), "IMC","GLM",`Fourth Corner`))
  
  # Get the correlations
  modelCors <- cor(modelPreds[,4:ncol(modelPreds)])
  
  if(plot)  pairs.panels(modelPreds[,4:ncol(modelPreds)], cex = 1, cex.labels = 1)
  return(modelCors)
}

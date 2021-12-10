#### Functions needed to sample interactions and create networks based on model predictions
#### Project: gallers
#### First edit: 2020.04.06
#### Latest edit: 2021.12.10
#### Author: Kate Wootton
#### R version 3.6.3

# Using the model prediction, sample a food web
sample_1_network <- function(a = 1, 
                             prediction){
  v_prediction <- as.vector(prediction)
  v_sampled <- ifelse(runif(length(v_prediction)) < v_prediction, 1, 0)
  nw_sampled <- matrix(v_sampled, 
                       ncol = ncol(prediction), 
                       nrow = nrow(prediction), 
                       dimnames = list(rownames(prediction),
                                       colnames(prediction)))

  return(nw_sampled)
}

# Sample n networks based on the model prediction
sample_many_networks <- function(prediction, n, return_form = "list"){
  sampled_nws <- vector("list", n) %>% 
    lapply(sample_1_network, prediction = prediction)

  return(sampled_nws)
}

# Calculate modularity
calc_modularity <- function(nw){
  mods <- computeModules(nw)
  return(mods@likelihood)
}

# Calculate how many interactions were correctly predicted
# "nw" is the sampled network, "observed_nw" is the empirical network
calc_correct1s <- function(nw, 
                           observed_nw){
  obs_long <- pivot_longer(as_tibble(observed_nw), 
                           cols = 1:ncol(observed_nw))
  smpl_long <- pivot_longer(as_tibble(nw, .name_repair = "minimal"), 
                            cols = 1:ncol(observed_nw))
  correct1s <- sum(smpl_long[which(obs_long$value == 1),"value"])/sum(obs_long$value, na.rm = TRUE)
  
  return(correct1s)
}

# Calculate how many non-interactions were correctly predicted
# "nw" is the sampled network, "observed_nw" is the empirical network
calc_correct0s <- function(nw, 
                           observed_nw){
  obs_long <- pivot_longer(as_tibble(observed_nw), 
                           cols = 1:ncol(observed_nw))
  smpl_long <- pivot_longer(as_tibble(nw, .name_repair = "minimal"), 
                            cols = 1:ncol(observed_nw))
  correct0s <- (length(which(obs_long$value == 0))-sum(smpl_long[which(obs_long$value == 0),"value"]))/length(which(obs_long$value==0))
  
  return(correct0s)
}

# Calculate metrics for each sampled network
calc_metrics_sampled_nw <- function(sampled_nws, observed_nw){
  metrics <- tibble(.rows=length(sampled_nws))

  metrics$Connectance <- sapply(sampled_nws, function(x) sum(x)/(nrow(x)*ncol(x)))
  metrics$Modularity <- sapply(sampled_nws, calc_modularity)
  metrics$Nestedness <- sapply(sampled_nws, function(x) nested(x, method = "NODF2"))
  metrics$DDskew_res <- sapply(sampled_nws, function(x) skewness(rowSums(x, na.rm = TRUE)))
  metrics$DDskew_cons <- sapply(sampled_nws, function(x) skewness(colSums(x, na.rm = TRUE)))
  metrics$DDspread_res <- sapply(sampled_nws, function(x) sd(rowSums(x, na.rm = TRUE)))
  metrics$DDspread_cons <- sapply(sampled_nws, function(x) sd(colSums(x, na.rm = TRUE)))
  metrics$correct1s <- sapply(sampled_nws, function(x) calc_correct1s(nw = x, observed_nw = observed_nw))
  metrics$correct0s <- sapply(sampled_nws, function(x) calc_correct0s(nw = x, observed_nw = observed_nw))
 
  return(metrics)
}

# Calculate the same metrics for the real network
calc_metrics_observed_nw <- function(nw){
  metrics <- tibble(.rows = 1)

  metrics$Connectance <- sum(nw)/(nrow(nw)*ncol(nw))
  mods <- computeModules(nw)
  metrics$Modularity <- mods@likelihood
  metrics$Nestedness <- nested(nw, method = "NODF2")
  metrics$DDskew_res <- skewness(rowSums(nw, na.rm = TRUE))
  metrics$DDskew_cons <- skewness(colSums(nw, na.rm = TRUE))
  metrics$DDspread_res <- sd(rowSums(nw, na.rm = TRUE))
  metrics$DDspread_cons <- sd(colSums(nw, na.rm = TRUE))
  # These two are 1 because we're assuming the data is correct (i.e. all interactions observed exist)
  metrics$correct1s <- 1
  metrics$correct0s <- 1
  
  return(metrics)
}

# Set combinations that never co-occur (and therefore are set to NA) to 0
set_non_cooccur_to_0 <- function(nw, obs_nw){
  nw0 <- nw
  nw0[which(is.na(obs_nw))] <- 0
  return(nw0)
}

  
# gallers-paper-code
Code used in the paper "Differential imprints of trait-matching in a tritrophic network: complementary insights from complementary methods"


## master_code.R
Run through this to run all analyses. Analyses/figures are labeled with the figure or table they correspond to.  Loads the following files which contain necessary functions.

## network_analyses_funcs.R
Contains functions needed to run network analyses (e.g. for Random Forest, GLM etc)
  - `get_RF_traits_formula` outputs the formula for traits when using Random Forest
  - `get_GLM_traits_formula` outputs the formula for traits when using GLM 
  - `get_RF_phylo_traits` gets the phylogeny "traits"
  - `get_RF_phylo_formula` outputs the formula for using phylogeny with Random Forest

## model_correlations.R
Contains functions used to compare model predictions.
  - `get_prediction_correlations` runs correlations on all model outputs, with option to print.
  - `make_prediction_df` makes a data frame of predictions from all models 
 
 ## fig_funcs2.R
 Contains functions to make interaction matrices figures showing model predictions.
  - `compare_predictions_fig` makes the figure
  - `melt_prb` turns the model prediction matrix into an edge list
  - `im2el` converts an incidence (adjacency) matrix into an edge list
  - `make_sp_names_table` outputs a latex table with the species names in the order they appear in the figure. Gallers are grouped by gall type
  
## plot_intdist_phylodist.R
Contains functions to calculate and plot the correlation of a community's phylogenetic, interaction, and trait dissimilarities/distances.
  - `calc_dist_ints` calculates distances between species based on shared interactions
  - `calc_dist_traits` calculates distances between species based on traits
  - `compare_distances` calculate (and plot) correlation between species' interaction, trait, and phylogenetic distances
  
## run_trait_combinations.R
Contains functions to run each model that uses traits (KNN, GLM, 4th Corner, RF) with all combination of n traits
  - `run_trait_combinations` runs each model with each combination of n traits
  - `scaleFUN` scales numbers to 2 decimal places (for plotting)
  - `plot_ntrait_fits` plots the fits of each model with different traits as a box plot faceted by model

  ## trim_traits.R
  Contains functions to remove correlated or un-useful traits
   - `trimTraits` removes traits that are correlated or just not that useful. See supplementary material for description of choices.
   - `removeTraits` is called by `trimTraits` and actually removes the traits.

#### Script to run analyses for the paper:
#### "Differential imprints of trait-matching in a tritrophic network: complementary insights from complementary methods"
#### Wootton et al 2021
#### Project: gallers
#### First edit: 2020.04.06
#### Latest edit: 2021.12.10
#### Author: Kate Wootton
#### R version 3.6.3


#####################################
## Set up 
#####################################

# load packages
library(ggpubr)        
library(bipartite)
library(e1071)
library(stargazer)
library(picante)
library(randomForest)

# Load Alien
#install.packages("remotes")
#remotes::install_github("TheoreticalEcosystemEcology/alien")
library(alien)


# Load scripts
setwd("/home/kate/Documents/PhDprojects/galler-traits")
source("submission_202105/code/network_analyses_funcs.R")
source("submission_202105/code/model_correlations.R")
source("submission_202105/code/fig_funcs2.R")
source("submission_202105/code/plot_intdist_phylodist.R")
source("submission_202105/code/sample_webs_funcs.R")
source("submission_202105/code/run_trait_combinations.R")
source("code/R.code/trim_traits.R")
source("code/R.code/train_test.R")


#####################################################
## Get data
#####################################################

# Load data
data(salixGal)
data(galPara)

# Clean traits to remove correlated or unuseful traits
salixGal <- trimTraits(salixGal)
galPara  <- trimTraits(galPara)

# And extract interaction matrices
gp_mat <- galPara$adjMat
sg_mat <- salixGal$adjMat


#################################################
## Run analyses
#################################################

# Set up lists to hold model output
models_sg <- list()
models_gp <- list()

#### Salix-galler ####

## KNN methods (will output warnings because not all species have 3 neighbours)
models_sg[["KNN"]] <- fitKNN(salixGal,
                           distFrom = "jaccard",
                           distTo = "jaccard",
                           nNeig = 3)

models_sg[["KNN_traits"]] <- fitKNN(salixGal,
                                distFrom = "jaccard",
                                distTo = "jaccard",
                                distTraitFrom = "euclidean",
                                distTraitTo = "euclidean",
                                nNeig = 3,
                                weight = 1)

models_sg[["KNN_phylo"]] <- fitKNN(salixGal,
                                distFrom = "jaccard",
                                distTo = "jaccard",
                                nNeig = 3,
                                weight = 1,
                                phylo=TRUE)

## RF methods
models_sg[["RF_traits"]] <- fitRF(salixGal,
                         formula = get_RF_traits_formula(salixGal),
                         ntree = 2000,
                         nodesize = 1)

salixGalPhylo <- get_RF_phylo_traits(salixGal)
models_sg[["RF_phylo"]] <- fitRF(salixGalPhylo,
                                 formula = get_RF_phylo_formula(salixGalPhylo),
                                 ntree = 2000,
                                 nodesize = 1)

## GLM & fourth-corner
# Algorithms don't converge with all traits, so lets use the most informative ones
gall_traits <- c("GALLTYPE", 
                 "BODYLENGTH.GAL", 
                 "PHENOLOGY.GAL", 
                 "OVIPOS.STRATEGY",
                 "GALL.WALL")

sal_traits  = c("TREE.VOLUME",
                "LEAF.THICKNESS",
                "LEAF.HAIRINESS",
                "GLUCOSIDES")

models_sg[["GLM"]] <- fitGLM(polyTrait(salixGal),
                             formula = get_GLM_traits_formula(salixGal,traits2use = c(gall_traits,sal_traits),which_model = "GLM"),
                             family = binomial(link = "logit"))

models_sg[["fourC"]] <- fitGLM(polyTrait(salixGal),
                               formula = get_GLM_traits_formula(salixGal, traits2use = c(gall_traits,sal_traits), which_model = "4C"),
                               family = binomial(link = "logit"))

## IMC (this takes a long time to run)
#models_sg[["IMC"]] <- fitIMC(salixGal,
#                             d = 1)
load("salixGalIMC.Rdata")
models_sg[["IMC"]] <- salixGalIMC

#### Galler-parasitoid ####

## KNN methods (will output warnings because not all species have 3 neighbours)
models_gp[["KNN"]] <- fitKNN(galPara,
                             distFrom = "jaccard",
                             distTo = "jaccard",
                             nNeig = 3)

models_gp[["KNN_traits"]] <- fitKNN(galPara,
                                    distFrom = "jaccard",
                                    distTo = "jaccard",
                                    distTraitFrom = "euclidean",
                                    distTraitTo = "euclidean",
                                    nNeig = 3,
                                    weight = 1)

models_gp[["KNN_phylo"]] <- fitKNN(galPara,
                                   distFrom = "jaccard",
                                   distTo = "jaccard",
                                   nNeig = 3,
                                   weight = 1,
                                   phylo=TRUE)

## RF methods
models_gp[["RF_traits"]] <- fitRF(galPara,
                                  formula = get_RF_traits_formula(galPara),
                                  ntree = 2000,
                                  nodesize = 1)

galParaPhylo <- get_RF_phylo_traits(galPara)
models_gp[["RF_phylo"]] <- fitRF(galParaPhylo,
                                 formula = get_RF_phylo_formula(galParaPhylo),
                                 ntree = 2000,
                                 nodesize = 1)


## GLM & fourth-corner

# Algorithms don't converge with all traits, so lets use the most informative ones
gall_traits <- c("GALLTYPE", 
                "BODYLENGTH.GAL", 
                "PHENOLOGY.GAL", 
                "OVIPOS.STRATEGY",
                "GALL.WALL")

para_traits <- c("P.I",  
                "OVIPOS.LNTH", 
                "ATTACK.STAGE",
                "PHENOLOGY.PAR", 
                "BODYLENGTH.PAR", 
                "ENDO.ECTO")

models_gp[["GLM"]] <- fitGLM(polyTrait(galPara),
                           formula = get_GLM_traits_formula(galPara,traits2use = c(gall_traits,para_traits),which_model = "GLM"),
                           family = binomial(link = "logit"))
models_gp[["fourC"]] <- fitGLM(polyTrait(galPara),
                               formula = get_GLM_traits_formula(galPara, traits2use = c(gall_traits,para_traits), which_model = "4C"),
                               family = binomial(link = "logit"))


## IMC (this takes a long time to run)
#models_gp[["IMC"]] <- fitIMC(galPara,
#                     d = 1)
load("galParaIMC.Rdata")
models_gp[["IMC"]] <- galParaIMC


##################################
## Check model fits
##################################

## Model fit (Table 1 of manuscript)

## Salix galler
ll_fits_sg <- list()
for(m in names(models_sg)){
  ll_fits_sg[[m]] <- logLik(models_sg[[m]], error = 0.001)
}

D_fits_sg <- list()
for(m in names(models_sg)){
  D_fits_sg[[m]] <- tjur(models_sg[[m]])
}

fits_sg <- rbind(unlist(ll_fits_sg, use.names = TRUE), unlist(D_fits_sg, use.names = TRUE)) 


## Galler parasitoid
ll_fits_gp <- list()
for(m in names(models_gp)){
  ll_fits_gp[[m]] <- logLik(models_gp[[m]], error = 0.001)
}

D_fits_gp <- list()
for(m in names(models_gp)){
  D_fits_gp[[m]] <- tjur(models_gp[[m]])
}

fits_gp <- rbind(unlist(ll_fits_gp, use.names = TRUE), unlist(D_fits_gp, use.names = TRUE)) 

####################################################################
## Correlations between model predictions (Tables 2 & 3 and Fig S2) 
####################################################################

png(filename = "submission_202105/figs/model_correlations.png", width = 600, height = 500)
  modelCors_gp <- get_prediction_correlations(models_gp, plot = TRUE)
dev.off()

png(filename = "submission_202105/figs/modelCors_sg.png", width = 600, height = 500)
  modelCors_sg <- get_prediction_correlations(models_sg, plot = TRUE)
dev.off()

# Create latex tables
stargazer(modelCors_gp, title="Correlation between each model's predictions when using all traits (or phylogeny) for the galler-parasitoid data")
stargazer(modelCors_sg, title="Correlation between each model's predictions when using all traits (or phylogeny) for the \textit{Salix}-galler data")

###################################################################
## Make interaction matrices figures (Fig 2 & Figs S3-S11)
###################################################################

# What traits should we group the resource/consumer by?
group_res_by = data.frame(RGALLER=rownames(galPara$traitFrom),
                          grouping_factor = plyr::mapvalues(galPara$traitFrom$GALLTYPE, from = levels(as.factor(galPara$traitFrom$GALLTYPE)), to = c("BG","LBBG","LBSG","LF","LMBG","Pt","LMPG","LR","LR", "Pt", "SG")))
group_con_by = data.frame(RGALLER=rownames(salixGal$traitTo),
                          grouping_factor = plyr::mapvalues(salixGal$traitTo$GALLTYPE, from = levels(as.factor(salixGal$traitTo$GALLTYPE)), to = c("BG","LBBG","LBSG","LF","LMBG","Pt","LMPG","LR","LR", "Pt", "SG")))

# No prediction (only interaction matrix)
pdf("submission_202105/figs/gp_interaction_matrix_showcoexist.pdf", height=8, width=9)
  compare_predictions_fig(models_gp[["KNN"]], galPara, "Galler", "Parasitoid", group_res_by = group_res_by,showpredict = FALSE, legpos = "none",showcoexist = TRUE)
dev.off()

pdf("submission_202105/figs/sg_interaction_matrix_showcoexist.pdf", height=8, width=9)
  compare_predictions_fig(models_sg[["KNN"]], salixGal, "Salix", "Galler", group_con_by = group_con_by, showpredict = FALSE, legpos = "none",showcoexist = TRUE)
dev.off()

model = names(models_gp)[9]
for(model in names(models_gp)){
  pdf(paste0("submission_202105/figs/",model,"_gp_showcoexist.pdf"), height=8, width=9)
    compare_predictions_fig(models_gp[[model]], galPara, "Galler", "Parasitoid", group_res_by = group_res_by,legpos = "none", showcoexist = TRUE)
  dev.off()
  
  pdf(paste0("submission_202105/figs/",model,"_sg_showcoexist.pdf"), height=8, width=9)
    compare_predictions_fig(models_sg[[model]], salixGal, "Salix", "Galler",  group_con_by = group_con_by, legpos = "none",showcoexist = TRUE)
  dev.off()
}

########################################################################
## Compare interaction distance vs phylogenetic distance (figs S12-S14)
########################################################################

png(filename = "submission_202105/figs/galler_phylo_int_traits_correlations.png", width = 300, height = 300)
  distCors_gp <- compare_distances(alienData = galPara, troph = "from", toPlot = TRUE)
dev.off()

png(filename = "submission_202105/figs/par_phylo_ints_traits_correlations.png", width = 300, height = 300)
  distCors_gp <- compare_distances(alienData = galPara, troph = "to", toPlot = TRUE)
dev.off()

png(filename = "submission_202105/figs/Gal_phy_int_trait_cor_wSal.png", width = 300, height = 300)
  distCors_sg <- compare_distances(alienData = salixGal, troph = "to", toPlot = TRUE)
dev.off()

png(filename = "submission_202105/figs/sal_phy_int_trait_cors.png", width = 300, height = 300)
  distCors_sg <- compare_distances(alienData = salixGal, troph = "from", toPlot = TRUE)
dev.off()

###########################################################
## Compare trait importance (from Random forest) (fig S17)
###########################################################

jpeg(filename = "submission_202105/figs/rf_trait_importance.jpeg", width = 300, height = 400)
  varImpPlot(attributes(models_gp[["RF_traits"]])$model, main =NA)
dev.off()

png(filename = "submission_202105/figs/rf_trait_importance_sg.png", width = 300, height = 400)
  varImpPlot(attributes(models_sg[["RF_traits"]])$model, main = NA)
dev.off()




##########################################################################
## Refitting the models with different combinations of traits (fig 3)
##########################################################################

#(Takes a long time to run)
ntraits = 4
#fits_sg <- run_trait_combinations(salixGal, ntraits, to_save = TRUE)
#fits_gp <- run_trait_combinations(galPara, ntraits, to_save = TRUE)

pdf(paste0("submission_202105/figs/traits_scores_",ntraits,"traits_sg_boxplot.pdf"), height = 6, width = 12)
  plot_ntrait_fits(salixGal, fits_sg, ntraits)
dev.off()

pdf(paste0("submission_202105/figs/traits_scores_",ntraits,"traits_gp_boxplot.pdf"), height = 6, width = 12)
  plot_ntrait_fits(galPara, fits_gp, ntraits = 4)
dev.off()


##########################################################################
## Comparing performance of GLM and RF (fig S15)
##########################################################################

pdf("submission_202105/figs/GLM_vs_RF_col_ntraits_sg.pdf", height=4, width=5)
ggplot(fits_sg) + 
  geom_point(aes(x = GLM, y = RF, col = total.ntraits), size = 1) +
  geom_abline(slope = 1, intercept = 0) + 
  scale_color_gradient(low = "yellow", high = "slateblue", name = "Total number of trait values") +
  xlab("Model performance of GLM") + 
  ylab("Model performance of RF") +
  theme_pubr() 
dev.off()


pdf("submission_202105/figs/GLM_vs_RF_col_ntraits_gp.pdf", height=4, width=5)
  ggplot(fits_gp) + 
    geom_point(aes(x = GLM, y = RF, col = total.ntraits), size = 1) +
    geom_abline(slope = 1, intercept = 0) + 
    scale_color_gradient(low = "yellow", high = "slateblue", name = "Total number of trait values") +
    xlab("Model performance of GLM") + 
    ylab("Model performance of RF") +
  theme_pubr()
dev.off()

##########################################################################
## Performance of RF with more trait values (Fig S16)
##########################################################################


pdf("submission_202105/figs/ntraits_vs_tjurD_sg.pdf", height=4, width=5)
  plot(fits_sg$total.ntraits,fits_sg$RF, xlab = "Number of unique trait values", ylab = "Model fit (Tjur's D)", col = rgb(red=0,green = 0,blue = 0,alpha=0.3),pch=16)
dev.off()

pdf("submission_202105/figs/ntraits_vs_tjurD_gp.pdf", height=4, width=5)
  plot(fits_gp$total.ntraits,fits_gp$RF, xlab = "Number of unique trait values", ylab = "Model fit (Tjur's D)", col = rgb(red=0,green = 0,blue = 0,alpha=0.3),pch=16)
dev.off()



#########################################################
## Filling holes in the web and calculating metrics 
#########################################################

## Can load using following commands, or re-run using lines 364-400
#load("data/alien/sampled_metrics_sg_0.Rdata")
#load("data/alien/sampled_metrics_sg.Rdata")
#load("data/alien/sampled_webs_from_predictions_sg_0.Rdata")
#load("data/alien/sampled_webs_from_predictions_sg.Rdata")

#load("data/alien/sampled_metrics_gp_0.Rdata")
#load("data/alien/sampled_metrics_gp.Rdata")
#load("data/alien/sampled_webs_from_predictions_gp_0.Rdata")
#load("data/alien/sampled_webs_from_predictions_gp.Rdata")


#### galler-parasitoid ####
sampled_from_predictions_gp <- list()
sampled_metrics_gp <- list()
sampled_from_predictions_gp_0 <- list()
sampled_metrics_gp_0 <- list()

for(m in names(models_gp)){
  ## Sample to construct networks from the models. Do n replicates
  sampled_from_predictions_gp[[m]] <- sample_many_networks(models_gp[[m]], n = 100)
  
  # Set non-co occurring species' interactions to zero
  sampled_from_predictions_gp_0[[m]] <- lapply(sampled_from_predictions_gp[[m]], set_non_cooccur_to_0, obs_nw = galPara$adjMat) #function(x) x[which(is.na(gp_mat))] = 0)
  
  # Calculate metrics (connectance etc) from the sampled networks 
  sampled_metrics_gp_0[[m]] <- calc_metrics_sampled_nw(sampled_from_predictions_gp_0[[m]], observed_nw =galPara$adjMat) %>% mutate(model = m)
  
}


save(sampled_metrics_gp_0, file = "submission_202105/data/sampled_metrics_gp_0.Rdata")
##save(sampled_metrics_gp, file = "data/alien/sampled_metrics_gp.Rdata")
save(sampled_from_predictions_gp_0, file = "submission_202105/data/sampled_webs_from_predictions_gp_0.Rdata")
save(sampled_from_predictions_gp, file = "submission_202105/data/sampled_webs_from_predictions_gp.Rdata")

#### salix-galler ####
sampled_from_predictions_sg <- list()
sampled_metrics_sg <- list()
sampled_from_predictions_sg_0 <- list()
sampled_metrics_sg_0 <- list()

for(m in names(models_sg)){
  ## Sample to construct networks from the models. Do n replicates
  #sampled_from_predictions_sg[[m]] <- sample_many_networks(models_sg[[m]], n = 100)

  # Set non-co occurring species' interactions to zero
  #sampled_from_predictions_sg_0[[m]] <- lapply(sampled_from_predictions_sg[[m]], set_non_cooccur_to_0, obs_nw = salixGal$adjMat)
  
  # Set non-co occurring species' interactions to zero
  sampled_metrics_sg_0[[m]] <- calc_metrics_sampled_nw(sampled_from_predictions_sg_0[[m]], observed_nw=salixGal$adjMat) %>% mutate(model = m)
}

save(sampled_metrics_sg_0, file = "submission_202105/data/sampled_metrics_sg_0.Rdata")
#save(sampled_metrics_sg, file = "data/alien/sampled_metrics_sg.Rdata")
save(sampled_from_predictions_sg_0, file = "submission_202105/data/sampled_webs_from_predictions_sg_0.Rdata")
save(sampled_from_predictions_sg, file = "submission_202105/data/sampled_webs_from_predictions_sg.Rdata")

## Need to run from here even if loading data

# Create an adjacency matrix with the NA values (no cooccurrence) set to 0
sg_mat_0 <- salixGal$adjMat
sg_mat_0[which(is.na(sg_mat_0))] <- 0
gp_mat_0 <- galPara$adjMat
gp_mat_0[which(is.na(gp_mat_0))] <- 0


# Calculate metrics on the observed network
observed_metrics_sg <- cbind(calc_metrics_observed_nw(sg_mat_0), model = "observed")
observed_metrics_gp <- cbind(calc_metrics_observed_nw(gp_mat_0), model = "observed")

# Reformat the data for figures
pivot_sampled_metrics_sg_0 <- do.call("rbind",sampled_metrics_sg_0) %>% #rbind(observed_metrics_gp[rep(1,n), ]) %>%
  pivot_longer(cols = colnames(.)[which(colnames(.) != "model")], values_to = "Value") %>%
  mutate(name = dplyr::recode(name, correct0s = "No interaction", correct1s = "Interaction", DDskew_cons = "Skew (G)", DDskew_res = "Skew (S)", DDspread_cons = "Spread (G)", DDspread_res = "Spread (S)")) %>%
  mutate(name = factor(name, ordered = TRUE, levels = c("Interaction", "No interaction", "Connectance", "Modularity", "Nestedness", "Skew (S)","Spread (S)", "Skew (G)", "Spread (G)"))) %>%
  mutate(Model = recode(model, GLM = "GLM", RF_traits = "RF", fourC = "4-corner", KNN_phylo = "KNN (phylo)", KNN_traits = "KNN (traits)", RF_phylo = "RF (phylo)")) 

pivot_observed_metrics_sg <- observed_metrics_sg %>% 
  pivot_longer(cols = colnames(.)[which(colnames(.) != "model")], values_to = "Value") %>%
  mutate(name = dplyr::recode(name, correct0s = "No interaction", correct1s = "Interaction", DDskew_cons = "Skew (G)", DDskew_res = "Skew (S)", DDspread_cons = "Spread (G)", DDspread_res = "Spread (S)")) %>%
  mutate(name = factor(name, ordered = TRUE, levels = c("Interaction", "No interaction", "Connectance", "Modularity", "Nestedness", "Skew (S)","Spread (S)", "Skew (G)", "Spread (G)"))) %>%
  mutate(Model = recode(model, GLM = "GLM", RF_traits = "RF", fourC = "4-corner", KNN_phylo = "KNN (phylo)", KNN_traits = "KNN (traits)", RF_phylo = "RF (phylo)")) 

pivot_sampled_metrics_gp_0 <- do.call("rbind",sampled_metrics_gp_0) %>% #rbind(observed_metrics_gp[rep(1,n), ]) %>%
  pivot_longer(cols = colnames(.)[which(colnames(.) != "model")], values_to = "Value") %>%
  mutate(name = dplyr::recode(name, correct0s = "No interaction", correct1s = "Interaction", DDskew_cons = "Skew (P)", DDskew_res = "Skew (G)", DDspread_cons = "Spread (P)", DDspread_res = "Spread (G)")) %>%
  mutate(name = factor(name, ordered = TRUE, levels = c("Interaction", "No interaction", "Connectance", "Modularity", "Nestedness", "Skew (G)","Spread (G)", "Skew (P)", "Spread (P)"))) %>%
  mutate(Model = recode(model, GLM = "GLM", RF_traits = "RF", fourC = "4-corner", KNN_phylo = "KNN (phylo)", KNN_traits = "KNN (traits)", RF_phylo = "RF (phylo)"))

pivot_observed_metrics_gp <- observed_metrics_gp %>% 
  pivot_longer(cols = colnames(.)[which(colnames(.) != "model")], values_to = "Value") %>%
  mutate(Model = recode(model, GLM = "GLM", RF_traits = "RF", fourC = "4-corner", KNN_phylo = "KNN (phylo)", KNN_traits = "KNN (traits)", RF_phylo = "RF (phylo)")) %>%
  mutate(name = dplyr::recode(name, correct0s = "No interaction", correct1s = "Interaction", DDskew_cons = "Skew (P)", DDskew_res = "Skew (G)", DDspread_cons = "Spread (P)", DDspread_res = "Spread (G)")) %>%
  mutate(name = factor(name, ordered = TRUE, levels = c("Interaction", "No interaction", "Connectance", "Modularity", "Nestedness", "Skew (G)","Spread (G)", "Skew (P)", "Spread (P)")))



## Make a density plot of the metrics from the sampled networks (fig 4)
p.gp0 <- ggdensity(pivot_sampled_metrics_gp_0, x = "Value", fill = "Model", ylab = "Density") +
  theme(legend.position = "none")
pdf("submission_202105/figs/sampled_gp0.pdf", height = 5.5, width = 9)
  facet(p.gp0, facet.by = c("name"), scales = "free") + geom_vline(aes(xintercept = Value), data = pivot_observed_metrics_gp, lty = "dashed")
dev.off()

p.sg0 <- ggdensity(pivot_sampled_metrics_sg_0, x = "Value", fill = "Model", ylab = "Density")  +
  theme(legend.position = "bottom")
pdf("submission_202105/figs/sampled_sg0.pdf", height = 6.3, width = 9)
  facet(p.sg0, facet.by = c("name"), scales = "free") + geom_vline(aes(xintercept = Value), data = pivot_observed_metrics_sg, lty = "dashed")
dev.off()




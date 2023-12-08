#### Script to run analyses for paper
#### "Layer-specific imprints of traits within a plant-herbivore-predator network - complementary insights from complementary methods"
#### Wootton et al 2024
#### Ecography
#### First edit: 2020.04.06
#### Latest edit: 2023.12.10
#### Author: Kate Wootton
#### R version 4.2.2


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

# Install Alien

install.packages("remotes")
remotes::install_github("TheoreticalEcosystemEcology/alien")
library(alien)


# Load scripts
source("network_analyses_funcs.R")
source("model_correlations.R")
source("fig_funcs2.R")
source("plot_intdist_phylodist.R")
source("run_trait_combinations.R")
source("trim_traits.R")

# Load data
data(salixGal)
data(galPara)

# Clean traits to remove correlated or un-useful traits
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

# IMC (this takes a long time to run, as in hours to days)
models_sg[["IMC"]] <- fitIMC(salixGal, d = 1)

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


## IMC (this takes a long time to run, as in hours to days)
models_gp[["IMC"]] <- fitIMC(galPara, d = 1)


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

AIC <- data.frame(matrix(data = 0, nrow = 3, ncol =2), row.names = c("GLM","fourC","IMC"))
colnames(AIC) <- c("sg","gp")
AIC["GLM","sg"] <- attributes(models_sg[["GLM"]])$model$aic
AIC["GLM","gp"] <- attributes(models_gp[["GLM"]])$model$aic
AIC["fourC","sg"] <- attributes(models_sg[["fourC"]])$model$aic
AIC["fourC","gp"] <- attributes(models_gp[["fourC"]])$model$aic
AIC["IMC","sg"] <- 2*8 - 2*fits_sg[1,"IMC"]
AIC["IMC","gp"] <- 2*8 - 2*fits_gp[1,"IMC"]

####################################################################
## Correlations between model predictions (Tables 2 & 3 and Fig S3) 
####################################################################

modelCors_gp <- get_prediction_correlations(models_gp, plot = TRUE)
modelCors_sg <- get_prediction_correlations(models_sg, plot = TRUE)


# Create latex tables
stargazer(modelCors_gp, title="Correlation between each model's predictions when using all traits (or phylogeny) for the galler-parasitoid data")
stargazer(modelCors_sg, title="Correlation between each model's predictions when using all traits (or phylogeny) for the \textit{Salix}-galler data")

###################################################################
## Make interaction matrices figures (Fig 2 & Figs S4-S12)
###################################################################

# What traits should we group the resource/consumer by?
group_res_by = data.frame(RGALLER=rownames(galPara$traitFrom),
                          grouping_factor = plyr::mapvalues(galPara$traitFrom$GALLTYPE, from = levels(as.factor(galPara$traitFrom$GALLTYPE)), to = c("BG","LBBG","LBSG","LF","LMBG","Pt","LMPG","LR","LR", "Pt", "SG")))
group_con_by = data.frame(RGALLER=rownames(salixGal$traitTo),
                          grouping_factor = plyr::mapvalues(salixGal$traitTo$GALLTYPE, from = levels(as.factor(salixGal$traitTo$GALLTYPE)), to = c("BG","LBBG","LBSG","LF","LMBG","Pt","LMPG","LR","LR", "Pt", "SG")))

# No prediction (only interaction matrix)
compare_predictions_fig(models_gp[["KNN"]], galPara, "Galler", "Parasitoid", group_res_by = group_res_by,showpredict = FALSE, legpos = "none",showcoexist = TRUE)
compare_predictions_fig(models_sg[["KNN"]], salixGal, "Salix", "Galler", group_con_by = group_con_by, showpredict = FALSE, legpos = "none",showcoexist = TRUE)

model = names(models_gp)[9]
for(model in names(models_gp)){
   print(compare_predictions_fig(models_gp[[model]], galPara, "Galler", "Parasitoid", group_res_by = group_res_by,legpos = "none", showcoexist = TRUE))
   print(compare_predictions_fig(models_sg[[model]], salixGal, "Salix", "Galler",  group_con_by = group_con_by, legpos = "none",showcoexist = TRUE))
}

################################################################################
## Make boxplot of predicted values for interactions vs no ints for each model (fig 3)
################################################################################

# Convert the adjacency matrix into an edge list
probs_gp <- im2el(gp_mat, weighted = FALSE) %>%
  mutate(web = "Galler-parasitoid")
probs_sg <- im2el(sg_mat, weighted = FALSE) %>%
  mutate(web = "Salix-galler")

for(m in names(models_gp)){
  gp <- melt_prb(models_gp[[m]])
  colnames(gp)[3] <- m
  probs_gp <- left_join(probs_gp, gp)
  
  sg <- melt_prb(models_sg[[m]])
  colnames(sg)[3] <- m
  probs_sg <- left_join(probs_sg, sg)
}

probs <- rbind(probs_gp, probs_sg)


probs_long <- probs %>% 
  pivot_longer(cols = KNN:IMC, names_to = "model", values_to = "Model-predicted probability") %>%
  filter(!is.na(IS)) %>%
  mutate(model = factor(model, 
                        levels = c("KNN","KNN_traits","KNN_phylo","RF_traits","RF_phylo","IMC","GLM","fourC"),
                        labels = c("KNN","KNN (traits)", "KNN (phylo)", "RF (traits)", "RF (phylo)", "IMC", "GLM", "Fourth Corner")))

ggplot(data = probs_long, aes(x = model, y = `Model-predicted probability`, fill = as.factor(IS))) +
  geom_boxplot() + 
  scale_x_discrete(limits = rev) +
  coord_flip() +
  facet_grid(cols = vars(web)) +
  scale_fill_manual(name = "Data", 
                    values = c("white","orange"), 
                    labels = c("No interaction","Interaction")) +
  theme_pubr() +
  theme(axis.title.y = element_blank()) 


########################################################################
## Compare interaction distance vs phylogenetic distance (figs S13-S15)
########################################################################

distCors_gp <- compare_distances(alienData = galPara, troph = "from", toPlot = TRUE)
distCors_gp <- compare_distances(alienData = galPara, troph = "to", toPlot = TRUE)

distCors_sg <- compare_distances(alienData = salixGal, troph = "to", toPlot = TRUE)
distCors_sg <- compare_distances(alienData = salixGal, troph = "from", toPlot = TRUE)

###########################################################
## Compare trait importance (from Random forest) (fig S18)
###########################################################

varImpPlot(attributes(models_gp[["RF_traits"]])$model, main =NA)
varImpPlot(attributes(models_sg[["RF_traits"]])$model, main = NA)


##########################################################################
## Refitting the models with different combinations of traits (fig 4)
##########################################################################

#(Takes a long time to run)
ntraits = 4
fits_sg_4traits <- run_trait_combinations(salixGal, ntraits)
fits_gp_4traits <- run_trait_combinations(galPara, ntraits)

plot_ntrait_fits(salixGal, fits_sg_4traits)
plot_ntrait_fits(galPara, fits_gp_4traits)


##########################################################################
## Comparing performance of GLM and RF (fig S16)
##########################################################################

ggplot(fits_sg_4traits) + 
  geom_point(aes(x = GLM, y = RF, col = total.ntraits), size = 1) +
  geom_abline(slope = 1, intercept = 0) + 
  scale_color_gradient(low = "yellow", high = "slateblue", name = "Total number of trait values") +
  xlab("Model performance of GLM") + 
  ylab("Model performance of RF") +
  theme_pubr() 


ggplot(fits_gp_4traits) + 
    geom_point(aes(x = GLM, y = RF, col = total.ntraits), size = 1) +
    geom_abline(slope = 1, intercept = 0) + 
    scale_color_gradient(low = "yellow", high = "slateblue", name = "Total number of trait values") +
    xlab("Model performance of GLM") + 
    ylab("Model performance of RF") +
  theme_pubr()

##########################################################################
## Performance of RF with more trait values (Fig S16)
##########################################################################


plot(fits_sg_4traits$total.ntraits,fits_sg_4traits$RF, xlab = "Number of unique trait values", ylab = "Model fit (Tjur's D)", col = rgb(red=0,green = 0,blue = 0,alpha=0.3),pch=16)
plot(fits_gp_4traits$total.ntraits,fits_gp_4traits$RF, xlab = "Number of unique trait values", ylab = "Model fit (Tjur's D)", col = rgb(red=0,green = 0,blue = 0,alpha=0.3),pch=16)


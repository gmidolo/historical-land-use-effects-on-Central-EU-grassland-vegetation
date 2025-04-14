### R code for: "Nineteenth-century land use shape the current occurrence of some plant species, but weakly affects richness and total composition of Central European grasslands"

### Authors: Midolo, G.; Skokanova; H.; Clark, A. T.; Vymazalova, M.; Chytry, M.; Dullinger, S.; Essl, F.; Sibik, J.; Keil, P.
### Year: 2024

### Objective: To model the response of species richness to land use history categories while accounting for other covariates

#1. Set working directory & load libraries ####
username = strsplit(getwd(), '/')[[1]][3]
setwd(paste0('C:/Users/',username,'/OneDrive - CZU v Praze/czu/abulu/scr/repo/')) # Set working directory (where the .RDS file is stored)
set.seed(22) # Set seed

# R packages:
suppressWarnings({
  suppressMessages({ 
    library(vegan)
    library(lme4)
    library(rsq)
    library(performance)
    library(randomForest)
    library(pdp)
  })
})

#2. Load and prepare data ####

input.data <- readRDS('Rdata.RDS')

dat = list() # Data input for modeling
for (i in c('R1','R2','R3')) { #repeat the operation for each grassland type (dry, mesic, wet grasslands)
  
  # Select plot data
  plot.data <- input.data$plot.data[[i]] 
  
  # Select species matrix data
  species.matrix <- input.data$species.data[[i]] 
  
  # Calculate richness
  plot.data$S <- specnumber(species.matrix) # Use vegan::specnumber() to calculate species richness (S)
  
  # Log-transform plot size
  plot.data$plot_size <- log(plot.data$plot_size) 
  
  # Scale unscaled numeric predictors (`bio` are already scaled)
  plot.data$plot_size <- as.numeric(scale(plot.data$plot_size))
  plot.data$soil_pH <- as.numeric(scale(plot.data$soil_pH))
  
  # Transform categorical predictors to factor
  plot.data$historic_land_use <- as.factor(plot.data$historic_land_use)
  plot.data$country <- as.factor(plot.data$country)
  
  # Subset data for modeling
  dat[[i]] <- plot.data[,c('S','historic_land_use','plot_size','bio2','bio3','bio8','bio9','bio15','soil_pH','country')] # Subset data 
  
}


#3. Generalized Linear Mixed-Effects Models (GLMM) with observation-level random effects (OLRE) ####

# Set control structure for GLMM
cntrl = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun=2e5)) 

# Define observation-level index for OLRE 
# For additional details concerning OLRE, check Harrison (2014): https://doi.org/10.7717/peerj.616
dat.GLMM = list()
for (i in names(dat)) {
  dat.GLMM[[i]] <- dat[[i]]
  dat.GLMM[[i]]$index <- as.character(1:nrow(dat.GLMM[[i]])) 
  
}

# First, we can fit full model with GLM:
full.GLM = list()
for (i in names(dat.GLMM)) {
  full.GLM[[i]] <- glm(S ~ historic_land_use + 
                         poly(plot_size,2) + # continuous predictors are fitted with quadratic term
                         poly(bio2,2) + 
                         poly(bio3,2) + 
                         poly(bio8,2) + 
                         poly(bio9,2) + 
                         poly(bio15,2) + 
                         poly(soil_pH,2) +
                         country,
                       family = poisson, # Poisson distribution is used for count data (species number)
                       data = dat.GLMM[[i]])
}

# lapply(full.GLM, check_overdispersion) # Overdispersion is always detected in poisson GLM

# Fit full model with GLMM-OLRE via lme4::glmer(), to control for overdispersion
full.GLMM = list()
for (i in names(dat.GLMM)) { # N.B: it takes approx 4 minutes to run
  full.GLMM[[i]] <- glmer(S ~ historic_land_use + 
                            poly(plot_size,2) +
                            poly(bio2,2) + 
                            poly(bio3,2) + 
                            poly(bio8,2) + 
                            poly(bio9,2) + 
                            poly(bio15,2) + 
                            poly(soil_pH,2) + 
                            country +
                            (1|index), # OLRE is fitted with observation-level random intercepts
                          family = poisson,
                          data = dat.GLMM[[i]], 
                          control = cntrl)
}

# summary(full.GLMM[['R1']]) # Inspect results

# lapply(full.GLMM, check_overdispersion) # Overdispersion is always detected in poisson GLM!

# N.B.: In this example, homogeneity of variance is present when using GLM, while fitted vs. residuals show clearly some non independent in GLMMs, possibly due to the OLRE random structure
# We overall used both modeling approach and both yielded similar results, but report results on GLMM in the manuscript
# In general, use performance::check_model() for posterior predictive check, homogeneity of variance, and normality of the residuals plots in both model types


#4. Random Forest ####

full.RF = list()
for (i in names(dat)) { #N.B: it takes approx 12 secs to run
  full.RF[[i]] <- randomForest(S ~., data=dat[[i]]) # Fit full model with default randomForest package; i.e., mtry=3; ntree=500 
}

# Plot variable importance in randomForest
par(mfrow=c(3,1))
for (i in names(full.RF)) {
  varImpPlot(full.RF[[i]], main=paste0('Importance: ', i))
}

#5. Perform Variation Partitioning ####

# Define function for variation partitioning
VarPart.2 <- function(A,B,AB){ # Adapted from Viana et al. (2022) (https://doi.org/10.1002/ecs2.4028)
  v1 <- AB-B
  v2 <- AB-A
  v12 <- A+B-AB
  resid <- 1-AB
  vp <- c(v1, v2, v12, resid)
  names(vp) <- c('LandUse','Environment','LandUse+Env.','residual')
  return(vp)
}


# Variation Partitioning for GLMM:

VP.GLMM = list()
for (i in names(full.GLMM)) {
  #N.B. we calculate R^2 of the fixed effect using rsq package
  R2.full.GLMM <- rsq.glmm(full.GLMM[[i]], adj=T)$fixed 
  R2.full.GLMM #R^2 of the fixed effect
  
  # Fit model for all predictors other than Historical land use
  E.GLMM <- update(full.GLMM[[i]], . ~ . -historic_land_use)
  R2.E.GLMM <- rsq.glmm(E.GLMM, adj=T)$fixed # Calculate R^2
  
  # Fit model for Historical land use only (but still accounting for plot size)
  L.GLMM <- update(full.GLMM[[i]], . ~ . - poly(bio2,2) - poly(bio3,2) - poly(bio8,2) - poly(bio9,2) - poly(bio15,2) - poly(soil_pH,2))
  R2.L.GLMM <- rsq.glmm(L.GLMM, adj=T)$fixed # Calculate R^2
  
  # Var partitioning
  vp <-  VarPart.2(R2.L.GLMM, R2.E.GLMM, R2.full.GLMM)
  vp[vp<0] <- 0 #set negative values to zero, if any
  
  # Variation partition results for GLMMs
  VP.GLMM[[i]] <- vp
}

#Vizualize results
par(mfrow=c(1,3))
for (i in names(VP.GLMM)) {
  barplot(height=(VP.GLMM[[i]]*100)[-4], names=names(VP.GLMM[[i]][-4]), 
          main=paste0('Var. Part. GLMM for ', i), cex.names = .9, xlab = 'Variation explained (%)',
          las=2, horiz = TRUE)
}


# Variation Partitioning for RandomForest:

VP.RF = list()
for (i in names(full.RF)) {
#N.B. we calculate R^2 of the fixed effect using rsq package
R2.full.RF <- full.RF[[i]]$rsq[full.RF[[i]]$ntree]
R2.full.RF #R^2 of the fixed effect model

# Fit model for all predictors other than Historical land use
E.RF <- randomForest(S ~., data=dat[[i]][c('S','plot_size','bio2','bio3','bio8','bio9','bio15','soil_pH','country')])
R2.E.RF <- E.RF$rsq[E.RF$ntree] # Calculate R^2

# Fit model for Historical land use only (but still accounting for plot size)
L.RF <- randomForest(S ~., data=dat[[i]][c('S','plot_size','historic_land_use')])
R2.L.RF <- L.RF$rsq[L.RF$ntree] # Calculate R^2

# Var partitioning
vp <- VarPart.2(R2.L.RF, R2.E.RF, R2.full.RF)
vp[vp<0] <- 0 #set negative values to zero, if any

VP.RF[[i]] <- vp
}

#Vizualize results
par(mfrow=c(3,1))
for (i in names(VP.RF)) {
  barplot(height=(VP.RF[[i]]*100)[-4], names=names(VP.RF[[i]][-4]), 
          main=paste0('Var. Part. randomForest for ', i), cex.names = .9, xlab = 'Variation explained (%)',
          las=2, horiz = TRUE)
}


#6. Partial Dependence Plots - examples ####

#Partial dependence plot for historic land use, for example for R1 (dry grasslands):
#GLMM
partial(full.GLMM$R1,
        pred.var = 'historic_land_use', type = 'regression',
        inv.link = exp, prob = T, plot = T, train = dat.GLMM$R1)
#RandomForest
partial(full.RF$R1,
        pred.var = 'historic_land_use', type = 'regression',
        prob = T, plot = T, train = dat$R1)

#Partial dependence plot for continuous variables, e.g. bio15 (= precipitation seasonality) in R1 (dry grasslands)
#GLMM
partial(full.GLMM$R1,
        pred.var = 'bio15', type = 'regression',
        inv.link = exp, prob = T, plot = T,  train = dat.GLMM$R1)
#RandomForest
partial(full.RF$R1,
        pred.var = 'bio15', type = 'regression',
        prob = T, plot = T,  train = dat$R1)
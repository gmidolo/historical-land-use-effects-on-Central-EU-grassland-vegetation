### R code for: "Nineteenth-century land use shape the current occurrence of some plant species, but weakly affects richness and total composition of Central European grasslands"

### Authors: Midolo, G.; Skokanova; H.; Clark, A. T.; Vymazalova, M.; Chytry, M.; Dullinger, S.; Essl, F.; Sibik, J.; Keil, P.
### Year: 2024

### Objective: To model the response of total species composition to land use history categories while accounting for other covariates

#1. Set working directory & load libraries ####
username = strsplit(getwd(), '/')[[1]][3]
setwd(paste0('C:/Users/',username,'/OneDrive - CZU v Praze/czu/abulu/scr/repo/')) # Set working directory (where the .RDS file is stored)
set.seed(22) # Set seed

# R packages:
suppressMessages({
  library(vegan)
  })

#2. Load and prepare data ####

input.data <- readRDS('Rdata.RDS')

# Prepare distance matrices
st = Sys.time()
dm <- list() # List of distance matrices, the input of CCA
for (i in names(input.data$species.data)) {
  
  # message('Preparing distance matrix for habitat: ', i)
  
  # Select site x species matrix for each habitat
  sp.matrix = input.data$species.data[[i]]
  
  # Remove genus-only taxa (e.g., 'Achillea species') from the community matrix:
  sp.matrix <- sp.matrix[,which(!grepl(' species', colnames(sp.matrix)))]
  # dim(sp.matrix) # species (=columns) have been removed
  
  # Calculate dissimilarity matrix
  dm[[i]] <- vegdist(sp.matrix, method = 'bray')
  
  # Apply Hellinger trasformation
  dm[[i]] <- decostand(as.matrix(dm[[i]]), method = 'hellinger') 
}
Sys.time()-st # total elapsed time; aprox. 45 secs!

# Prepare input data
plot.data <- list() # List of distance matrices, the input of CCA
for (i in names(input.data$plot.data)) {
  # Select plot data
  plot.data[[i]] <- input.data$plot.data[[i]] 

  # Scale unscaled numeric predictors (`bio` are already scaled)
  plot.data[[i]]$soil_pH <- as.numeric(scale(plot.data[[i]]$soil_pH))
  
  # Transform categorical predictors to factor
  plot.data[[i]]$historic_land_use <- as.factor(plot.data[[i]]$historic_land_use)
  plot.data[[i]]$country <- as.factor(plot.data[[i]]$country)
}


#3. Fit the CCA and perform variation partitioning ####

# It takes approx. 20 minutes in total!

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

# Define the number of permutations to use when computing the adjusted R-squared for a cca
perm.R2.cca = 100 

# Fit full models and calculate r-squared

full.CCA <- list() # List to store CCA including all predictors analyzed
VP.CCA   <- list() # List to store variation partitioning via CCA
st.tot = Sys.time()
for (i in names(dm)) {
  
  st=Sys.time()
  message('Running CCA on habitat ', paste0(i, ':'))
  
  # Fit environmental-only model and calculate r-squared
  message('...fitting land-use history + environment model...')
  full.CCA.mod <- cca(dm[[i]] ~ historic_land_use + bio15 + bio2 + bio3 + bio8 + bio9 + soil_pH + country, data=plot.data[[i]])
  full.CCA.mod.R2 <- RsquareAdj(full.CCA.mod, permutations=perm.R2.cca)$adj.r.squared

  # Fit environmental-only model and calculate r-squared
  message('...fitting environment-only model... ')
  E.CCA.mod <- cca(dm[[i]] ~  bio15 + bio2 + bio3 + bio8 + bio9 + soil_pH + country, data=plot.data[[i]])
  E.CCA.mod.R2 <- RsquareAdj(E.CCA.mod, permutations=perm.R2.cca)$adj.r.squared

  # Fit historical land use -only model and calculate r-squared
  message('...fitting land-use history -only model...')
  L.CCA.mod <- cca(dm[[i]] ~ historic_land_use, data=plot.data[[i]])
  L.CCA.mod.R2 <- RsquareAdj(L.CCA.mod, permutations=perm.R2.cca)$adj.r.squared
  
  # Caclulate variation partitioning and store the results
  VP.CCA[[i]] <- VarPart.2(L.CCA.mod.R2, E.CCA.mod.R2, full.CCA.mod.R2)
  VP.CCA[[i]][VP.CCA[[i]]<0] <- 0 #set negative values to zero, if any
  full.CCA[[i]] <- full.CCA.mod
  
  message(paste0('Habitat ', i, ': DONE.'))
  print(Sys.time()-st)
}
Sys.time() - st.tot # total elapsed time


#4. Visualize results ####

# plot CCA (with vegan)
par(mfrow=c(3,1), mar=c(1,2,1,3))
for (i in names(VP.CCA)) {
  plot(full.CCA[[i]])
  mtext(paste0('CCA ', i), side=4, padj=.5)
}

# Vizualize Variation partitioning results
par(mfrow=c(3,1), mar = rep(5,4))
for (i in names(VP.CCA)) {
  barplot(height=(VP.CCA[[i]]*100)[-4], 
          names= c('land','env','env+land'), 
          main=paste0('Var. Part. CCA for ', i), 
          ylab = 'Variation explained (%)', 
          las=2)
}



#5. Fit PERMANOVA via adonis2 ####

perm.adonis2 = 100 # Define the number of permutations in adonis2

PRMNV <- list() # List to store results of adonis2
st.tot = Sys.time()
for (i in names(dm)) {
  st=Sys.time()
  message('Running PERMANOVA on habitat ', paste0(i, '...'))
  PRMNV[[i]] <- adonis2(dm[[i]] ~ historic_land_use + bio15 + bio2 + bio3 + bio8 + bio9 + soil_pH + country, data=plot.data[[i]], permutations = perm.adonis2)
  message(paste0('Habitat ', i, ': DONE.'))
  print(Sys.time()-st)
}
Sys.time() - st.tot # total elapsed time

# lapply(PRMNV.mod, print) # inspect adonis2 results
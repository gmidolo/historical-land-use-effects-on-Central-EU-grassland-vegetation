### R code for: "Nineteenth-century land use shape the current occurrence of some plant species, but weakly affects richness and total composition of Central European grasslands"

### Authors: Midolo, G.; Skokanova; H.; Clark, A. T.; Vymazalova, M.; Chytry, M.; Dullinger, S.; Essl, F.; Sibik, J.; Keil, P.
### Year: 2024

### Objective: Assign species to historical land use categories based on IndVal statistic and check differences in ecological and disturbance indicator values

#1. Set working directory & load libraries ####
username = strsplit(getwd(), '/')[[1]][3]
setwd(paste0('C:/Users/',username,'/OneDrive - CZU v Praze/czu/abulu/scr/repo/')) # Set working directory (where the .RDS file is stored)
set.seed(22) # Set seed

# R packages:
suppressMessages({
  library(indicspecies)
  library(pracma)
  library(multcompView)
})


#2. Load data ####
input.data <- readRDS('Rdata.RDS')


#3. Calculate IndVal for each habitat####
metric = 'IndVal.g' # Metric employed in indicspecies::multipatt()
nperm = 999 # Number of permutations in indicspecies::multipatt()
res.IndVal <- list()
for (i in names(input.data$species.data)) {
  
  message('Calculate IndVal for ',i)
  
  # Select site x species matrix for each habitat
  species.matrix <- input.data$species.data[[i]]
  
  # Remove genus-only taxa (e.g., 'Achillea species') from the community matrix:
  species.matrix <- species.matrix[,which(!grepl('.species', colnames(species.matrix)))]
  
  # Select plot data
  plot.data <- input.data$plot.data[[i]]
  
  all(plot.data$plot_id == as.numeric(rownames(species.matrix))) #all vegetation plots id are in the same order?
  
  st = Sys.time()
  res.IndVal[[i]] <- multipatt(species.matrix, # Use the multipatt function from the `indicspecies` package
                               cluster = plot.data$historic_land_use,
                               func = metric,
                               max.order = 1,
                               control = how(nperm=nperm))
  print(Sys.time() - st)
  
}


#4. Tidy the results of IndVal from multipatt and merge indicator values ####

# Select indicator data (PCA axes of ecological and disturbance indicator values)
ind.pca.axis <- input.data$indic.val[,c('Species', paste0('RC', 1:4))]
head(ind.pca.axis)

# Tidy IndVal data
res = list()
for (i in names(res.IndVal)) {
  # Extract restults from multipatt()
  res[[i]] <- res.IndVal[[i]]$sign 
  # Identify species
  res[[i]]$Species <- rownames(res[[i]])
  # Identify historical land-use categories present in the data
  hlu.names <- names(res[[i]])[grepl('s\\.',names(res[[i]]))] 
  # Reshape the dataframe
  res[[i]] <- reshape(res[[i]], 
                     varying = hlu.names, 
                     timevar = 'hlu.code', v.names = 'value', direction = 'long') 
  # Assign historical land-use categories names
  conversion.land.use <- data.frame(hlu.code = 1:length(hlu.names), Historic_land_use=sub('s\\.', '', hlu.names)) 
  res[[i]] <- merge(res[[i]], conversion.land.use, by='hlu.code', all.x=TRUE)
  # Filter only the historical land-use categories with the highest IndVal for each species
  res[[i]] <- res[[i]][res[[i]]$value == 1, ] 
  # Select columns needed
  res[[i]] <- res[[i]][,c('Species', 'Historic_land_use', 'stat', 'p.value')] 
  # Set colnames
  names(res[[i]]) <- c('Species', 'Historic_land_use', metric,' p.value')
  # Remove row names
  rownames(res[[i]]) = NULL 
  # Merge the rotated component values with the IndVal result table
  res[[i]] <- merge(res[[i]], ind.pca.axis, by='Species', all.x=TRUE)
  # Only retain species with data on indicator values (24 species are dropped):
  res[[i]] <- na.omit(res[[i]]) 
}

head(res$R1) # head(), for e.g. Dry grasslands


#6. Perform ANOVA and Tukey test on the indicator values ####

weighted.ANOVA = T # weight the anova by IndVal ?
par(mfrow=c(2,2), mar = c(4,4,4,4))
col.palette = data.frame(
  name = c('arable.land','forest','grassland','permanent.crop','settlement','water.body'),
  col = c('darkgoldenrod1', 'springgreen4', 'chartreuse', 'cornsilk3', 'darkred', 'deepskyblue1')
)

for (i in names(res)) {

  for (k in paste0('RC',1:4)) {
    dfk <- res[[i]]
    names(dfk)[which(names(dfk)%in%k)] <- 'y'
    
    if(weighted.ANOVA){
      anova <- aov(y ~ Historic_land_use, data = dfk, weights = IndVal.g) # Perform ANOVA weighted by IndVal
    } else {
      anova <- aov(y ~ Historic_land_use, data = dfk) # Perform unweighted standard ANOVA
    }
    
    tukey <- TukeyHSD(anova) # Perform TukeyHSD
    cld <- multcompLetters4(anova, tukey) # Extracting the compact letter to display in boxplots
    cld <- as.data.frame.list(cld[[1]]) # Transoform this to data
    cld$Historic_land_use <- rownames(cld) # Create names of historical land use
    cld <- cld[, c('Historic_land_use','Letters')] # Select columns needed
    cld <- cld[order(rownames(cld)), ] # sort alphabetically
    
    # Calculate standard quantile for text positioning in the plots
    dfklist <- split(dfk, dfk$Historic_land_use)
    qntls = list()
    for (q in names(dfklist)) {
      qntls[[q]] <- data.frame(Historic_land_use = q, quant = quantile(dfklist[[q]]$y, probs = .85))
    }
    
    boxplot(y ~ Historic_land_use, data = dfk, notch = T,
            xlab='', ylab = k, las=2, main=paste0('Habitat: ', i),
            col = col.palette$col[col.palette$name %in% sort(cld$Historic_land_use)])
    text(x=1:length(qntls), y=do.call('rbind', qntls)$quant, labels=cld$Letters, font = 2, adj=c(-.75,-.5))
    
  }
}

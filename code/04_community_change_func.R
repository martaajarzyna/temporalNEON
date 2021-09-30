# Load packages
require(rhdf5)
require(here)
require(tidyverse)
require(BAT)
require(codyn)
require(data.table)
select <- dplyr::select
summarize <- dplyr::summarize
# functional diversity packages
require(FD)
require(ape)
require(phytools)
require(picante)


################### CHANGE IN FUNCTIONAL COMMUNITY COMPOSITION

################### small mammals
data.mam.bout <- readRDS(file="data/mammal_abund-bout.rds")
sites.id <- unique(data.mam.bout$siteID)

###### Get trait data
spnames <- as_tibble(colnames(data.mam.bout[,8:ncol(data.mam.bout)]))
colnames(spnames) <- "scientificName"
# Read in taxonomy from Elton traits
spnames.et <- read.csv(file="traits/MammalTraits_taxonomy.csv", header=TRUE, sep=",")
spnames.et <- spnames.et %>%
  mutate(Scientific_name = paste0(Genus," ",Species))

# check for taxonomic agreement
spnames.et.sel <- spnames %>%
  filter(!scientificName %in% spnames.et$Scientific_name) #select names from neon data that are absent from EltonTrait taxonomy
spnames.et.sel
# "Ictidomys tridecemlineatus" - this species is listed in EltonTraits taxonomy as "Spermophilus tridecemlineatus"
# Replace "Spermophilus tridecemlineatus" in EltonTraits with "Ictidomys tridecemlineatus"
spnames.et[spnames.et$Scientific_name == "Spermophilus tridecemlineatus",5] <- "Ictidomys tridecemlineatus"

# Read in EltonTraits
# Data from: https://esajournals.onlinelibrary.wiley.com/doi/10.1890/13-1917.1 
traits <- read.csv(file="traits/MammalTraits_clean.csv", header=TRUE, sep=",")
head(traits)

traits <- traits %>%
  mutate(Scientific_name = paste0(Genus," ",Species)) %>%
  select(MSW3_ID,Scientific_name,Diet.Inv,Diet.Vend,Diet.Vect,Diet.Vfish,Diet.Vunk,Diet.Scav,Diet.Fruit,Diet.Nect,
         Diet.Seed,Diet.PlantO,ForStrat.Value,Activity.Nocturnal,Activity.Crepuscular,Activity.Diurnal,BodyMass.Value) %>%
  arrange(Scientific_name)
# Replace "Spermophilus tridecemlineatus" in EltonTraits with "Ictidomys tridecemlineatus"
traits[traits$Scientific_name == "Spermophilus tridecemlineatus",2] <- "Ictidomys tridecemlineatus"

traits.sel <- traits %>%
  filter(Scientific_name %in% spnames$scientificName) %>%
  arrange(Scientific_name)

spnames <- spnames %>%
  arrange(scientificName)

# Order data table in order of the trait table
data.mam.bout <- setcolorder(data.mam.bout, 
                             as.character(c("domainID","siteID","lat","lon","year","month","bout", 
                                            spnames$scientificName)))

### Compute functional diversity
sp <- traits.sel[,2]
func.dat <- subset(traits.sel, select = -c(MSW3_ID, Scientific_name))
Weights=c(rep(0.25/10,10),0.25,rep(0.25/3,3),0.25) 
# for Diet(num,10), ForagingStratum(num,1),ActivityTime(num,3),Mass(num,1)
# Create functional dissimilarity matrix using gower distance and weights and perform PCoA
func.dist <- gowdis(func.dat, ord = "podani", w = Weights)
# Perform clustering using UPGMA
clust.obj <- hclust(func.dist,  method = "average")
clust.obj$labels <- sp
# Calculate tree based on UPGMA clustering
tree.mam <- as.phylo(clust.obj)
ape::write.tree(tree.mam, file = "data/mammal_func-tree")

# Quantify correlation between functional distances and cophenetic distances on a dendrogram to ensure that functional space is properly captured by the dendrogram
cor(c(as.matrix(func.dist)), c(cophenetic.phylo(tree.mam))) #correlation >0.93
# Quanitfy mSD metric from Maire
dist1 <- as.matrix(func.dist)
dist2 <- cophenetic.phylo(tree.mam)
dist1[upper.tri(dist1)] <- NA
dist2[upper.tri(dist2)] <- NA
dist.all <- cbind(melt(dist1), melt(dist2))
colnames(dist.all) <- c("sp1.1","sp2.1","dist_coph","sp1.2","sp2.2","dist_gower")
dist.all <- dist.all %>%
  dplyr::filter(dist_coph != "NA") %>%
  dplyr::filter(dist_gower != "NA")
dist.all <- dist.all %>%
  mutate(diff = (dist_gower - dist_coph)^2)
msd <- sum(dist.all$diff)/((length(unique(dist.all$sp1.2))*(length(unique(dist.all$sp1.2))-1))/2) #msd=0.022


###### Functional community composition change: Package BAT
for (i in 1:length(sites.id)){
  mi <- data.mam.bout %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  yi <- nrow(mi)-1
  yi_ls <- list()
  for (p in 1:yi){
    ps <- seq(1,p,1)
    yi_ls[[(yi+1)-p]] <- ps
  }
  lags <- melt(yi_ls)
  
  betadiv <- beta(mi[,8:ncol(mi)], tree=tree.mam, abund = TRUE, raref = 0)
  out.bout.bat <- matrix(NA, nrow(lags), 8)
  out.bout.bat <- as.data.frame(out.bout.bat)
  jacc <- as.matrix(betadiv$Btotal)
  jacc[upper.tri(jacc)] <- NA
  diag(jacc) <- NA
  jacc <- melt(jacc)
  jacc <- jacc[complete.cases(jacc),]
  
  jaccrepl <- as.matrix(betadiv$Brepl)
  jaccrepl[upper.tri(jaccrepl)] <- NA
  diag(jaccrepl) <- NA
  jaccrepl <- melt(jaccrepl)
  jaccrepl <- jaccrepl[complete.cases(jaccrepl),]
  
  jaccrich <- as.matrix(betadiv$Brich)
  jaccrich[upper.tri(jaccrich)] <- NA
  diag(jaccrich) <- NA
  jaccrich <- melt(jaccrich)
  jaccrich <- jaccrich[complete.cases(jaccrich),]
  
  out.bout.bat[,1] <- c(domain)
  out.bout.bat[,2] <- sites.id[i]
  out.bout.bat[,3] <- lat
  out.bout.bat[,4] <- lon
  out.bout.bat[,5] <- lags[,1]
  out.bout.bat[,6] <- jacc[,3]
  out.bout.bat[,7] <- jaccrepl[,3]
  out.bout.bat[,8] <- jaccrich[,3]
  
  if (i == 1) {
    out.bout.bat.all <- out.bout.bat
  } else {
    out.bout.bat.all <- out.bout.bat.all %>%
      bind_rows(out.bout.bat)
  }
}

colnames(out.bout.bat.all) <- c("domainID","siteID","lat","lon","bout_lag","bout_consec_jacc","bout_consec_jaccrepl","bout_consec_jaccrich")
out.bout.bat.all <- out.bout.bat.all %>% 
  mutate(bout_consec_contr_jaccrepl = bout_consec_jaccrepl/bout_consec_jacc, bout_consec_contr_jaccrich = bout_consec_jaccrich/bout_consec_jacc)
saveRDS(out.bout.bat.all, file="data/mammal_outBAT-FD-bout-consecutive-all.rds")


########read in all data 
stab.bout.cons.all <- readRDS(file="data/mammal_outStability-bout-consecutive-all.rds")
stab.bout.cons.all$bout_lag <- as.numeric(stab.bout.cons.all$bout_lag)
bat.bout.cons.all <- readRDS(file="data/mammal_outBAT-FD-bout-consecutive-all.rds")
bat.bout.cons.all$bout_lag <- as.numeric(bat.bout.cons.all$bout_lag)

all.bout.cons.all <- cbind(stab.bout.cons.all, bat.bout.cons.all[,6:ncol(bat.bout.cons.all)])
saveRDS(all.bout.cons.all, file="data/mammal_outAll-FD-bout-consecutive-all.rds")




################### fish
data.fish.bout <- readRDS(file="data/fish_abund-bout.rds")
sites.id <- unique(data.fish.bout$siteID)

###### Get trait data
spnames <- as_tibble(colnames(data.fish.bout[,8:ncol(data.fish.bout)]))
colnames(spnames) <- "scientificName"

# Read in Traits and taxonomy
# Data from: https://esajournals.onlinelibrary.wiley.com/doi/10.1890/13-1917.1 
spnames.et <- read.csv(file="traits/FishTraits_clean.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
#some species names are misspelled
spnames.et[69,2] <- "Notropis suttkusi"
spnames.et[79,2] <- "Oncorhynchus clarki"
spnames.et[105,2] <- "Sicydium plumieri"
spnames.et <- spnames.et %>%
  dplyr::rename(scientificName = species)

# check for taxonomic agreement
spnames.et.sel <- spnames %>%
  filter(!scientificName %in% spnames.et$scientificName) #select names from neon data that are absent from EltonTrait taxonomy
spnames.et.sel
# all names are aligned 

# traits
traits <- spnames.et

#filter traits data to only those present in neon
traits <- traits %>%
  filter(scientificName %in% spnames$scientificName) #seem like all names are in agreement
# select traits of interest
traits.sel <- traits %>%
  select(rowid,scientificName,algphyto,macvascu,detritus,invlvfsh,fshcrcrb,maxtl,matuage,longevity,fecundity,muck,claysilt,
         sand,gravel,cobble,boulder,bedrock,vegetat,debrdetr,lwd,potanadr) %>%
  arrange(scientificName)

spnames <- spnames %>%
  arrange(scientificName)

# Order data table in order of the trait table
data.fish.bout <- setcolorder(data.fish.bout, 
                              as.character(c("domainID","siteID","lat","lon","year","month","bout", 
                                             spnames$scientificName)))

### Compute functional diversity
sp <- traits.sel[,2]
func.dat <- subset(traits.sel, select = -c(rowid, scientificName))
func.dat[,c(1:5,10:ncol(func.dat))] <- lapply(func.dat[,c(1:5,10:ncol(func.dat))], as.numeric)
Weights=c(rep(1/5,5),1,1,1,1,rep(1/10,10),1) 
# because some traits are within the same category Weights=c(rep(1,ncol(funcdat)))
# because guilds are not mutually exclusive, we assign all traits a weight of 1 
func.dist <- gowdis(func.dat, ord = "podani", w = Weights)
# Perform clustering using UPGMA
clust.obj <- hclust(func.dist,  method = "average")
clust.obj$labels <- sp
# Calculate tree based on UPGMA clustering
tree.fish <- as.phylo(clust.obj)
ape::write.tree(tree.fish, file = "data/fish_func-tree")

# Quantify correlation between functional distances and cophenetic distances on a dendrogram to ensure that functional space is properly captured by the dendrogram
cor(c(as.matrix(func.dist)), c(cophenetic.phylo(tree.fish))) #correlation >0.93
# Quanitfy mSD metric from Maire
dist1 <- as.matrix(func.dist)
dist2 <- cophenetic.phylo(tree.fish)
dist1[upper.tri(dist1)] <- NA
dist2[upper.tri(dist2)] <- NA
dist.all <- cbind(melt(dist1), melt(dist2))
colnames(dist.all) <- c("sp1.1","sp2.1","dist_coph","sp1.2","sp2.2","dist_gower")
dist.all <- dist.all %>%
  dplyr::filter(dist_coph != "NA") %>%
  dplyr::filter(dist_gower != "NA")
dist.all <- dist.all %>%
  mutate(diff = (dist_gower - dist_coph)^2)
msd <- sum(dist.all$diff)/((length(unique(dist.all$sp1.2))*(length(unique(dist.all$sp1.2))-1))/2) #msd=0.022


###### Functional community composition change: Package BAT
for (i in 1:length(sites.id)){
  mi <- data.fish.bout %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  yi <- nrow(mi)-1
  yi_ls <- list()
  for (p in 1:yi){
    ps <- seq(1,p,1)
    yi_ls[[(yi+1)-p]] <- ps
  }
  lags <- melt(yi_ls)
  
  betadiv <- beta(mi[,8:ncol(mi)], tree=tree.fish, abund = TRUE, raref = 0)
  out.bout.bat <- matrix(NA, nrow(lags), 8)
  out.bout.bat <- as.data.frame(out.bout.bat)
  jacc <- as.matrix(betadiv$Btotal)
  jacc[upper.tri(jacc)] <- NA
  diag(jacc) <- NA
  jacc <- melt(jacc)
  jacc <- jacc[complete.cases(jacc),]
  
  jaccrepl <- as.matrix(betadiv$Brepl)
  jaccrepl[upper.tri(jaccrepl)] <- NA
  diag(jaccrepl) <- NA
  jaccrepl <- melt(jaccrepl)
  jaccrepl <- jaccrepl[complete.cases(jaccrepl),]
  
  jaccrich <- as.matrix(betadiv$Brich)
  jaccrich[upper.tri(jaccrich)] <- NA
  diag(jaccrich) <- NA
  jaccrich <- melt(jaccrich)
  jaccrich <- jaccrich[complete.cases(jaccrich),]
  
  out.bout.bat[,1] <- c(domain)
  out.bout.bat[,2] <- sites.id[i]
  out.bout.bat[,3] <- lat
  out.bout.bat[,4] <- lon
  out.bout.bat[,5] <- lags[,1]
  out.bout.bat[,6] <- jacc[,3]
  out.bout.bat[,7] <- jaccrepl[,3]
  out.bout.bat[,8] <- jaccrich[,3]
  
  if (i == 1) {
    out.bout.bat.all <- out.bout.bat
  } else {
    out.bout.bat.all <- out.bout.bat.all %>%
      bind_rows(out.bout.bat)
  }
}

colnames(out.bout.bat.all) <- c("domainID","siteID","lat","lon","bout_lag","bout_consec_jacc","bout_consec_jaccrepl","bout_consec_jaccrich")
out.bout.bat.all <- out.bout.bat.all %>% 
  mutate(bout_consec_contr_jaccrepl = bout_consec_jaccrepl/bout_consec_jacc, bout_consec_contr_jaccrich = bout_consec_jaccrich/bout_consec_jacc)
saveRDS(out.bout.bat.all, file="data/fish_outBAT-FD-bout-consecutive-all.rds")


########read in all data 
stab.bout.cons.all <- readRDS(file="data/fish_outStability-bout-consecutive-all.rds")
stab.bout.cons.all$bout_lag <- as.numeric(stab.bout.cons.all$bout_lag)
bat.bout.cons.all <- readRDS(file="data/fish_outBAT-FD-bout-consecutive-all.rds")
bat.bout.cons.all$bout_lag <- as.numeric(bat.bout.cons.all$bout_lag)

all.bout.cons.all <- cbind(stab.bout.cons.all, bat.bout.cons.all[,6:ncol(bat.bout.cons.all)])
saveRDS(all.bout.cons.all, file="data/fish_outAll-FD-bout-consecutive-all.rds")

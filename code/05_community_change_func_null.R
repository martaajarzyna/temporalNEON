# Load packages
require(rhdf5)
require(here)
require(tidyverse)
require(BAT)
require(codyn)
require(data.table)
select <- dplyr::select
# functional diversity packages
require(FD)
require(ape)
require(phytools)
require(picante)


################### CHANGE IN FUNCTIONAL COMMUNITY COMPOSITION: NULL EXPECTATION

################### small mammals
data.mam.bout <- readRDS(file="data/mammal_abund-bout.rds")
data.mam.year <- readRDS(file="data/mammal_abund-year.rds")
sites.id <- unique(data.mam.bout$siteID)

### Reshuffle (1000x) species abundances, separately for each domain because domain will be used as a delineation of the regional species pool
# first order by species name (column)
cols <- sort(colnames(data.mam.bout[,8:ncol(data.mam.bout)]))
cols.bout <- c(colnames(data.mam.bout[,1:7]),cols)
cols.year <- c(colnames(data.mam.year[,1:5]),cols)
data.mam.bout <- data.mam.bout[cols.bout]
data.mam.year <- data.mam.year[cols.year]

# get all species ever recorded and then get list for each domain
data_small_mammal <- readRDS(file = "data/mammal_processed.rds")
domain.sp <- data_small_mammal %>%
  select(domainID, siteID, scientificName) %>%
  group_by(domainID, scientificName) %>%
  tally(name = "count") %>% 
  ungroup()

domains <- unique(domain.sp$domainID)

# bout-level
rand.bout <- list()

for (i in 1:1000){
  cat(i)
  rand.j <- list()
  for (j in 1:length(domains)){
    data.j <- data.mam.bout %>% 
      filter(domainID == domains[j])
    domain.j <- domain.sp %>% 
      filter(domainID == domains[j])
    
    ids <- data.j %>%
      select(domainID, siteID, lat, lon, year, month, bout)
    species <- data.j %>%
      select(-c(domainID, siteID, lat, lon, year, month, bout))
    species.dom <- species[,colnames(species) %in% domain.j$scientificName]
    species.ndom <- species[,!(colnames(species) %in% domain.j$scientificName)]

colnames(species.dom) <- sample(colnames(species.dom)) #sampling randomly
species.all <- species.ndom %>%
  bind_cols(species.dom)
species.all <- species.all[,sort(names(species.all))]
species.all <- ids %>%
  bind_cols(species.all)
rand.j <- rand.j %>%
  bind_rows(species.all)
}
rand.bout[[i]] <- rand.j
}

saveRDS(rand.bout, file="data/mammal_abund-bout_null.rds")

# year-level
rand.year <- list()

for (i in 1:1000){
  cat(i)
  rand.j <- list()
  for (j in 1:length(domains)){
    data.j <- data.mam.year %>% 
      filter(domainID == domains[j])
    domain.j <- domain.sp %>% 
      filter(domainID == domains[j])
    
    ids <- data.j %>%
      select(domainID, siteID, lat, lon, year)
    species <- data.j %>%
      select(-c(domainID, siteID, lat, lon, year))
    species.dom <- species[,colnames(species) %in% domain.j$scientificName]
    species.ndom <- species[,!(colnames(species) %in% domain.j$scientificName)]
    
    colnames(species.dom) <- sample(colnames(species.dom)) #sampling randomly
    species.all <- species.ndom %>%
      bind_cols(species.dom)
    species.all <- species.all[,sort(names(species.all))]
    species.all <- ids %>%
      bind_cols(species.all)
    rand.j <- rand.j %>%
      bind_rows(species.all)
  }
  rand.year[[i]] <- rand.j
}


saveRDS(rand.year, file="data/mammal_abund-year_null.rds")


###### Null functional community composition change: Package BAT
### Get trait data (same code as in 04_community_change_func.R)
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
data.mam.year <- setcolorder(data.mam.year, 
                             as.character(c("domainID","siteID","lat","lon","year", 
                                            spnames$scientificName)))

### Compute functional diversity
sp <- traits.sel[,2]
func.dat <- subset(traits.sel, select = -c(MSW3_ID, Scientific_name))
Weights=c(rep(0.25/10,10),0.25,rep(0.25/3,3),0.25) # for Diet(num,10), ForagingStratum(num,1),ActivityTime(num,3),Mass(num,1)
# Create functional dissimilarity matrix using gower distance and weights and perform PCoA
func.dist <- gowdis(func.dat, ord = "podani", w = Weights)
# Perform clustering using UPGMA
clust.obj <- hclust(func.dist,  method = "average")
clust.obj$labels <- sp
# Calculate tree based on UPGMA clustering
tree.mam <- as.phylo(clust.obj)

### Read in species list randomizations
rand.bout <- readRDS(file="data/mammal_abund-bout_null.rds")
rand.year <- readRDS(file="data/mammal_abund-year_null.rds")

### Function to retrieve the diagonal of a distance matrix
### courtesy of https://stackoverflow.com/questions/39231961/extract-diagonals-from-a-distance-matrix-in-r
### retain only diagonal because we only want boout-to-boout and year-to-year comparison not all pair-wise combinations, to avoid biases associated with length of time series
subdiag <- function (dist_obj, d) {
  if (!inherits(dist_obj, "dist")) stop("please provide a 'dist' object!")
  n <- attr(dist_obj, "Size")
  if (d < 1 || d > (n - 1)) stop(sprintf("'d' needs be between 1 and %d", n - 1L))
  j_1 <- c(0, seq.int(from = n - 1, by = -1, length = n - d - 1))
  subdiag_ind <- d + cumsum(j_1)
  dist_obj[subdiag_ind]
}

null.bout <- list()
null.year <- list()

### Intra-annual (i.e., across bouts) change in community functional composition
### Loop over randomizations
for (k in 501:1000){
  cat(k)
  data.mam.bout <- rand.bout[[k]]
  sites.id <- unique(data.mam.bout$siteID)
  
  for (i in 1:length(sites.id)) {
    mi <- data.mam.bout %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  yi <- unique(mi$year)
  
  for (j in 1:length(yi)){
    mij <- mi %>%
      filter(year == yi[j])
    
    out.intra.bat <- matrix(NA, (nrow(mij)-1), 12)
    out.intra.bat <- as.data.frame(out.intra.bat)
    out.intra.bat[,1] <- c(domain)
    out.intra.bat[,2] <- sites.id[i]
    out.intra.bat[,3] <- lat
    out.intra.bat[,4] <- lon
    out.intra.bat[,5] <- yi[j]
    betadiv <- beta(mij[,8:ncol(mij)], tree=tree.mam, abund = TRUE, raref = 0)
    bouts <- mij[,7]
    out.intra.bat[,6] <- bouts[1:(nrow(bouts)-1),] #bout 1
    out.intra.bat[,7] <- bouts[2:nrow(bouts),] #bout 2
    out.intra.bat[,8] <- subdiag(betadiv$Btotal, 1)
    out.intra.bat[,9] <- subdiag(betadiv$Brepl, 1)
    out.intra.bat[,10] <- subdiag(betadiv$Brich, 1)
    out.intra.bat <- out.intra.bat %>% 
      mutate(V11 = V9/V8, V12 = V10/V8)
    
    if (j == 1) {
      out.intra.bat.ally <- out.intra.bat
    } else {
      out.intra.bat.ally <- out.intra.bat.ally %>%
        bind_rows(out.intra.bat)
    } 
  }
  
  if (i == 1) {
    out.intra.bat.all <- out.intra.bat.ally
  } else {
    out.intra.bat.all <- out.intra.bat.all %>%
      bind_rows(out.intra.bat.ally)
  } 
}

colnames(out.intra.bat.all) <- c("domainID","siteID","lat","lon","year","bout_from","bout_to","bout_jacc","bout_jaccrepl","bout_jaccrich","bout_contr_jaccrepl","bout_contr_jaccrich")
null.bout[[k]] <- out.intra.bat.all
}

saveRDS(null.bout, file="data/mammal_outBAT-FD-bout_null.rds")

### Inter-annual (i.e., across years) change in community functional composition
### Loop over randomizations
for (k in 1:1000){
  cat(k)
  data.mam.year <- rand.year[[k]]
  sites.id <- unique(data.mam.year$siteID)
  
  for (i in 1:length(sites.id)){
    mi <- data.mam.year %>%
      filter(siteID == sites.id[i])
    domain <- unique(mi$domainID)
    lat <- mean(mi$lat)
    lon <- mean(mi$lon)
    yi <- unique(mi$year)
    
    betadiv <- beta(mi[,6:ncol(mi)], tree=tree.mam, abund = TRUE, raref = 0)
    out.inter.bat <- matrix(NA, (nrow(mi)-1), 11)
    out.inter.bat <- as.data.frame(out.inter.bat)
    out.inter.bat[,1] <- c(domain)
    out.inter.bat[,2] <- sites.id[i]
    out.inter.bat[,3] <- lat
    out.inter.bat[,4] <- lon
    out.inter.bat[,5] <- yi[1:(length(yi)-1)]
    out.inter.bat[,6] <- yi[2:length(yi)]
    out.inter.bat[,7] <- subdiag(betadiv$Btotal, 1)
    out.inter.bat[,8] <- subdiag(betadiv$Brepl, 1)
    out.inter.bat[,9] <- subdiag(betadiv$Brich, 1)
    out.inter.bat <- out.inter.bat %>% 
      mutate(V10 = V8/V7, V11 = V9/V7)
    
    if (i == 1) {
      out.inter.bat.all <- out.inter.bat
    } else {
      out.inter.bat.all <- out.inter.bat.all %>%
        bind_rows(out.inter.bat)
    }
  }
  
colnames(out.inter.bat.all) <- c("domainID","siteID","lat","lon","year_from","year_to","year_jacc","year_jaccrepl","year_jaccrich","year_contr_jaccrepl","year_contr_jaccrich")
null.year[[k]] <- out.inter.bat.all
}

saveRDS(null.year, file="data/mammal_outBAT-FD-year_null.rds")



###### Quantify the standardized effect size (SES)
obs.bout <- readRDS(file="data/mammal_outBAT-FD-bout.rds")
obs.year <- readRDS(file="data/mammal_outBAT-FD-year.rds")
null.bout <- readRDS(file="data/mammal_outBAT-FD-bout_null.rds")
null.year <- readRDS(file="data/mammal_outBAT-FD-year_null.rds")
obs.bout <- obs.bout %>% 
  mutate(bout_contr_jaccrepl = bout_jaccrepl/bout_jacc, bout_contr_jaccrich = bout_jaccrich/bout_jacc)
obs.year <- obs.year %>% 
  mutate(year_contr_jaccrepl = year_jaccrepl/year_jacc, year_contr_jaccrich = year_jaccrich/year_jacc)

### bout-level
# Jaccard dissimilarity
null.jacc.bout <- obs.bout %>%
  select(domainID, siteID, lat, lon, year, bout_from, bout_to, bout_jacc)
for (i in 1:1000){
  null.jacc.i <- null.bout[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(bout_jacc)
  null.jacc.bout <- null.jacc.bout %>%
    bind_cols(null.jacc.i)
}

null.jacc.bout$bout_jaccses <- NA
null.jacc.bout$bout_jaccrank <- NA
null.jacc.bout$bout_jaccpval <- NA
for (i in 1:nrow(null.jacc.bout)){
  null.jacc.bout[i,1009] <- ((as.numeric(null.jacc.bout[i,8]) - mean(as.numeric(null.jacc.bout[i,9:1008])))/sd(as.numeric(null.jacc.bout[i,9:1008])))
  rank.sim <- apply(null.jacc.bout[i,8:1008],1,rank)[,1]
  null.jacc.bout[i,1010] <- rank.sim[1]
  null.jacc.bout[i,1011] <- null.jacc.bout[i,1010]/1001
}

# Jaccard dissimilarity due to replacement 
null.jaccrepl.bout <- obs.bout %>%
  select(domainID, siteID, lat, lon, year, bout_from, bout_to, bout_jaccrepl)
for (i in 1:1000){
  null.jacc.i <- null.bout[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(bout_jaccrepl)
  null.jaccrepl.bout <- null.jaccrepl.bout %>%
    bind_cols(null.jacc.i)
}

null.jaccrepl.bout$bout_jaccreplses <- NA
null.jaccrepl.bout$bout_jaccreplrank <- NA
null.jaccrepl.bout$bout_jaccreplpval <- NA
for (i in 1:nrow(null.jaccrepl.bout)){
  null.jaccrepl.bout[i,1009] <- ((as.numeric(null.jaccrepl.bout[i,8]) - mean(as.numeric(null.jaccrepl.bout[i,9:1008])))/sd(as.numeric(null.jaccrepl.bout[i,9:1008])))
  rank.sim <- apply(null.jaccrepl.bout[i,8:1008],1,rank)[,1]
  null.jaccrepl.bout[i,1010] <- rank.sim[1]
  null.jaccrepl.bout[i,1011] <- null.jaccrepl.bout[i,1010]/1001
}

# Jaccard dissimilarity due to richness 
null.jaccrich.bout <- obs.bout %>%
  select(domainID, siteID, lat, lon, year, bout_from, bout_to, bout_jaccrich)
for (i in 1:1000){
  null.jacc.i <- null.bout[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(bout_jaccrich)
  null.jaccrich.bout <- null.jaccrich.bout %>%
    bind_cols(null.jacc.i)
}

null.jaccrich.bout$bout_jaccrichses <- NA
null.jaccrich.bout$bout_jaccrichrank <- NA
null.jaccrich.bout$bout_jaccrichpval <- NA
for (i in 1:nrow(null.jaccrich.bout)){
  null.jaccrich.bout[i,1009] <- ((as.numeric(null.jaccrich.bout[i,8]) - mean(as.numeric(null.jaccrich.bout[i,9:1008])))/sd(as.numeric(null.jaccrich.bout[i,9:1008])))
  rank.sim <- apply(null.jaccrich.bout[i,8:1008],1,rank)[,1]
  null.jaccrich.bout[i,1010] <- rank.sim[1]
  null.jaccrich.bout[i,1011] <- null.jaccrich.bout[i,1010]/1001
}

# Contribution of replacement to Jaccard dissimilarity
null.jaccreplcontr.bout <- obs.bout %>%
  select(domainID, siteID, lat, lon, year, bout_from, bout_to, bout_contr_jaccrepl)
for (i in 1:1000){
  null.jacc.i <- null.bout[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(bout_contr_jaccrepl)
  null.jaccreplcontr.bout <- null.jaccreplcontr.bout %>%
    bind_cols(null.jacc.i)
}

null.jaccreplcontr.bout$bout_contr_jaccreplses <- NA
null.jaccreplcontr.bout$bout_contr_jaccreplrank <- NA
null.jaccreplcontr.bout$bout_contr_jaccreplpval <- NA
for (i in 1:nrow(null.jaccreplcontr.bout)){
  null.jaccreplcontr.bout[i,1009] <- ((as.numeric(null.jaccreplcontr.bout[i,8]) - mean(as.numeric(null.jaccreplcontr.bout[i,9:1008])))/sd(as.numeric(null.jaccreplcontr.bout[i,9:1008])))
  rank.sim <- apply(null.jaccreplcontr.bout[i,8:1008],1,rank)[,1]
  null.jaccreplcontr.bout[i,1010] <- rank.sim[1]
  null.jaccreplcontr.bout[i,1011] <- null.jaccreplcontr.bout[i,1010]/1001
}

# Contribution of richness to Jaccard dissimilarity
null.jaccrichcontr.bout <- obs.bout %>%
  select(domainID, siteID, lat, lon, year, bout_from, bout_to, bout_contr_jaccrich)
for (i in 1:1000){
  null.jacc.i <- null.bout[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(bout_contr_jaccrich)
  null.jaccrichcontr.bout <- null.jaccrichcontr.bout %>%
    bind_cols(null.jacc.i)
}

null.jaccrichcontr.bout$bout_contr_jaccrichses <- NA
null.jaccrichcontr.bout$bout_contr_jaccrichrank <- NA
null.jaccrichcontr.bout$bout_contr_jaccrichpval <- NA
for (i in 1:nrow(null.jaccrichcontr.bout)){
  null.jaccrichcontr.bout[i,1009] <- ((as.numeric(null.jaccrichcontr.bout[i,8]) - mean(as.numeric(null.jaccrichcontr.bout[i,9:1008])))/sd(as.numeric(null.jaccrichcontr.bout[i,9:1008])))
  rank.sim <- apply(null.jaccrichcontr.bout[i,8:1008],1,rank)[,1]
  null.jaccrichcontr.bout[i,1010] <- rank.sim[1]
  null.jaccrichcontr.bout[i,1011] <- null.jaccrichcontr.bout[i,1010]/1001
}



null.jaccrepl.bout <- null.jaccrepl.bout %>%
  select(bout_jaccrepl, bout_jaccreplses, bout_jaccreplrank, bout_jaccreplpval) 
null.jaccrich.bout <- null.jaccrich.bout %>%
  select(bout_jaccrich, bout_jaccrichses, bout_jaccrichrank, bout_jaccrichpval) 
null.jaccreplcontr.bout <- null.jaccreplcontr.bout %>%
  select(bout_contr_jaccrepl, bout_contr_jaccreplses, bout_contr_jaccreplrank, bout_contr_jaccreplpval) 
null.jaccrichcontr.bout <- null.jaccrichcontr.bout %>%
  select(bout_contr_jaccrich, bout_contr_jaccrichses, bout_contr_jaccrichrank, bout_contr_jaccrichpval) 

mam.ses.bout <- null.jacc.bout %>%
  select(domainID, siteID, lat, lon, year, bout_from, bout_to, bout_jacc, bout_jaccses, bout_jaccrank, bout_jaccpval) %>%
  bind_cols(null.jaccrepl.bout) %>%
  bind_cols(null.jaccrich.bout) %>%
  bind_cols(null.jaccreplcontr.bout) %>%
  bind_cols(null.jaccrichcontr.bout)
mam.ses.bout <- as_tibble(mam.ses.bout)

saveRDS(mam.ses.bout, file="data/mammal_outBAT-FD-SES-bout_null.rds")


### year-level
# Jaccard dissimilarity
null.jacc.year <- obs.year %>%
  select(domainID, siteID, lat, lon, year_from, year_to, year_jacc)
for (i in 1:1000){
  null.jacc.i <- null.year[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(year_jacc)
  null.jacc.year <- null.jacc.year %>%
    bind_cols(null.jacc.i)
}

null.jacc.year$year_jaccses <- NA
null.jacc.year$year_jaccrank <- NA
null.jacc.year$year_jaccpval <- NA
for (i in 1:nrow(null.jacc.year)){
  null.jacc.year[i,1008] <- ((as.numeric(null.jacc.year[i,7]) - mean(as.numeric(null.jacc.year[i,8:1007])))/sd(as.numeric(null.jacc.year[i,8:1007])))
  rank.sim <- apply(null.jacc.year[i,7:1007],1,rank)[,1]
  null.jacc.year[i,1009] <- rank.sim[1]
  null.jacc.year[i,1010] <- null.jacc.year[i,1009]/1001
}

# Jaccard dissimilarity due to replacement 
null.jaccrepl.year <- obs.year %>%
  select(domainID, siteID, lat, lon, year_from, year_to, year_jaccrepl)
for (i in 1:1000){
  null.jacc.i <- null.year[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(year_jaccrepl)
  null.jaccrepl.year <- null.jaccrepl.year %>%
    bind_cols(null.jacc.i)
}

null.jaccrepl.year$year_jaccreplses <- NA
null.jaccrepl.year$year_jaccreplrank <- NA
null.jaccrepl.year$year_jaccreplpval <- NA
for (i in 1:nrow(null.jaccrepl.year)){
  null.jaccrepl.year[i,1008] <- ((as.numeric(null.jaccrepl.year[i,7]) - mean(as.numeric(null.jaccrepl.year[i,8:1007])))/sd(as.numeric(null.jaccrepl.year[i,8:1007])))
  rank.sim <- apply(null.jaccrepl.year[i,7:1007],1,rank)[,1]
  null.jaccrepl.year[i,1009] <- rank.sim[1]
  null.jaccrepl.year[i,1010] <- null.jaccrepl.year[i,1009]/1001
}

# Jaccard dissimilarity due to richness 
null.jaccrich.year <- obs.year %>%
  select(domainID, siteID, lat, lon, year_from, year_to, year_jaccrich)
for (i in 1:1000){
  null.jacc.i <- null.year[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(year_jaccrich)
  null.jaccrich.year <- null.jaccrich.year %>%
    bind_cols(null.jacc.i)
}

null.jaccrich.year$year_jaccrichses <- NA
null.jaccrich.year$year_jaccrichrank <- NA
null.jaccrich.year$year_jaccrichpval <- NA
for (i in 1:nrow(null.jaccrich.year)){
  null.jaccrich.year[i,1008] <- ((as.numeric(null.jaccrich.year[i,7]) - mean(as.numeric(null.jaccrich.year[i,8:1007])))/sd(as.numeric(null.jaccrich.year[i,8:1007])))
  rank.sim <- apply(null.jaccrich.year[i,7:1007],1,rank)[,1]
  null.jaccrich.year[i,1009] <- rank.sim[1]
  null.jaccrich.year[i,1010] <- null.jaccrich.year[i,1009]/1001
}


# Contribution of replacement to Jaccard dissimilarity
null.jaccreplcontr.year <- obs.year %>%
  select(domainID, siteID, lat, lon, year_from, year_to, year_contr_jaccrepl)
for (i in 1:1000){
  null.jacc.i <- null.year[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(year_contr_jaccrepl)
  null.jaccreplcontr.year <- null.jaccreplcontr.year %>%
    bind_cols(null.jacc.i)
}

null.jaccreplcontr.year$year_contr_jaccreplses <- NA
null.jaccreplcontr.year$year_contr_jaccreplrank <- NA
null.jaccreplcontr.year$year_contr_jaccreplpval <- NA
for (i in 1:nrow(null.jaccreplcontr.year)){
  null.jaccreplcontr.year[i,1008] <- ((as.numeric(null.jaccreplcontr.year[i,7]) - mean(as.numeric(null.jaccreplcontr.year[i,8:1007])))/sd(as.numeric(null.jaccreplcontr.year[i,8:1007])))
  rank.sim <- apply(null.jaccreplcontr.year[i,7:1007],1,rank)[,1]
  null.jaccreplcontr.year[i,1009] <- rank.sim[1]
  null.jaccreplcontr.year[i,1010] <- null.jaccreplcontr.year[i,1009]/1001
}

# Contribution of richness to Jaccard dissimilarity
null.jaccrichcontr.year <- obs.year %>%
  select(domainID, siteID, lat, lon, year_from, year_to, year_contr_jaccrich)
for (i in 1:1000){
  null.jacc.i <- null.year[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(year_contr_jaccrich)
  null.jaccrichcontr.year <- null.jaccrichcontr.year %>%
    bind_cols(null.jacc.i)
}

null.jaccrichcontr.year$year_contr_jaccrichses <- NA
null.jaccrichcontr.year$year_contr_jaccrichrank <- NA
null.jaccrichcontr.year$year_contr_jaccrichpval <- NA
for (i in 1:nrow(null.jaccrichcontr.year)){
  null.jaccrichcontr.year[i,1008] <- ((as.numeric(null.jaccrichcontr.year[i,7]) - mean(as.numeric(null.jaccrichcontr.year[i,8:1007])))/sd(as.numeric(null.jaccrichcontr.year[i,8:1007])))
  rank.sim <- apply(null.jaccrichcontr.year[i,7:1007],1,rank)[,1]
  null.jaccrichcontr.year[i,1009] <- rank.sim[1]
  null.jaccrichcontr.year[i,1010] <- null.jaccrichcontr.year[i,1009]/1001
}


null.jaccrepl.year <- null.jaccrepl.year %>%
  select(year_jaccrepl, year_jaccreplses, year_jaccreplrank, year_jaccreplpval) 
null.jaccrich.year <- null.jaccrich.year %>%
  select(year_jaccrich, year_jaccrichses, year_jaccrichrank, year_jaccrichpval) 
null.jaccreplcontr.year <- null.jaccreplcontr.year %>%
  select(year_contr_jaccrepl, year_contr_jaccreplses, year_contr_jaccreplrank, year_contr_jaccreplpval) 
null.jaccrichcontr.year <- null.jaccrichcontr.year %>%
  select(year_contr_jaccrich, year_contr_jaccrichses, year_contr_jaccrichrank, year_contr_jaccrichpval) 

mam.ses.year <- null.jacc.year %>%
  select(domainID, siteID, lat, lon, year_from, year_to, year_jacc, year_jaccses, year_jaccrank, year_jaccpval) %>%
  bind_cols(null.jaccrepl.year) %>%
  bind_cols(null.jaccrich.year) %>%
  bind_cols(null.jaccreplcontr.year) %>%
  bind_cols(null.jaccrichcontr.year)
mam.ses.year <- as_tibble(mam.ses.year)

saveRDS(mam.ses.year, file="data/mammal_outBAT-FD-SES-year_null.rds")







################### fish
data.fish.bout <- readRDS(file="data/fish_abund-bout.rds")
data.fish.year <- readRDS(file="data/fish_abund-year.rds")
sites.id <- unique(data.fish.bout$siteID)

### Reshuffle (1000x) species abundances, separately for each domain because domain will be used as a delineation of the regional species pool
# first order by species name (column)
cols <- sort(colnames(data.fish.bout[,8:ncol(data.fish.bout)]))
cols.bout <- c(colnames(data.fish.bout[,1:7]),cols)
cols.year <- c(colnames(data.fish.year[,1:5]),cols)
data.fish.bout <- data.fish.bout[cols.bout]
data.fish.year <- data.fish.year[cols.year]

# get all species ever recorded and then get list for each domain
data_small_fish <- readRDS(file = "data/fish_processed.rds")
domain.sp <- data_small_fish %>%
  select(domainID, siteID, scientificName) %>%
  group_by(domainID, scientificName) %>%
  tally(name = "count") %>% 
  ungroup()

domains <- unique(domain.sp$domainID)

# bout-level
rand.bout <- list()

for (i in 1:1000){
  cat(i)
  rand.j <- list()
  for (j in 1:length(domains)){
    data.j <- data.fish.bout %>% 
      filter(domainID == domains[j])
    domain.j <- domain.sp %>% 
      filter(domainID == domains[j])
    
    ids <- data.j %>%
      select(domainID, siteID, lat, lon, year, month, bout)
    species <- data.j %>%
      select(-c(domainID, siteID, lat, lon, year, month, bout))
    species.dom <- species[,colnames(species) %in% domain.j$scientificName]
    species.ndom <- species[,!(colnames(species) %in% domain.j$scientificName)]
    
    colnames(species.dom) <- sample(colnames(species.dom)) #sampling randomly
    species.all <- species.ndom %>%
      bind_cols(species.dom)
    species.all <- species.all[,sort(names(species.all))]
    species.all <- ids %>%
      bind_cols(species.all)
    rand.j <- rand.j %>%
      bind_rows(species.all)
  }
  rand.bout[[i]] <- rand.j
}

saveRDS(rand.bout, file="data/fish_abund-bout_null.rds")

# year-level
rand.year <- list()

for (i in 1:1000){
  cat(i)
  rand.j <- list()
  for (j in 1:length(domains)){
    data.j <- data.fish.year %>% 
      filter(domainID == domains[j])
    domain.j <- domain.sp %>% 
      filter(domainID == domains[j])
    
    ids <- data.j %>%
      select(domainID, siteID, lat, lon, year)
    species <- data.j %>%
      select(-c(domainID, siteID, lat, lon, year))
    species.dom <- species[,colnames(species) %in% domain.j$scientificName]
    species.ndom <- species[,!(colnames(species) %in% domain.j$scientificName)]
    
    colnames(species.dom) <- sample(colnames(species.dom)) #sampling randomly
    species.all <- species.ndom %>%
      bind_cols(species.dom)
    species.all <- species.all[,sort(names(species.all))]
    species.all <- ids %>%
      bind_cols(species.all)
    rand.j <- rand.j %>%
      bind_rows(species.all)
  }
  rand.year[[i]] <- rand.j
}


saveRDS(rand.year, file="data/fish_abund-year_null.rds")


###### Null functional community composition change: Package BAT
###### Get trait data
spnames <- as_tibble(colnames(data.fish.bout[,8:ncol(data.fish.bout)]))
colnames(spnames) <- "scientificName"

# Read in Traits and taxonomy
spnames.et <- read.csv(file="traits/FishTraits_clean.csv", header=TRUE, sep=",")
#some species names are misspelled
spnames.et[69,2] <- "Notropis suttkusi"
spnames.et[79,2] <- "Oncorhynchus clarki"
spnames.et[105,2] <- "Sicydium plumieri"
spnames.et <- spnames.et %>%
  rename(scientificName = species)

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
data.fish.year <- setcolorder(data.fish.year, 
                              as.character(c("domainID","siteID","lat","lon","year", 
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

### Read in species list randomizations
rand.bout <- readRDS(file="data/fish_abund-bout_null.rds")
rand.year <- readRDS(file="data/fish_abund-year_null.rds")

### Function to retrieve the diagonal of a distance matrix
### courtesy of https://stackoverflow.com/questions/39231961/extract-diagonals-from-a-distance-matrix-in-r
### retain only diagonal because we only want boout-to-boout and year-to-year comparison not all pair-wise combinations, to avoid biases associated with length of time series
subdiag <- function (dist_obj, d) {
  if (!inherits(dist_obj, "dist")) stop("please provide a 'dist' object!")
  n <- attr(dist_obj, "Size")
  if (d < 1 || d > (n - 1)) stop(sprintf("'d' needs be between 1 and %d", n - 1L))
  j_1 <- c(0, seq.int(from = n - 1, by = -1, length = n - d - 1))
  subdiag_ind <- d + cumsum(j_1)
  dist_obj[subdiag_ind]
}

null.bout <- list()
null.year <- list()

### Intra-annual (i.e., across bouts) change in community functional composition
### Loop over randomizations
for (k in 101:1000){
  cat(k)
  data.fish.bout <- rand.bout[[k]]
  sites.id <- unique(data.fish.bout$siteID)
  
  for (i in 1:length(sites.id)) {
    mi <- data.fish.bout %>%
      filter(siteID == sites.id[i])
    domain <- unique(mi$domainID)
    lat <- mean(mi$lat)
    lon <- mean(mi$lon)
    yi <- unique(mi$year)
    
    for (j in 1:length(yi)){
      mij <- mi %>%
        filter(year == yi[j])
      
      out.intra.bat <- matrix(NA, (nrow(mij)-1), 12)
      out.intra.bat <- as.data.frame(out.intra.bat)
      out.intra.bat[,1] <- c(domain)
      out.intra.bat[,2] <- sites.id[i]
      out.intra.bat[,3] <- lat
      out.intra.bat[,4] <- lon
      out.intra.bat[,5] <- yi[j]
      betadiv <- beta(mij[,8:ncol(mij)], tree=tree.fish, abund = TRUE, raref = 0)
      bouts <- mij[,7]
      out.intra.bat[,6] <- bouts[1:(nrow(bouts)-1),] #bout 1
      out.intra.bat[,7] <- bouts[2:nrow(bouts),] #bout 2
      out.intra.bat[,8] <- subdiag(betadiv$Btotal, 1)
      out.intra.bat[,9] <- subdiag(betadiv$Brepl, 1)
      out.intra.bat[,10] <- subdiag(betadiv$Brich, 1)
      out.intra.bat <- out.intra.bat %>% 
        mutate(V11 = V9/V8, V12 = V10/V8)
      
      if (j == 1) {
        out.intra.bat.ally <- out.intra.bat
      } else {
        out.intra.bat.ally <- out.intra.bat.ally %>%
          bind_rows(out.intra.bat)
      } 
    }
    
    if (i == 1) {
      out.intra.bat.all <- out.intra.bat.ally
    } else {
      out.intra.bat.all <- out.intra.bat.all %>%
        bind_rows(out.intra.bat.ally)
    } 
  }
  
  colnames(out.intra.bat.all) <- c("domainID","siteID","lat","lon","year","bout_from","bout_to","bout_jacc","bout_jaccrepl","bout_jaccrich","bout_contr_jaccrepl","bout_contr_jaccrich")
  null.bout[[k]] <- out.intra.bat.all
}

saveRDS(null.bout, file="data/fish_outBAT-FD-bout_null.rds")

### Inter-annual (i.e., across years) change in community functional composition
### Loop over randomizations
for (k in 1:1000){
  cat(k)
  data.fish.year <- rand.year[[k]]
  sites.id <- unique(data.fish.year$siteID)
  
  for (i in 1:length(sites.id)){
    mi <- data.fish.year %>%
      filter(siteID == sites.id[i])
    domain <- unique(mi$domainID)
    lat <- mean(mi$lat)
    lon <- mean(mi$lon)
    yi <- unique(mi$year)
    
    betadiv <- beta(mi[,6:ncol(mi)], tree=tree.fish, abund = TRUE, raref = 0)
    out.inter.bat <- matrix(NA, (nrow(mi)-1), 11)
    out.inter.bat <- as.data.frame(out.inter.bat)
    out.inter.bat[,1] <- c(domain)
    out.inter.bat[,2] <- sites.id[i]
    out.inter.bat[,3] <- lat
    out.inter.bat[,4] <- lon
    out.inter.bat[,5] <- yi[1:(length(yi)-1)]
    out.inter.bat[,6] <- yi[2:length(yi)]
    out.inter.bat[,7] <- subdiag(betadiv$Btotal, 1)
    out.inter.bat[,8] <- subdiag(betadiv$Brepl, 1)
    out.inter.bat[,9] <- subdiag(betadiv$Brich, 1)
    out.inter.bat <- out.inter.bat %>% 
      mutate(V10 = V8/V7, V11 = V9/V7)
    
    if (i == 1) {
      out.inter.bat.all <- out.inter.bat
    } else {
      out.inter.bat.all <- out.inter.bat.all %>%
        bind_rows(out.inter.bat)
    }
  }
  
  colnames(out.inter.bat.all) <- c("domainID","siteID","lat","lon","year_from","year_to","year_jacc","year_jaccrepl","year_jaccrich","year_contr_jaccrepl","year_contr_jaccrich")
  null.year[[k]] <- out.inter.bat.all
}

saveRDS(null.year, file="data/fish_outBAT-FD-year_null.rds")



###### Quantify the standardized effect size (SES)
obs.bout <- readRDS(file="data/fish_outBAT-FD-bout.rds")
obs.year <- readRDS(file="data/fish_outBAT-FD-year.rds")
null.bout <- readRDS(file="data/fish_outBAT-FD-bout_null.rds")
null.year <- readRDS(file="data/fish_outBAT-FD-year_null.rds")
obs.bout <- obs.bout %>% 
  mutate(bout_contr_jaccrepl = bout_jaccrepl/bout_jacc, bout_contr_jaccrich = bout_jaccrich/bout_jacc)
obs.year <- obs.year %>% 
  mutate(year_contr_jaccrepl = year_jaccrepl/year_jacc, year_contr_jaccrich = year_jaccrich/year_jacc)

### bout-level
# Jaccard dissimilarity
null.jacc.bout <- obs.bout %>%
  select(domainID, siteID, lat, lon, year, bout_from, bout_to, bout_jacc)
for (i in 1:1000){
  null.jacc.i <- null.bout[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(bout_jacc)
  null.jacc.bout <- null.jacc.bout %>%
    bind_cols(null.jacc.i)
}

null.jacc.bout$bout_jaccses <- NA
null.jacc.bout$bout_jaccrank <- NA
null.jacc.bout$bout_jaccpval <- NA
for (i in 1:nrow(null.jacc.bout)){
  null.jacc.bout[i,1009] <- ((as.numeric(null.jacc.bout[i,8]) - mean(as.numeric(null.jacc.bout[i,9:1008])))/sd(as.numeric(null.jacc.bout[i,9:1008])))
  rank.sim <- apply(null.jacc.bout[i,8:1008],1,rank)[,1]
  null.jacc.bout[i,1010] <- rank.sim[1]
  null.jacc.bout[i,1011] <- null.jacc.bout[i,1010]/1001
}

# Jaccard dissimilarity due to replacement 
null.jaccrepl.bout <- obs.bout %>%
  select(domainID, siteID, lat, lon, year, bout_from, bout_to, bout_jaccrepl)
for (i in 1:1000){
  null.jacc.i <- null.bout[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(bout_jaccrepl)
  null.jaccrepl.bout <- null.jaccrepl.bout %>%
    bind_cols(null.jacc.i)
}

null.jaccrepl.bout$bout_jaccreplses <- NA
null.jaccrepl.bout$bout_jaccreplrank <- NA
null.jaccrepl.bout$bout_jaccreplpval <- NA
for (i in 1:nrow(null.jaccrepl.bout)){
  null.jaccrepl.bout[i,1009] <- ((as.numeric(null.jaccrepl.bout[i,8]) - mean(as.numeric(null.jaccrepl.bout[i,9:1008])))/sd(as.numeric(null.jaccrepl.bout[i,9:1008])))
  rank.sim <- apply(null.jaccrepl.bout[i,8:1008],1,rank)[,1]
  null.jaccrepl.bout[i,1010] <- rank.sim[1]
  null.jaccrepl.bout[i,1011] <- null.jaccrepl.bout[i,1010]/1001
}

# Jaccard dissimilarity due to richness 
null.jaccrich.bout <- obs.bout %>%
  select(domainID, siteID, lat, lon, year, bout_from, bout_to, bout_jaccrich)
for (i in 1:1000){
  null.jacc.i <- null.bout[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(bout_jaccrich)
  null.jaccrich.bout <- null.jaccrich.bout %>%
    bind_cols(null.jacc.i)
}

null.jaccrich.bout$bout_jaccrichses <- NA
null.jaccrich.bout$bout_jaccrichrank <- NA
null.jaccrich.bout$bout_jaccrichpval <- NA
for (i in 1:nrow(null.jaccrich.bout)){
  null.jaccrich.bout[i,1009] <- ((as.numeric(null.jaccrich.bout[i,8]) - mean(as.numeric(null.jaccrich.bout[i,9:1008])))/sd(as.numeric(null.jaccrich.bout[i,9:1008])))
  rank.sim <- apply(null.jaccrich.bout[i,8:1008],1,rank)[,1]
  null.jaccrich.bout[i,1010] <- rank.sim[1]
  null.jaccrich.bout[i,1011] <- null.jaccrich.bout[i,1010]/1001
}

# Contribution of replacement to Jaccard dissimilarity
null.jaccreplcontr.bout <- obs.bout %>%
  select(domainID, siteID, lat, lon, year, bout_from, bout_to, bout_contr_jaccrepl)
for (i in 1:1000){
  null.jacc.i <- null.bout[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(bout_contr_jaccrepl)
  null.jaccreplcontr.bout <- null.jaccreplcontr.bout %>%
    bind_cols(null.jacc.i)
}

null.jaccreplcontr.bout$bout_contr_jaccreplses <- NA
null.jaccreplcontr.bout$bout_contr_jaccreplrank <- NA
null.jaccreplcontr.bout$bout_contr_jaccreplpval <- NA
for (i in 1:nrow(null.jaccreplcontr.bout)){
  null.jaccreplcontr.bout[i,1009] <- ((as.numeric(null.jaccreplcontr.bout[i,8]) - mean(as.numeric(null.jaccreplcontr.bout[i,9:1008])))/sd(as.numeric(null.jaccreplcontr.bout[i,9:1008])))
  rank.sim <- apply(null.jaccreplcontr.bout[i,8:1008],1,rank)[,1]
  null.jaccreplcontr.bout[i,1010] <- rank.sim[1]
  null.jaccreplcontr.bout[i,1011] <- null.jaccreplcontr.bout[i,1010]/1001
}

# Contribution of richness to Jaccard dissimilarity
null.jaccrichcontr.bout <- obs.bout %>%
  select(domainID, siteID, lat, lon, year, bout_from, bout_to, bout_contr_jaccrich)
for (i in 1:1000){
  null.jacc.i <- null.bout[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(bout_contr_jaccrich)
  null.jaccrichcontr.bout <- null.jaccrichcontr.bout %>%
    bind_cols(null.jacc.i)
}

null.jaccrichcontr.bout$bout_contr_jaccrichses <- NA
null.jaccrichcontr.bout$bout_contr_jaccrichrank <- NA
null.jaccrichcontr.bout$bout_contr_jaccrichpval <- NA
for (i in 1:nrow(null.jaccrichcontr.bout)){
  null.jaccrichcontr.bout[i,1009] <- ((as.numeric(null.jaccrichcontr.bout[i,8]) - mean(as.numeric(null.jaccrichcontr.bout[i,9:1008])))/sd(as.numeric(null.jaccrichcontr.bout[i,9:1008])))
  rank.sim <- apply(null.jaccrichcontr.bout[i,8:1008],1,rank)[,1]
  null.jaccrichcontr.bout[i,1010] <- rank.sim[1]
  null.jaccrichcontr.bout[i,1011] <- null.jaccrichcontr.bout[i,1010]/1001
}



null.jaccrepl.bout <- null.jaccrepl.bout %>%
  select(bout_jaccrepl, bout_jaccreplses, bout_jaccreplrank, bout_jaccreplpval) 
null.jaccrich.bout <- null.jaccrich.bout %>%
  select(bout_jaccrich, bout_jaccrichses, bout_jaccrichrank, bout_jaccrichpval) 
null.jaccreplcontr.bout <- null.jaccreplcontr.bout %>%
  select(bout_contr_jaccrepl, bout_contr_jaccreplses, bout_contr_jaccreplrank, bout_contr_jaccreplpval) 
null.jaccrichcontr.bout <- null.jaccrichcontr.bout %>%
  select(bout_contr_jaccrich, bout_contr_jaccrichses, bout_contr_jaccrichrank, bout_contr_jaccrichpval) 

fish.ses.bout <- null.jacc.bout %>%
  select(domainID, siteID, lat, lon, year, bout_from, bout_to, bout_jacc, bout_jaccses, bout_jaccrank, bout_jaccpval) %>%
  bind_cols(null.jaccrepl.bout) %>%
  bind_cols(null.jaccrich.bout) %>%
  bind_cols(null.jaccreplcontr.bout) %>%
  bind_cols(null.jaccrichcontr.bout)
fish.ses.bout <- as_tibble(fish.ses.bout)

saveRDS(fish.ses.bout, file="data/fish_outBAT-FD-SES-bout_null.rds")


### year-level
# Jaccard dissimilarity
null.jacc.year <- obs.year %>%
  select(domainID, siteID, lat, lon, year_from, year_to, year_jacc)
for (i in 1:1000){
  null.jacc.i <- null.year[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(year_jacc)
  null.jacc.year <- null.jacc.year %>%
    bind_cols(null.jacc.i)
}

null.jacc.year$year_jaccses <- NA
null.jacc.year$year_jaccrank <- NA
null.jacc.year$year_jaccpval <- NA
for (i in 1:nrow(null.jacc.year)){
  null.jacc.year[i,1008] <- ((as.numeric(null.jacc.year[i,7]) - mean(as.numeric(null.jacc.year[i,8:1007])))/sd(as.numeric(null.jacc.year[i,8:1007])))
  rank.sim <- apply(null.jacc.year[i,7:1007],1,rank)[,1]
  null.jacc.year[i,1009] <- rank.sim[1]
  null.jacc.year[i,1010] <- null.jacc.year[i,1009]/1001
}

# Jaccard dissimilarity due to replacement 
null.jaccrepl.year <- obs.year %>%
  select(domainID, siteID, lat, lon, year_from, year_to, year_jaccrepl)
for (i in 1:1000){
  null.jacc.i <- null.year[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(year_jaccrepl)
  null.jaccrepl.year <- null.jaccrepl.year %>%
    bind_cols(null.jacc.i)
}

null.jaccrepl.year$year_jaccreplses <- NA
null.jaccrepl.year$year_jaccreplrank <- NA
null.jaccrepl.year$year_jaccreplpval <- NA
for (i in 1:nrow(null.jaccrepl.year)){
  null.jaccrepl.year[i,1008] <- ((as.numeric(null.jaccrepl.year[i,7]) - mean(as.numeric(null.jaccrepl.year[i,8:1007])))/sd(as.numeric(null.jaccrepl.year[i,8:1007])))
  rank.sim <- apply(null.jaccrepl.year[i,7:1007],1,rank)[,1]
  null.jaccrepl.year[i,1009] <- rank.sim[1]
  null.jaccrepl.year[i,1010] <- null.jaccrepl.year[i,1009]/1001
}

# Jaccard dissimilarity due to richness 
null.jaccrich.year <- obs.year %>%
  select(domainID, siteID, lat, lon, year_from, year_to, year_jaccrich)
for (i in 1:1000){
  null.jacc.i <- null.year[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(year_jaccrich)
  null.jaccrich.year <- null.jaccrich.year %>%
    bind_cols(null.jacc.i)
}

null.jaccrich.year$year_jaccrichses <- NA
null.jaccrich.year$year_jaccrichrank <- NA
null.jaccrich.year$year_jaccrichpval <- NA
for (i in 1:nrow(null.jaccrich.year)){
  null.jaccrich.year[i,1008] <- ((as.numeric(null.jaccrich.year[i,7]) - mean(as.numeric(null.jaccrich.year[i,8:1007])))/sd(as.numeric(null.jaccrich.year[i,8:1007])))
  rank.sim <- apply(null.jaccrich.year[i,7:1007],1,rank)[,1]
  null.jaccrich.year[i,1009] <- rank.sim[1]
  null.jaccrich.year[i,1010] <- null.jaccrich.year[i,1009]/1001
}


# Contribution of replacement to Jaccard dissimilarity
null.jaccreplcontr.year <- obs.year %>%
  select(domainID, siteID, lat, lon, year_from, year_to, year_contr_jaccrepl)
for (i in 1:1000){
  null.jacc.i <- null.year[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(year_contr_jaccrepl)
  null.jaccreplcontr.year <- null.jaccreplcontr.year %>%
    bind_cols(null.jacc.i)
}

null.jaccreplcontr.year$year_contr_jaccreplses <- NA
null.jaccreplcontr.year$year_contr_jaccreplrank <- NA
null.jaccreplcontr.year$year_contr_jaccreplpval <- NA
for (i in 1:nrow(null.jaccreplcontr.year)){
  null.jaccreplcontr.year[i,1008] <- ((as.numeric(null.jaccreplcontr.year[i,7]) - mean(as.numeric(null.jaccreplcontr.year[i,8:1007])))/sd(as.numeric(null.jaccreplcontr.year[i,8:1007])))
  rank.sim <- apply(null.jaccreplcontr.year[i,7:1007],1,rank)[,1]
  null.jaccreplcontr.year[i,1009] <- rank.sim[1]
  null.jaccreplcontr.year[i,1010] <- null.jaccreplcontr.year[i,1009]/1001
}

# Contribution of richness to Jaccard dissimilarity
null.jaccrichcontr.year <- obs.year %>%
  select(domainID, siteID, lat, lon, year_from, year_to, year_contr_jaccrich)
for (i in 1:1000){
  null.jacc.i <- null.year[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(year_contr_jaccrich)
  null.jaccrichcontr.year <- null.jaccrichcontr.year %>%
    bind_cols(null.jacc.i)
}

null.jaccrichcontr.year$year_contr_jaccrichses <- NA
null.jaccrichcontr.year$year_contr_jaccrichrank <- NA
null.jaccrichcontr.year$year_contr_jaccrichpval <- NA
for (i in 1:nrow(null.jaccrichcontr.year)){
  null.jaccrichcontr.year[i,1008] <- ((as.numeric(null.jaccrichcontr.year[i,7]) - mean(as.numeric(null.jaccrichcontr.year[i,8:1007])))/sd(as.numeric(null.jaccrichcontr.year[i,8:1007])))
  rank.sim <- apply(null.jaccrichcontr.year[i,7:1007],1,rank)[,1]
  null.jaccrichcontr.year[i,1009] <- rank.sim[1]
  null.jaccrichcontr.year[i,1010] <- null.jaccrichcontr.year[i,1009]/1001
}


null.jaccrepl.year <- null.jaccrepl.year %>%
  select(year_jaccrepl, year_jaccreplses, year_jaccreplrank, year_jaccreplpval) 
null.jaccrich.year <- null.jaccrich.year %>%
  select(year_jaccrich, year_jaccrichses, year_jaccrichrank, year_jaccrichpval) 
null.jaccreplcontr.year <- null.jaccreplcontr.year %>%
  select(year_contr_jaccrepl, year_contr_jaccreplses, year_contr_jaccreplrank, year_contr_jaccreplpval) 
null.jaccrichcontr.year <- null.jaccrichcontr.year %>%
  select(year_contr_jaccrich, year_contr_jaccrichses, year_contr_jaccrichrank, year_contr_jaccrichpval) 

fish.ses.year <- null.jacc.year %>%
  select(domainID, siteID, lat, lon, year_from, year_to, year_jacc, year_jaccses, year_jaccrank, year_jaccpval) %>%
  bind_cols(null.jaccrepl.year) %>%
  bind_cols(null.jaccrich.year) %>%
  bind_cols(null.jaccreplcontr.year) %>%
  bind_cols(null.jaccrichcontr.year)
fish.ses.year <- as_tibble(fish.ses.year)

saveRDS(fish.ses.year, file="data/fish_outBAT-FD-SES-year_null.rds")



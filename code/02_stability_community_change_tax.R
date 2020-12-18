#Load packages
require(rhdf5)
require(here)
require(tidyverse)
require(BAT)
require(codyn)
select <- dplyr::select


################### ECOSYSTEM STABILITY & (TAXONOMIC) COMMUNITY CHANGE USING FUNCTIONS FROM PACKAGES BAT AND CODYN

################### small mammals
data.mam.bout <- readRDS(file="data/mammal_abund-bout.rds")
data.mam.year <- readRDS(file="data/mammal_abund-year.rds")
sites.id <- unique(data.mam.bout$siteID)

###### Ecosystem stability
### Intra-annual (i.e., across bouts, within a year) stability in community abundance
for (i in 1:length(sites.id)){
  mi <- data.mam.bout %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  yi <- unique(mi$year)
  
  out.intra.stab <- matrix(NA, length(yi), 7)
  out.intra.stab <- as.data.frame(out.intra.stab)
  for (j in 1:length(yi)){
    mij <- mi %>%
      filter(year == yi[j])
    
    cv <- sd(apply(mij[,8:ncol(mij)],1,sum))/mean(apply(mij[,8:ncol(mij)],1,sum))
    stab <- 1/cv
    out.intra.stab[j,1] <- c(domain)
    out.intra.stab[j,2] <- sites.id[i]
    out.intra.stab[j,3] <- lat
    out.intra.stab[j,4] <- lon
    out.intra.stab[j,5] <- yi[j]
    out.intra.stab[j,6] <- cv
    out.intra.stab[j,7] <- stab
  } 
  
  if (i == 1) {
    out.intra.stab.all <- out.intra.stab
  } else {
    out.intra.stab.all <- out.intra.stab.all %>%
      bind_rows(out.intra.stab)
  } 
}

colnames(out.intra.stab.all) <- c("domainID","siteID","lat", "lon", "year","bout_cv","bout_stability")
saveRDS(out.intra.stab.all, file="data/mammal_outStability-bout.rds")


### Inter-annual (i.e., across years) stability in community abundance
out.inter.stab <- matrix(NA, length(sites.id), 6)
out.inter.stab <- as.data.frame(out.inter.stab)

for (i in 1:length(sites.id)){
  mi <- data.mam.year %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  
  cv <- sd(apply(mi[,6:ncol(mi)],1,sum))/mean(apply(mi[,6:ncol(mi)],1,sum))
  stab <- 1/cv
  out.inter.stab[i,1] <- c(domain)
  out.inter.stab[i,2] <- sites.id[i]
  out.inter.stab[i,3] <- lat
  out.inter.stab[i,4] <- lon
  out.inter.stab[i,5] <- cv
  out.inter.stab[i,6] <- stab
}

colnames(out.inter.stab) <- c("domainID","siteID","lat", "lon", "year_cv","year_stability")
saveRDS(out.inter.stab, file="data/mammal_outStability-year.rds")



###### Community change: Package BAT
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


### Intra-annual (i.e., across bouts) change in community composition
for (i in 1:length(sites.id)){
  mi <- data.mam.bout %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  yi <- unique(mi$year)
  
  for (j in 1:length(yi)){
    mij <- mi %>%
      filter(year == yi[j])
    
      out.intra.bat <- matrix(NA, (nrow(mij)-1), 10)
      out.intra.bat <- as.data.frame(out.intra.bat)
      out.intra.bat[,1] <- c(domain)
      out.intra.bat[,2] <- sites.id[i]
      out.intra.bat[,3] <- lat
      out.intra.bat[,4] <- lon
      out.intra.bat[,5] <- yi[j]
      betadiv <- beta(mij[,8:ncol(mij)], abund = TRUE, raref = 0)
      bouts <- mij[,7]
      out.intra.bat[,6] <- bouts[1:(nrow(bouts)-1),] #bout 1
      out.intra.bat[,7] <- bouts[2:nrow(bouts),] #bout 2
      out.intra.bat[,8] <- subdiag(betadiv$Btotal, 1)
      out.intra.bat[,9] <- subdiag(betadiv$Brepl, 1)
      out.intra.bat[,10] <- subdiag(betadiv$Brich, 1)
    
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


colnames(out.intra.bat.all) <- c("domainID","siteID","lat","lon","year","bout_from","bout_to","bout_jacc","bout_jaccrepl","bout_jaccrich")
saveRDS(out.intra.bat.all, file="data/mammal_outBAT-bout.rds")


### Inter-annual (i.e., across year) change in community composition
for (i in 1:length(sites.id)){
  mi <- data.mam.year %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  yi <- unique(mi$year)
  
  betadiv <- beta(mi[,6:ncol(mi)], abund = TRUE, raref = 0)
  out.inter.bat <- matrix(NA, (nrow(mi)-1), 9)
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
  
  if (i == 1) {
    out.inter.bat.all <- out.inter.bat
  } else {
    out.inter.bat.all <- out.inter.bat.all %>%
      bind_rows(out.inter.bat)
  }
}

colnames(out.inter.bat.all) <- c("domainID","siteID","lat","lon","year_from","year_to","year_jacc","year_jaccrepl","year_jaccrich")
saveRDS(out.inter.bat.all, file="data/mammal_outBAT-year.rds")



###### Community change: Package codyn
mean.coord <- data.mam.bout %>%
  group_by(domainID, siteID) %>%
  summarize(lat = mean(lat), lon = mean(lon, na.rm=TRUE)) 

data.mam.bout.piv <- data.mam.bout %>%
  pivot_longer(-c(domainID, siteID, lat, lon, year, month, bout), names_to = "scientificName", values_to = "abundance") 

data.mam.year.piv <- data.mam.year %>%
  pivot_longer(-c(domainID, siteID, lat, lon, year), names_to = "scientificName", values_to = "abundance") 

### Intra-annual (i.e., across bout) change in community composition
out.intra.codyn <- RAC_change(data.mam.bout.piv, 
                               time.var = "bout", species.var = "scientificName",
                               abundance.var = "abundance", replicate.var = "siteID") %>%
  #remove bout comparisons that cross years
  filter(substr(bout, 1, 4) == substr(bout2, 1, 4)) %>%
  left_join(mean.coord, by="siteID") %>%
  mutate(bout_from = bout, bout_to = bout2, year = as.numeric(substr(bout_from, start = 1, stop = 4))) %>%
  select(domainID, siteID, lat, lon, year, bout_from, bout_to, richness_change, evenness_change, rank_change, gains, losses)

saveRDS(out.intra.codyn, file="data/mammal_outCodyn-bout.rds")

### Inter-annual (i.e., across year) change in community composition
out.inter.codyn <- RAC_change(data.mam.year.piv, 
                               time.var = "year", species.var = "scientificName",
                               abundance.var = "abundance", replicate.var = "siteID") %>%
  left_join(mean.coord, by="siteID") %>%
  mutate(year_from = year, year_to = year2) %>%
  select(domainID, siteID, lat, lon, year_from, year_to, richness_change, evenness_change, rank_change, gains, losses)

saveRDS(out.inter.codyn, file="data/mammal_outCodyn-year.rds")




################### ground beetles
data.beetle.bout <- readRDS(file="data/beetle_abund-bout.rds")
data.beetle.year <- readRDS(file="data/beetle_abund-year.rds")
sites.id <- unique(data.beetle.bout$siteID)

###### Ecosystem stability
### Intra-annual (i.e., across bouts, within a year) stability in community abundance
for (i in 1:length(sites.id)){
  mi <- data.beetle.bout %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  yi <- unique(mi$year)
  
  out.intra.stab <- matrix(NA, length(yi), 7)
  out.intra.stab <- as.data.frame(out.intra.stab)
  for (j in 1:length(yi)){
    mij <- mi %>%
      filter(year == yi[j])
    
    cv <- sd(apply(mij[,7:ncol(mij)],1,sum))/mean(apply(mij[,7:ncol(mij)],1,sum))
    stab <- 1/cv
    out.intra.stab[j,1] <- c(domain)
    out.intra.stab[j,2] <- sites.id[i]
    out.intra.stab[j,3] <- lat
    out.intra.stab[j,4] <- lon
    out.intra.stab[j,5] <- yi[j]
    out.intra.stab[j,6] <- cv
    out.intra.stab[j,7] <- stab
  } 
  
  if (i == 1) {
    out.intra.stab.all <- out.intra.stab
  } else {
    out.intra.stab.all <- out.intra.stab.all %>%
      bind_rows(out.intra.stab)
  } 
}

colnames(out.intra.stab.all) <- c("domainID","siteID","lat", "lon", "year","bout_cv","bout_stability")
out.intra.stab.all <- out.intra.stab.all %>%
  mutate(year = as.numeric(year))
saveRDS(out.intra.stab.all, file="data/beetle_outStability-bout.rds")


### Inter-annual (i.e., across years) stability in community abundance
out.inter.stab <- matrix(NA, length(sites.id), 6)
out.inter.stab <- as.data.frame(out.inter.stab)

for (i in 1:length(sites.id)){
  mi <- data.beetle.year %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  
  cv <- sd(apply(mi[,6:ncol(mi)],1,sum))/mean(apply(mi[,6:ncol(mi)],1,sum))
  stab <- 1/cv
  out.inter.stab[i,1] <- c(domain)
  out.inter.stab[i,2] <- sites.id[i]
  out.inter.stab[i,3] <- lat
  out.inter.stab[i,4] <- lon
  out.inter.stab[i,5] <- cv
  out.inter.stab[i,6] <- stab
}

colnames(out.inter.stab) <- c("domainID","siteID","lat", "lon", "year_cv","year_stability")
saveRDS(out.inter.stab, file="data/beetle_outStability-year.rds")



###### Community change: Package BAT
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


### Intra-annual (i.e., across bouts) change in community composition
for (i in 1:length(sites.id)){
  mi <- data.beetle.bout %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  yi <- unique(mi$year)
  
  for (j in 1:length(yi)){
    mij <- mi %>%
      filter(year == yi[j])
    
    out.intra.bat <- matrix(NA, (nrow(mij)-1), 10)
    out.intra.bat <- as.data.frame(out.intra.bat)
    out.intra.bat[,1] <- c(domain)
    out.intra.bat[,2] <- sites.id[i]
    out.intra.bat[,3] <- lat
    out.intra.bat[,4] <- lon
    out.intra.bat[,5] <- yi[j]
    betadiv <- beta(mij[,7:ncol(mij)], abund = TRUE, raref = 0)
    bouts <- mij[,6]
    out.intra.bat[,6] <- bouts[1:(nrow(bouts)-1),] #bout 1
    out.intra.bat[,7] <- bouts[2:nrow(bouts),] #bout 2
    out.intra.bat[,8] <- subdiag(betadiv$Btotal, 1)
    out.intra.bat[,9] <- subdiag(betadiv$Brepl, 1)
    out.intra.bat[,10] <- subdiag(betadiv$Brich, 1)
    
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


colnames(out.intra.bat.all) <- c("domainID","siteID","lat","lon","year","bout_from","bout_to","bout_jacc","bout_jaccrepl","bout_jaccrich")
out.intra.bat.all <- out.intra.bat.all %>%
  mutate(year = as.numeric(year))
saveRDS(out.intra.bat.all, file="data/beetle_outBAT-bout.rds")


### Inter-annual (i.e., across year) change in community composition
for (i in 1:length(sites.id)){
  mi <- data.beetle.year %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  yi <- unique(mi$year)
  
  betadiv <- beta(mi[,6:ncol(mi)], abund = TRUE, raref = 0)
  out.inter.bat <- matrix(NA, (nrow(mi)-1), 9)
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
  
  if (i == 1) {
    out.inter.bat.all <- out.inter.bat
  } else {
    out.inter.bat.all <- out.inter.bat.all %>%
      bind_rows(out.inter.bat)
  }
}

colnames(out.inter.bat.all) <- c("domainID","siteID","lat","lon","year_from","year_to","year_jacc","year_jaccrepl","year_jaccrich")
saveRDS(out.inter.bat.all, file="data/beetle_outBAT-year.rds")



###### Community change: Package codyn
mean.coord <- data.beetle.bout %>%
  group_by(domainID, siteID) %>%
  summarize(lat = mean(lat), lon = mean(lon, na.rm=TRUE)) 

data.beetle.bout.piv <- data.beetle.bout %>%
  pivot_longer(-c(domainID, siteID, lat, lon, year, bout), names_to = "scientificName", values_to = "abundance") 

data.beetle.year.piv <- data.beetle.year %>%
  pivot_longer(-c(domainID, siteID, lat, lon, year), names_to = "scientificName", values_to = "abundance") 

### Intra-annual (i.e., across bout) change in community composition
out.intra.codyn <- RAC_change(data.beetle.bout.piv, 
                              time.var = "bout", species.var = "scientificName",
                              abundance.var = "abundance", replicate.var = "siteID") %>%
  #remove bout comparisons that cross years
  filter(substr(bout, 6, 9) == substr(bout2, 6, 9)) %>%
  left_join(mean.coord, by="siteID") %>%
  mutate(bout_from = bout, bout_to = bout2, year = as.numeric(substr(bout_from, start = 6, stop = 9))) %>%
  select(domainID, siteID, lat, lon, year, bout_from, bout_to, richness_change, evenness_change, rank_change, gains, losses)

saveRDS(out.intra.codyn, file="data/beetle_outCodyn-bout.rds")

### Inter-annual (i.e., across year) change in community composition
out.inter.codyn <- RAC_change(data.beetle.year.piv, 
                              time.var = "year", species.var = "scientificName",
                              abundance.var = "abundance", replicate.var = "siteID") %>%
  left_join(mean.coord, by="siteID") %>%
  mutate(year_from = year, year_to = year2) %>%
  select(domainID, siteID, lat, lon, year_from, year_to, richness_change, evenness_change, rank_change, gains, losses)

saveRDS(out.inter.codyn, file="data/beetle_outCodyn-year.rds")




################### fish
data.fish.bout <- readRDS(file="data/fish_abund-bout.rds")
data.fish.year <- readRDS(file="data/fish_abund-year.rds")
sites.id <- unique(data.fish.bout$siteID)

###### Ecosystem stability
### Intra-annual (i.e., across bouts, within a year) stability in community abundance
for (i in 1:length(sites.id)){
  mi <- data.fish.bout %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  yi <- unique(mi$year)
  
  out.intra.stab <- matrix(NA, length(yi), 7)
  out.intra.stab <- as.data.frame(out.intra.stab)
  for (j in 1:length(yi)){
    mij <- mi %>%
      filter(year == yi[j])
    
    cv <- sd(apply(mij[,8:ncol(mij)],1,sum))/mean(apply(mij[,8:ncol(mij)],1,sum))
    stab <- 1/cv
    out.intra.stab[j,1] <- c(domain)
    out.intra.stab[j,2] <- sites.id[i]
    out.intra.stab[j,3] <- lat
    out.intra.stab[j,4] <- lon
    out.intra.stab[j,5] <- yi[j]
    out.intra.stab[j,6] <- cv
    out.intra.stab[j,7] <- stab
  } 
  
  if (i == 1) {
    out.intra.stab.all <- out.intra.stab
  } else {
    out.intra.stab.all <- out.intra.stab.all %>%
      bind_rows(out.intra.stab)
  } 
}

colnames(out.intra.stab.all) <- c("domainID","siteID","lat", "lon", "year","bout_cv","bout_stability")
out.intra.stab.all <- out.intra.stab.all %>%
  mutate(year = as.numeric(year))
saveRDS(out.intra.stab.all, file="data/fish_outStability-bout.rds")


### Inter-annual (i.e., across years) stability in community abundance
out.inter.stab <- matrix(NA, length(sites.id), 6)
out.inter.stab <- as.data.frame(out.inter.stab)

for (i in 1:length(sites.id)){
  mi <- data.fish.year %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  
  cv <- sd(apply(mi[,6:ncol(mi)],1,sum))/mean(apply(mi[,6:ncol(mi)],1,sum))
  stab <- 1/cv
  out.inter.stab[i,1] <- c(domain)
  out.inter.stab[i,2] <- sites.id[i]
  out.inter.stab[i,3] <- lat
  out.inter.stab[i,4] <- lon
  out.inter.stab[i,5] <- cv
  out.inter.stab[i,6] <- stab
}

colnames(out.inter.stab) <- c("domainID","siteID","lat", "lon", "year_cv","year_stability")
saveRDS(out.inter.stab, file="data/fish_outStability-year.rds")



###### Community change: Package BAT
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


### Intra-annual (i.e., across bouts) change in community composition
for (i in 1:length(sites.id)){
  mi <- data.fish.bout %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  yi <- unique(mi$year)
  
  for (j in 1:length(yi)){
    mij <- mi %>%
      filter(year == yi[j])
    
    out.intra.bat <- matrix(NA, (nrow(mij)-1), 10)
    out.intra.bat <- as.data.frame(out.intra.bat)
    out.intra.bat[,1] <- c(domain)
    out.intra.bat[,2] <- sites.id[i]
    out.intra.bat[,3] <- lat
    out.intra.bat[,4] <- lon
    out.intra.bat[,5] <- yi[j]
    betadiv <- beta(mij[,8:ncol(mij)], abund = TRUE, raref = 0)
    bouts <- mij[,7]
    out.intra.bat[,6] <- bouts[1:(nrow(bouts)-1),] #bout 1
    out.intra.bat[,7] <- bouts[2:nrow(bouts),] #bout 2
    out.intra.bat[,8] <- subdiag(betadiv$Btotal, 1)
    out.intra.bat[,9] <- subdiag(betadiv$Brepl, 1)
    out.intra.bat[,10] <- subdiag(betadiv$Brich, 1)
    
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


colnames(out.intra.bat.all) <- c("domainID","siteID","lat","lon","year","bout_from","bout_to","bout_jacc","bout_jaccrepl","bout_jaccrich")
out.intra.bat.all <- out.intra.bat.all %>%
  mutate(year = as.numeric(year))
saveRDS(out.intra.bat.all, file="data/fish_outBAT-bout.rds")


### Inter-annual (i.e., across year) change in community composition
for (i in 1:length(sites.id)){
  mi <- data.fish.year %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  yi <- unique(mi$year)
  
  betadiv <- beta(mi[,6:ncol(mi)], abund = TRUE, raref = 0)
  out.inter.bat <- matrix(NA, (nrow(mi)-1), 9)
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
  
  if (i == 1) {
    out.inter.bat.all <- out.inter.bat
  } else {
    out.inter.bat.all <- out.inter.bat.all %>%
      bind_rows(out.inter.bat)
  }
}

colnames(out.inter.bat.all) <- c("domainID","siteID","lat","lon","year_from","year_to","year_jacc","year_jaccrepl","year_jaccrich")
saveRDS(out.inter.bat.all, file="data/fish_outBAT-year.rds")



###### Community change: Package codyn
mean.coord <- data.fish.bout %>%
  group_by(domainID, siteID) %>%
  summarize(lat = mean(lat), lon = mean(lon, na.rm=TRUE)) 

data.fish.bout.piv <- data.fish.bout %>%
  pivot_longer(-c(domainID, siteID, lat, lon, year, month, bout), names_to = "scientificName", values_to = "abundance") 

data.fish.year.piv <- data.fish.year %>%
  pivot_longer(-c(domainID, siteID, lat, lon, year), names_to = "scientificName", values_to = "abundance") 

### Intra-annual (i.e., across bout) change in community composition
out.intra.codyn <- RAC_change(data.fish.bout.piv, 
                              time.var = "bout", species.var = "scientificName",
                              abundance.var = "abundance", replicate.var = "siteID") %>%
  #remove bout comparisons that cross years
  filter(substr(bout, 1, 4) == substr(bout2, 1, 4)) %>%
  left_join(mean.coord, by="siteID") %>%
  mutate(bout_from = bout, bout_to = bout2, year = as.numeric(substr(bout_from, start = 1, stop = 4))) %>%
  select(domainID, siteID, lat, lon, year, bout_from, bout_to, richness_change, evenness_change, rank_change, gains, losses)

saveRDS(out.intra.codyn, file="data/fish_outCodyn-bout.rds")

### Inter-annual (i.e., across year) change in community composition
out.inter.codyn <- RAC_change(data.fish.year.piv, 
                              time.var = "year", species.var = "scientificName",
                              abundance.var = "abundance", replicate.var = "siteID") %>%
  left_join(mean.coord, by="siteID") %>%
  mutate(year_from = year, year_to = year2) %>%
  select(domainID, siteID, lat, lon, year_from, year_to, richness_change, evenness_change, rank_change, gains, losses)

saveRDS(out.inter.codyn, file="data/fish_outCodyn-year.rds")






################### aquatic macroinvertebrates
data.macroinv.bout <- readRDS(file="data/macroinv_abund-bout.rds")
data.macroinv.year <- readRDS(file="data/macroinv_abund-year.rds")
sites.id <- unique(data.macroinv.bout$siteID)

###### Ecosystem stability
### Intra-annual (i.e., across bouts, within a year) stability in community abundance
for (i in 1:length(sites.id)){
  mi <- data.macroinv.bout %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  yi <- unique(mi$year)
  
  out.intra.stab <- matrix(NA, length(yi), 7)
  out.intra.stab <- as.data.frame(out.intra.stab)
  for (j in 1:length(yi)){
    mij <- mi %>%
      filter(year == yi[j])
    
    cv <- sd(apply(mij[,8:ncol(mij)],1,sum))/mean(apply(mij[,8:ncol(mij)],1,sum))
    stab <- 1/cv
    out.intra.stab[j,1] <- c(domain)
    out.intra.stab[j,2] <- sites.id[i]
    out.intra.stab[j,3] <- lat
    out.intra.stab[j,4] <- lon
    out.intra.stab[j,5] <- yi[j]
    out.intra.stab[j,6] <- cv
    out.intra.stab[j,7] <- stab
  } 
  
  if (i == 1) {
    out.intra.stab.all <- out.intra.stab
  } else {
    out.intra.stab.all <- out.intra.stab.all %>%
      bind_rows(out.intra.stab)
  } 
}

colnames(out.intra.stab.all) <- c("domainID","siteID","lat", "lon", "year","bout_cv","bout_stability")
out.intra.stab.all <- out.intra.stab.all %>%
  mutate(year = as.numeric(year))
saveRDS(out.intra.stab.all, file="data/macroinv_outStability-bout.rds")


### Inter-annual (i.e., across years) stability in community abundance
out.inter.stab <- matrix(NA, length(sites.id), 6)
out.inter.stab <- as.data.frame(out.inter.stab)

for (i in 1:length(sites.id)){
  mi <- data.macroinv.year %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  
  cv <- sd(apply(mi[,6:ncol(mi)],1,sum))/mean(apply(mi[,6:ncol(mi)],1,sum))
  stab <- 1/cv
  out.inter.stab[i,1] <- c(domain)
  out.inter.stab[i,2] <- sites.id[i]
  out.inter.stab[i,3] <- lat
  out.inter.stab[i,4] <- lon
  out.inter.stab[i,5] <- cv
  out.inter.stab[i,6] <- stab
}

colnames(out.inter.stab) <- c("domainID","siteID","lat", "lon", "year_cv","year_stability")
saveRDS(out.inter.stab, file="data/macroinv_outStability-year.rds")



###### Community change: Package BAT
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


### Intra-annual (i.e., across bouts) change in community composition
for (i in 1:length(sites.id)){
  mi <- data.macroinv.bout %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  yi <- unique(mi$year)
  
  for (j in 1:length(yi)){
    mij <- mi %>%
      filter(year == yi[j])
    
    out.intra.bat <- matrix(NA, (nrow(mij)-1), 10)
    out.intra.bat <- as.data.frame(out.intra.bat)
    out.intra.bat[,1] <- c(domain)
    out.intra.bat[,2] <- sites.id[i]
    out.intra.bat[,3] <- lat
    out.intra.bat[,4] <- lon
    out.intra.bat[,5] <- yi[j]
    betadiv <- beta(mij[,8:ncol(mij)], abund = TRUE, raref = 0)
    bouts <- mij[,7]
    out.intra.bat[,6] <- bouts[1:(nrow(bouts)-1),] #bout 1
    out.intra.bat[,7] <- bouts[2:nrow(bouts),] #bout 2
    out.intra.bat[,8] <- subdiag(betadiv$Btotal, 1)
    out.intra.bat[,9] <- subdiag(betadiv$Brepl, 1)
    out.intra.bat[,10] <- subdiag(betadiv$Brich, 1)
    
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


colnames(out.intra.bat.all) <- c("domainID","siteID","lat","lon","year","bout_from","bout_to","bout_jacc","bout_jaccrepl","bout_jaccrich")
out.intra.bat.all <- out.intra.bat.all %>%
  mutate(year = as.numeric(year))
saveRDS(out.intra.bat.all, file="data/macroinv_outBAT-bout.rds")


### Inter-annual (i.e., across year) change in community composition
for (i in 1:length(sites.id)){
  mi <- data.macroinv.year %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  yi <- unique(mi$year)
  
  betadiv <- beta(mi[,6:ncol(mi)], abund = TRUE, raref = 0)
  out.inter.bat <- matrix(NA, (nrow(mi)-1), 9)
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
  
  if (i == 1) {
    out.inter.bat.all <- out.inter.bat
  } else {
    out.inter.bat.all <- out.inter.bat.all %>%
      bind_rows(out.inter.bat)
  }
}

colnames(out.inter.bat.all) <- c("domainID","siteID","lat","lon","year_from","year_to","year_jacc","year_jaccrepl","year_jaccrich")
saveRDS(out.inter.bat.all, file="data/macroinv_outBAT-year.rds")



###### Community change: Package codyn
mean.coord <- data.macroinv.bout %>%
  group_by(domainID, siteID) %>%
  summarize(lat = mean(lat), lon = mean(lon, na.rm=TRUE)) 

data.macroinv.bout.piv <- data.macroinv.bout %>%
  pivot_longer(-c(domainID, siteID, lat, lon, year, month, bout), names_to = "scientificName", values_to = "abundance") 

data.macroinv.year.piv <- data.macroinv.year %>%
  pivot_longer(-c(domainID, siteID, lat, lon, year), names_to = "scientificName", values_to = "abundance") 

### Intra-annual (i.e., across bout) change in community composition
out.intra.codyn <- RAC_change(data.macroinv.bout.piv, 
                              time.var = "bout", species.var = "scientificName",
                              abundance.var = "abundance", replicate.var = "siteID") %>%
  #remove bout comparisons that cross years
  filter(substr(bout, 1, 4) == substr(bout2, 1, 4)) %>%
  left_join(mean.coord, by="siteID") %>%
  mutate(bout_from = bout, bout_to = bout2, year = as.numeric(substr(bout_from, start = 1, stop = 4))) %>%
  select(domainID, siteID, lat, lon, year, bout_from, bout_to, richness_change, evenness_change, rank_change, gains, losses)

saveRDS(out.intra.codyn, file="data/macroinv_outCodyn-bout.rds")

### Inter-annual (i.e., across year) change in community composition
out.inter.codyn <- RAC_change(data.macroinv.year.piv, 
                              time.var = "year", species.var = "scientificName",
                              abundance.var = "abundance", replicate.var = "siteID") %>%
  left_join(mean.coord, by="siteID") %>%
  mutate(year_from = year, year_to = year2) %>%
  select(domainID, siteID, lat, lon, year_from, year_to, richness_change, evenness_change, rank_change, gains, losses)

saveRDS(out.inter.codyn, file="data/macroinv_outCodyn-year.rds")









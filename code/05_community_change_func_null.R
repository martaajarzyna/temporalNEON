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
sites.id <- unique(data.mam.bout$siteID)

### Reshuffle (1000x) species abundances, separately for each domain because domain will be used as a delineation of the regional species pool
# first order by species name (column)
cols <- sort(colnames(data.mam.bout[,8:ncol(data.mam.bout)]))
cols.bout <- c(colnames(data.mam.bout[,1:7]),cols)
data.mam.bout <- data.mam.bout[cols.bout]

# get all species ever recorded and then get list for each domain
data_small_mammal <- readRDS(file = "data/mammal_processed.rds")
domain.sp <- data_small_mammal %>%
  select(domainID, siteID, scientificName) %>%
  group_by(domainID, scientificName) %>%
  tally(name = "count") %>% 
  ungroup()

domains <- unique(domain.sp$domainID)


null.bout <- list()

system.time(
  for (k in 1:1000){
    cat(k)
    data.mam.bout <- rand.bout[[k]]
    sites.id <- unique(data.mam.bout$siteID)
    
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
    null.bout[[k]] <- out.bout.bat.all
  }
)

saveRDS(null.bout, file="data/mammal_outBAT-FD-bout_null-consecutive-all.rds")


###### Quantify the standardized effect size (SES)
obs.bout.tax <- readRDS(file="data/mammal_outBAT-bout-consecutive-all.rds")
obs.bout <- readRDS(file="data/mammal_outBAT-FD-bout-consecutive-all.rds")
null.bout <- readRDS(file="data/mammal_outBAT-FD-bout_null-consecutive-all.rds")

# Jaccard dissimilarity
null.jacc.bout <- obs.bout %>%
  select(domainID, siteID, lat, lon, bout_lag, bout_consec_jacc)
for (i in 1:1000){
  null.jacc.i <- null.bout[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(bout_consec_jacc)
  null.jacc.bout <- null.jacc.bout %>%
    bind_cols(null.jacc.i)
}

null.jacc.bout$bout_consec_jaccses <- NA
null.jacc.bout$bout_consec_jaccrank <- NA
null.jacc.bout$bout_consec_jaccpval <- NA
for (i in 1:nrow(null.jacc.bout)){
  null.jacc.bout[i,1007] <- ((as.numeric(null.jacc.bout[i,6]) - mean(as.numeric(null.jacc.bout[i,7:1006])))/sd(as.numeric(null.jacc.bout[i,7:1006])))
  rank.sim <- apply(null.jacc.bout[i,7:1006],1,rank)[,1]
  null.jacc.bout[i,1008] <- rank.sim[1]
  null.jacc.bout[i,1009] <- null.jacc.bout[i,1008]/1001
}

# Jaccard dissimilarity due to replacement 
null.jaccrepl.bout <- obs.bout %>%
  select(domainID, siteID, lat, lon, bout_lag, bout_consec_jaccrepl)
for (i in 1:1000){
  null.jaccrepl.i <- null.bout[[i]]
  null.jaccrepl.i <- null.jaccrepl.i %>%
    select(bout_consec_jaccrepl)
  null.jaccrepl.bout <- null.jaccrepl.bout %>%
    bind_cols(null.jaccrepl.i)
}

null.jaccrepl.bout$bout_consec_jaccreplses <- NA
null.jaccrepl.bout$bout_consec_jaccreplrank <- NA
null.jaccrepl.bout$bout_consec_jaccreplpval <- NA
for (i in 1:nrow(null.jaccrepl.bout)){
  null.jaccrepl.bout[i,1007] <- ((as.numeric(null.jaccrepl.bout[i,6]) - mean(as.numeric(null.jaccrepl.bout[i,7:1006])))/sd(as.numeric(null.jaccrepl.bout[i,7:1006])))
  rank.sim <- apply(null.jaccrepl.bout[i,7:1006],1,rank)[,1]
  null.jaccrepl.bout[i,1008] <- rank.sim[1]
  null.jaccrepl.bout[i,1009] <- null.jaccrepl.bout[i,1008]/1001
}

# Jaccard dissimilarity due to richness 
null.jaccrich.bout <- obs.bout %>%
  select(domainID, siteID, lat, lon, bout_lag, bout_consec_jaccrich)
for (i in 1:1000){
  null.jaccrich.i <- null.bout[[i]]
  null.jaccrich.i <- null.jaccrich.i %>%
    select(bout_consec_jaccrich)
  null.jaccrich.bout <- null.jaccrich.bout %>%
    bind_cols(null.jaccrich.i)
}

null.jaccrich.bout$bout_consec_jaccrichses <- NA
null.jaccrich.bout$bout_consec_jaccrichrank <- NA
null.jaccrich.bout$bout_consec_jaccrichpval <- NA
for (i in 1:nrow(null.jaccrich.bout)){
  null.jaccrich.bout[i,1007] <- ((as.numeric(null.jaccrich.bout[i,6]) - mean(as.numeric(null.jaccrich.bout[i,7:1006])))/sd(as.numeric(null.jaccrich.bout[i,7:1006])))
  rank.sim <- apply(null.jaccrich.bout[i,7:1006],1,rank)[,1]
  null.jaccrich.bout[i,1008] <- rank.sim[1]
  null.jaccrich.bout[i,1009] <- null.jaccrich.bout[i,1008]/1001
}

# Contribution of replacement to Jaccard dissimilarity
null.jaccrepl.contr.bout <- obs.bout %>%
  select(domainID, siteID, lat, lon, bout_lag, bout_consec_contr_jaccrepl)
for (i in 1:1000){
  null.jaccrepl.i <- null.bout[[i]]
  null.jaccrepl.i <- null.jaccrepl.i %>%
    select(bout_consec_contr_jaccrepl)
  null.jaccrepl.contr.bout <- null.jaccrepl.contr.bout %>%
    bind_cols(null.jaccrepl.i)
}

null.jaccrepl.contr.bout$bout_consec_contr_jaccreplses <- NA
null.jaccrepl.contr.bout$bout_consec_contr_jaccreplrank <- NA
null.jaccrepl.contr.bout$bout_consec_contr_jaccreplpval <- NA
for (i in 1:nrow(null.jaccrepl.bout)){
  null.jaccrepl.contr.bout[i,1007] <- ((as.numeric(null.jaccrepl.contr.bout[i,6]) - mean(as.numeric(null.jaccrepl.contr.bout[i,7:1006])))/sd(as.numeric(null.jaccrepl.contr.bout[i,7:1006])))
  rank.sim <- apply(null.jaccrepl.contr.bout[i,7:1006],1,rank)[,1]
  null.jaccrepl.contr.bout[i,1008] <- rank.sim[1]
  null.jaccrepl.contr.bout[i,1009] <- null.jaccrepl.contr.bout[i,1008]/1001
}


# Contribution of richness to Jaccard dissimilarity
null.jaccrich.contr.bout <- obs.bout %>%
  select(domainID, siteID, lat, lon, bout_lag, bout_consec_contr_jaccrich)
for (i in 1:1000){
  null.jaccrich.i <- null.bout[[i]]
  null.jaccrich.i <- null.jaccrich.i %>%
    select(bout_consec_contr_jaccrich)
  null.jaccrich.contr.bout <- null.jaccrich.contr.bout %>%
    bind_cols(null.jaccrich.i)
}

null.jaccrich.contr.bout$bout_consec_contr_jaccrichses <- NA
null.jaccrich.contr.bout$bout_consec_contr_jaccrichrank <- NA
null.jaccrich.contr.bout$bout_consec_contr_jaccrichpval <- NA
for (i in 1:nrow(null.jaccrich.bout)){
  null.jaccrich.contr.bout[i,1007] <- ((as.numeric(null.jaccrich.contr.bout[i,6]) - mean(as.numeric(null.jaccrich.contr.bout[i,7:1006])))/sd(as.numeric(null.jaccrich.contr.bout[i,7:1006])))
  rank.sim <- apply(null.jaccrich.contr.bout[i,7:1006],1,rank)[,1]
  null.jaccrich.contr.bout[i,1008] <- rank.sim[1]
  null.jaccrich.contr.bout[i,1009] <- null.jaccrich.contr.bout[i,1008]/1001
}

##subset and bind
null.jaccrepl.bout <- null.jaccrepl.bout %>%
  select(bout_consec_jaccrepl, bout_consec_jaccreplses, bout_consec_jaccreplrank, bout_consec_jaccreplpval) 
null.jaccrich.bout <- null.jaccrich.bout %>%
  select(bout_consec_jaccrich, bout_consec_jaccrichses, bout_consec_jaccrichrank, bout_consec_jaccrichpval) 
null.jaccrepl.contr.bout <- null.jaccrepl.contr.bout %>%
  select(bout_consec_contr_jaccrepl, bout_consec_contr_jaccreplses, bout_consec_contr_jaccreplrank, bout_consec_contr_jaccreplpval) 
null.jaccrich.contr.bout <- null.jaccrich.contr.bout %>%
  select(bout_consec_contr_jaccrich, bout_consec_contr_jaccrichses, bout_consec_contr_jaccrichrank, bout_consec_contr_jaccrichpval) 

mam.ses.bout <- null.jacc.bout %>%
  select(domainID, siteID, lat, lon, bout_lag, bout_consec_jacc, bout_consec_jaccses, bout_consec_jaccrank, bout_consec_jaccpval) %>%
  bind_cols(null.jaccrepl.bout) %>%
  bind_cols(null.jaccrich.bout) %>%
  bind_cols(null.jaccrepl.contr.bout) %>%
  bind_cols(null.jaccrich.contr.bout)
mam.ses.bout <- as_tibble(mam.ses.bout)

saveRDS(mam.ses.bout, file="data/mammal_outBAT-FD-SES-bout_null-consecutive-all.rds")



########read in all data 
stab.bout.cons.all <- readRDS(file="data/mammal_outStability-bout-consecutive-all.rds")
stab.bout.cons.all$bout_lag <- as.numeric(stab.bout.cons.all$bout_lag)
bat.bout.cons.all <- readRDS(file="data/mammal_outBAT-FD-SES-bout_null-consecutive-all.rds")
bat.bout.cons.all$bout_lag <- as.numeric(bat.bout.cons.all$bout_lag)

all.bout.cons.all <- cbind(stab.bout.cons.all, bat.bout.cons.all[,6:ncol(bat.bout.cons.all)])
saveRDS(all.bout.cons.all, file="data/mammal_outAll-FD-SES-bout_null-consecutive-all.rds")





################### fish
data.fish.bout <- readRDS(file="Data/fish_abund-bout.rds")
sites.id <- unique(data.fish.bout$siteID)

### Reshuffle (1000x) species abundances, separately for each domain because domain will be used as a delineation of the regional species pool
# first order by species name (column)
cols <- sort(colnames(data.fish.bout[,8:ncol(data.fish.bout)]))
cols.bout <- c(colnames(data.fish.bout[,1:7]),cols)
data.fish.bout <- data.fish.bout[cols.bout]

# get all species ever recorded and then get list for each domain
data_small_fish <- readRDS(file = "Data/fish_processed.rds")
domain.sp <- data_small_fish %>%
  select(domainID, siteID, scientificName) %>%
  group_by(domainID, scientificName) %>%
  tally(name = "count") %>% 
  ungroup()

domains <- unique(domain.sp$domainID)

null.bout <- list()

system.time(
  for (k in 1:1000){
    cat(k)
    data.fish.bout <- rand.bout[[k]]
    sites.id <- unique(data.fish.bout$siteID)
    
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
    null.bout[[k]] <- out.bout.bat.all
  }
)

saveRDS(null.bout, file="data/fish_outBAT-FD-bout_null-consecutive-all.rds")



###### Quantify the standardized effect size (SES)
obs.bout.tax <- readRDS(file="Data/fish_outBAT-bout-consecutive-all.rds")
obs.bout <- readRDS(file="Data/fish_outBAT-FD-bout-consecutive-all.rds")
null.bout <- readRDS(file="Data/fish_outBAT-FD-bout_null-consecutive-all.rds")

# Jaccard dissimilarity
null.jacc.bout <- obs.bout %>%
  select(domainID, siteID, lat, lon, bout_lag, bout_consec_jacc)
for (i in 1:1000){
  null.jacc.i <- null.bout[[i]]
  null.jacc.i <- null.jacc.i %>%
    select(bout_consec_jacc)
  null.jacc.bout <- null.jacc.bout %>%
    bind_cols(null.jacc.i)
}

null.jacc.bout$bout_consec_jaccses <- NA
null.jacc.bout$bout_consec_jaccrank <- NA
null.jacc.bout$bout_consec_jaccpval <- NA
for (i in 1:nrow(null.jacc.bout)){
  null.jacc.bout[i,1007] <- ((as.numeric(null.jacc.bout[i,6]) - mean(as.numeric(null.jacc.bout[i,7:1006])))/sd(as.numeric(null.jacc.bout[i,7:1006])))
  rank.sim <- apply(null.jacc.bout[i,7:1006],1,rank)[,1]
  null.jacc.bout[i,1008] <- rank.sim[1]
  null.jacc.bout[i,1009] <- null.jacc.bout[i,1008]/1001
}

# Jaccard dissimilarity due to replacement 
null.jaccrepl.bout <- obs.bout %>%
  select(domainID, siteID, lat, lon, bout_lag, bout_consec_jaccrepl)
for (i in 1:1000){
  null.jaccrepl.i <- null.bout[[i]]
  null.jaccrepl.i <- null.jaccrepl.i %>%
    select(bout_consec_jaccrepl)
  null.jaccrepl.bout <- null.jaccrepl.bout %>%
    bind_cols(null.jaccrepl.i)
}

null.jaccrepl.bout$bout_consec_jaccreplses <- NA
null.jaccrepl.bout$bout_consec_jaccreplrank <- NA
null.jaccrepl.bout$bout_consec_jaccreplpval <- NA
for (i in 1:nrow(null.jaccrepl.bout)){
  null.jaccrepl.bout[i,1007] <- ((as.numeric(null.jaccrepl.bout[i,6]) - mean(as.numeric(null.jaccrepl.bout[i,7:1006])))/sd(as.numeric(null.jaccrepl.bout[i,7:1006])))
  rank.sim <- apply(null.jaccrepl.bout[i,7:1006],1,rank)[,1]
  null.jaccrepl.bout[i,1008] <- rank.sim[1]
  null.jaccrepl.bout[i,1009] <- null.jaccrepl.bout[i,1008]/1001
}

# Jaccard dissimilarity due to richness 
null.jaccrich.bout <- obs.bout %>%
  select(domainID, siteID, lat, lon, bout_lag, bout_consec_jaccrich)
for (i in 1:1000){
  null.jaccrich.i <- null.bout[[i]]
  null.jaccrich.i <- null.jaccrich.i %>%
    select(bout_consec_jaccrich)
  null.jaccrich.bout <- null.jaccrich.bout %>%
    bind_cols(null.jaccrich.i)
}

null.jaccrich.bout$bout_consec_jaccrichses <- NA
null.jaccrich.bout$bout_consec_jaccrichrank <- NA
null.jaccrich.bout$bout_consec_jaccrichpval <- NA
for (i in 1:nrow(null.jaccrich.bout)){
  null.jaccrich.bout[i,1007] <- ((as.numeric(null.jaccrich.bout[i,6]) - mean(as.numeric(null.jaccrich.bout[i,7:1006])))/sd(as.numeric(null.jaccrich.bout[i,7:1006])))
  rank.sim <- apply(null.jaccrich.bout[i,7:1006],1,rank)[,1]
  null.jaccrich.bout[i,1008] <- rank.sim[1]
  null.jaccrich.bout[i,1009] <- null.jaccrich.bout[i,1008]/1001
}

# Contribution of replacement to Jaccard dissimilarity
null.jaccrepl.contr.bout <- obs.bout %>%
  select(domainID, siteID, lat, lon, bout_lag, bout_consec_contr_jaccrepl)
for (i in 1:1000){
  null.jaccrepl.i <- null.bout[[i]]
  null.jaccrepl.i <- null.jaccrepl.i %>%
    select(bout_consec_contr_jaccrepl)
  null.jaccrepl.contr.bout <- null.jaccrepl.contr.bout %>%
    bind_cols(null.jaccrepl.i)
}

null.jaccrepl.contr.bout$bout_consec_contr_jaccreplses <- NA
null.jaccrepl.contr.bout$bout_consec_contr_jaccreplrank <- NA
null.jaccrepl.contr.bout$bout_consec_contr_jaccreplpval <- NA
for (i in 1:nrow(null.jaccrepl.bout)){
  null.jaccrepl.contr.bout[i,1007] <- ((as.numeric(null.jaccrepl.contr.bout[i,6]) - mean(as.numeric(null.jaccrepl.contr.bout[i,7:1006])))/sd(as.numeric(null.jaccrepl.contr.bout[i,7:1006])))
  rank.sim <- apply(null.jaccrepl.contr.bout[i,7:1006],1,rank)[,1]
  null.jaccrepl.contr.bout[i,1008] <- rank.sim[1]
  null.jaccrepl.contr.bout[i,1009] <- null.jaccrepl.contr.bout[i,1008]/1001
}


# Contribution of richness to Jaccard dissimilarity
null.jaccrich.contr.bout <- obs.bout %>%
  select(domainID, siteID, lat, lon, bout_lag, bout_consec_contr_jaccrich)
for (i in 1:1000){
  null.jaccrich.i <- null.bout[[i]]
  null.jaccrich.i <- null.jaccrich.i %>%
    select(bout_consec_contr_jaccrich)
  null.jaccrich.contr.bout <- null.jaccrich.contr.bout %>%
    bind_cols(null.jaccrich.i)
}

null.jaccrich.contr.bout$bout_consec_contr_jaccrichses <- NA
null.jaccrich.contr.bout$bout_consec_contr_jaccrichrank <- NA
null.jaccrich.contr.bout$bout_consec_contr_jaccrichpval <- NA
for (i in 1:nrow(null.jaccrich.bout)){
  null.jaccrich.contr.bout[i,1007] <- ((as.numeric(null.jaccrich.contr.bout[i,6]) - mean(as.numeric(null.jaccrich.contr.bout[i,7:1006])))/sd(as.numeric(null.jaccrich.contr.bout[i,7:1006])))
  rank.sim <- apply(null.jaccrich.contr.bout[i,7:1006],1,rank)[,1]
  null.jaccrich.contr.bout[i,1008] <- rank.sim[1]
  null.jaccrich.contr.bout[i,1009] <- null.jaccrich.contr.bout[i,1008]/1001
}

##subset and bind
null.jaccrepl.bout <- null.jaccrepl.bout %>%
  select(bout_consec_jaccrepl, bout_consec_jaccreplses, bout_consec_jaccreplrank, bout_consec_jaccreplpval) 
null.jaccrich.bout <- null.jaccrich.bout %>%
  select(bout_consec_jaccrich, bout_consec_jaccrichses, bout_consec_jaccrichrank, bout_consec_jaccrichpval) 
null.jaccrepl.contr.bout <- null.jaccrepl.contr.bout %>%
  select(bout_consec_contr_jaccrepl, bout_consec_contr_jaccreplses, bout_consec_contr_jaccreplrank, bout_consec_contr_jaccreplpval) 
null.jaccrich.contr.bout <- null.jaccrich.contr.bout %>%
  select(bout_consec_contr_jaccrich, bout_consec_contr_jaccrichses, bout_consec_contr_jaccrichrank, bout_consec_contr_jaccrichpval) 

mam.ses.bout <- null.jacc.bout %>%
  select(domainID, siteID, lat, lon, bout_lag, bout_consec_jacc, bout_consec_jaccses, bout_consec_jaccrank, bout_consec_jaccpval) %>%
  bind_cols(null.jaccrepl.bout) %>%
  bind_cols(null.jaccrich.bout) %>%
  bind_cols(null.jaccrepl.contr.bout) %>%
  bind_cols(null.jaccrich.contr.bout)
mam.ses.bout <- as_tibble(mam.ses.bout)

saveRDS(mam.ses.bout, file="data/fish_outBAT-FD-SES-bout_null-consecutive-all.rds")



########read in all data 
stab.bout.cons.all <- readRDS(file="data/fish_outStability-bout-consecutive-all.rds")
stab.bout.cons.all$bout_lag <- as.numeric(stab.bout.cons.all$bout_lag)
bat.bout.cons <- readRDS(file="data/fish_outBAT-FD-SES-bout_null-consecutive.rds")
bat.bout.cons$bout_lag <- as.numeric(bat.bout.cons$bout_lag)
bat.bout.cons.all <- readRDS(file="data/fish_outBAT-FD-SES-bout_null-consecutive-all.rds")
bat.bout.cons.all$bout_lag <- as.numeric(bat.bout.cons.all$bout_lag)

all.bout.cons.all <- cbind(stab.bout.cons.all, bat.bout.cons.all[,6:ncol(bat.bout.cons.all)])
saveRDS(all.bout.cons.all, file="data/fish_outAll-FD-SES-bout_null-consecutive-all.rds")




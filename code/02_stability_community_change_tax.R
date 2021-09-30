#Load packages
require(rhdf5)
require(here)
require(tidyverse)
require(BAT)
require(codyn)
require(tidyverse)
select <- dplyr::select


################### COMMUNITY STABILITY & (TAXONOMIC) COMMUNITY CHANGE USING FUNCTIONS FROM PACKAGES BAT AND CODYN

################### small mammals
data.mam.bout <- readRDS(file="data/mammal_abund-bout.rds")
sites.id <- unique(data.mam.bout$siteID)

###### Community stability
### Pair-wise community stability at bout-resolution
for (i in 1:length(sites.id)){
  mi <- data.mam.bout %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  yi <- nrow(mi)-1
  yi_ls <- list()
  for (p in 1:yi){
    ps <- seq((p+1),(yi+1),1)
    yi_ls[[p]] <- ps
  }
  yi_vec <- seq(2,(yi+1),1)
  yi_ls2 <- list()
  for (p in 1:yi){
    ps <- seq(1,p,1)
    yi_ls2[[(yi+1)-p]] <- ps
  }
  yi_vec2 <- seq(1,yi,1)
  

  for (k in 1:(nrow(mi)-1)){
    mik <- mi[k:nrow(mi),]
    yik <- nrow(mik)-1
    yik_vec <- seq(2,(yik+1),1)
    out.intra.stab <- matrix(NA, yik, 7)
    out.intra.stab <- as.data.frame(out.intra.stab)
  for (j in 1:length(yik_vec)){
    mikj <- mik[1:yik_vec[j],]
    
    cv <- sd(apply(mikj[,8:ncol(mikj)],1,sum))/mean(apply(mikj[,8:ncol(mikj)],1,sum))
    stab <- 1/cv
    out.intra.stab[j,1] <- c(domain)
    out.intra.stab[j,2] <- sites.id[i]
    out.intra.stab[j,3] <- lat
    out.intra.stab[j,4] <- lon
    out.intra.stab[j,5] <- yi_ls2[[k]][j]
    out.intra.stab[j,6] <- cv
    out.intra.stab[j,7] <- stab
  }
    if (k == 1) {
      out.intra.stab.all <- out.intra.stab
    } else {
      out.intra.stab.all <- out.intra.stab.all %>%
        bind_rows(out.intra.stab)
    }} 
  
  if (i == 1) {
    out.intra.stab.consecutive <- out.intra.stab.all
  } else {
    out.intra.stab.consecutive <- out.intra.stab.consecutive %>%
      bind_rows(out.intra.stab.all)
  }}

colnames(out.intra.stab.consecutive) <- c("domainID","siteID","lat", "lon", "bout_lag","bout_cv","bout_stability")
saveRDS(out.intra.stab.consecutive, file="Data/mammal_outStability-bout-consecutive-all.rds")
#out.intra.stab.consecutive <- readRDS(file="Data/mammal_outStability-bout-consecutive-all.rds")


###### Community change: Package BAT
##### pair-wise comparisons
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
  
  betadiv <- beta(mi[,8:ncol(mi)], abund = TRUE, raref = 0)
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
saveRDS(out.bout.bat.all, file="data/mammal_outBAT-bout-consecutive-all.rds")


###### Community change: Package BAT
##### pair-wise comparisons
mean.coord <- data.mam.bout %>%
  group_by(domainID, siteID) %>%
  summarize(lat = mean(lat), lon = mean(lon, na.rm=TRUE)) 

data.mam.bout.piv <- data.mam.bout %>%
  pivot_longer(-c(domainID, siteID, lat, lon, year, month, bout), names_to = "scientificName", values_to = "abundance") 

siteids <- unique(data.mam.bout.piv$siteID)

for (i in 1:length(siteids)){
  sel.site <- data.mam.bout.piv %>% 
    filter(siteID == siteids[i]) %>% 
    arrange(year,month)
  
  sel.time <- unique(sel.site$bout)
  sel.time <- as.data.frame(cbind(sel.time, seq(1,length(sel.time),1)))
  colnames(sel.time) <- c("bout","bout_no")
  sel.time$bout_no <- as.numeric(as.character(sel.time$bout_no))
  sel.site <- sel.site %>% left_join(sel.time, by="bout")
  yi <- nrow(sel.time)
  yi_ls <- list()
  for (p in 1:yi){
    ps <- seq(p,yi,1)
    yi_ls[[p]] <- ps
  }
  
  for (j in 1:(yi-1)){
    sel.site.j <- sel.site %>% filter(bout_no %in% yi_ls[[j]])
      ref.time <- as.numeric(yi_ls[[j]][1])
      out.intra.codyn <- RAC_change(sel.site.j, 
                                time.var = "bout_no", reference.time = ref.time, species.var = "scientificName",
                                abundance.var = "abundance", replicate.var = "siteID") %>%
    #remove bout comparisons that cross years
    mutate(bout_lag = bout_no2-bout_no) %>%
    left_join(mean.coord, by="siteID") %>%
    select(domainID, siteID, lat, lon, bout_no, bout_no2, bout_lag, richness_change, evenness_change, rank_change, gains, losses)
    
      if (j == 1) {
      out.intra.codyn.all <- out.intra.codyn
    } else {
      out.intra.codyn.all <- out.intra.codyn.all %>%
        bind_rows(out.intra.codyn)
    }}
    
    if (i == 1) {
      out.intra.codyn.cons.all <- out.intra.codyn.all
      } else {
        out.intra.codyn.cons.all <- out.intra.codyn.cons.all %>%
          bind_rows(out.intra.codyn.all)
  }}

saveRDS(out.intra.codyn.cons.all, file="Data/mammal_outCodyn-bout-consecutive-all.rds")

########read in all data
stab.bout.cons.all <- readRDS(file="data/mammal_outStability-bout-consecutive-all.rds")
stab.bout.cons.all$bout_lag <- as.numeric(stab.bout.cons.all$bout_lag)
bat.bout.cons.all <- readRDS(file="data/mammal_outBAT-bout-consecutive-all.rds")
bat.bout.cons.all$bout_lag <- as.numeric(bat.bout.cons.all$bout_lag)
codyn.bout.cons.all <- readRDS(file="data/mammal_outCodyn-bout-consecutive-all.rds")
codyn.bout.cons.all$bout_lag <- as.numeric(codyn.bout.cons.all$bout_lag)

all.bout.cons.all <- cbind(stab.bout.cons.all, bat.bout.cons.all[,6:ncol(bat.bout.cons.all)], codyn.bout.cons.all[,8:ncol(codyn.bout.cons.all)])
saveRDS(all.bout.cons.all, file="data/mammal_outAll-bout-consecutive-all.rds")


################### ground beetles
data.beetle.bout <- readRDS(file="data/beetle_abund-bout.rds")
sites.id <- unique(data.beetle.bout$siteID)

###### Community stability
### Pair-wise community stability at bout-resolution
for (i in 1:length(sites.id)){
  mi <- data.beetle.bout %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  yi <- nrow(mi)-1
  yi_ls <- list()
  for (p in 1:yi){
    ps <- seq((p+1),(yi+1),1)
    yi_ls[[p]] <- ps
  }
  yi_vec <- seq(2,(yi+1),1)
  yi_ls2 <- list()
  for (p in 1:yi){
    ps <- seq(1,p,1)
    yi_ls2[[(yi+1)-p]] <- ps
  }
  yi_vec2 <- seq(1,yi,1)
  
  
  for (k in 1:(nrow(mi)-1)){
    mik <- mi[k:nrow(mi),]
    yik <- nrow(mik)-1
    yik_vec <- seq(2,(yik+1),1)
    out.intra.stab <- matrix(NA, yik, 7)
    out.intra.stab <- as.data.frame(out.intra.stab)
    for (j in 1:length(yik_vec)){
      mikj <- mik[1:yik_vec[j],]
      
      cv <- sd(apply(mikj[,7:ncol(mikj)],1,sum))/mean(apply(mikj[,7:ncol(mikj)],1,sum))
      stab <- 1/cv
      out.intra.stab[j,1] <- c(domain)
      out.intra.stab[j,2] <- sites.id[i]
      out.intra.stab[j,3] <- lat
      out.intra.stab[j,4] <- lon
      out.intra.stab[j,5] <- yi_ls2[[k]][j]
      out.intra.stab[j,6] <- cv
      out.intra.stab[j,7] <- stab
    }
    if (k == 1) {
      out.intra.stab.all <- out.intra.stab
    } else {
      out.intra.stab.all <- out.intra.stab.all %>%
        bind_rows(out.intra.stab)
    }} 
  
  if (i == 1) {
    out.intra.stab.consecutive <- out.intra.stab.all
  } else {
    out.intra.stab.consecutive <- out.intra.stab.consecutive %>%
      bind_rows(out.intra.stab.all)
  }}

colnames(out.intra.stab.consecutive) <- c("domainID","siteID","lat", "lon", "bout_lag","bout_cv","bout_stability")
saveRDS(out.intra.stab.consecutive, file="data/beetle_outStability-bout-consecutive-all.rds")


###### Community change: Package BAT
### pair-wise comparisons
for (i in 1:length(sites.id)){
  mi <- data.beetle.bout %>%
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
  
  betadiv <- beta(mi[,7:ncol(mi)], abund = TRUE, raref = 0)
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
saveRDS(out.bout.bat.all, file="data/beetle_outBAT-bout-consecutive-all.rds")


###### Community change: Package codyn
### pair-wise comparisons
data.beetle.bout <- data.beetle.bout %>%
  mutate(year = as.numeric(substr(bout, start = 6, stop = 9)), month = as.numeric(substr(bout, start = 11, stop = 12)), day = as.numeric(substr(bout, start = 14, stop = 15)))

mean.coord <- data.beetle.bout %>%
  group_by(domainID, siteID) %>%
  summarize(lat = mean(lat), lon = mean(lon, na.rm=TRUE)) 

data.beetle.bout.piv <- data.beetle.bout %>%
  pivot_longer(-c(domainID, siteID, lat, lon, year, month, day, bout), names_to = "scientificName", values_to = "abundance") 

siteids <- unique(data.beetle.bout.piv$siteID)

for (i in 1:length(siteids)){
  sel.site <- data.beetle.bout.piv %>% 
    filter(siteID == siteids[i]) %>% 
    arrange(year,month,day)
  
  sel.time <- unique(sel.site$bout)
  sel.time <- as.data.frame(cbind(sel.time, seq(1,length(sel.time),1)))
  colnames(sel.time) <- c("bout","bout_no")
  sel.time$bout_no <- as.numeric(as.character(sel.time$bout_no))
  sel.site <- sel.site %>% left_join(sel.time, by="bout")
  yi <- nrow(sel.time)
  yi_ls <- list()
  for (p in 1:yi){
    ps <- seq(p,yi,1)
    yi_ls[[p]] <- ps
  }
  
  for (j in 1:(yi-1)){
    sel.site.j <- sel.site %>% filter(bout_no %in% yi_ls[[j]])
    ref.time <- as.numeric(yi_ls[[j]][1])
    out.intra.codyn <- RAC_change(sel.site.j, 
                                  time.var = "bout_no", reference.time = ref.time, species.var = "scientificName",
                                  abundance.var = "abundance", replicate.var = "siteID") %>%
      #remove bout comparisons that cross years
      mutate(bout_lag = bout_no2-bout_no) %>%
      left_join(mean.coord, by="siteID") %>%
      select(domainID, siteID, lat, lon, bout_no, bout_no2, bout_lag, richness_change, evenness_change, rank_change, gains, losses)
    
    if (j == 1) {
      out.intra.codyn.all <- out.intra.codyn
    } else {
      out.intra.codyn.all <- out.intra.codyn.all %>%
        bind_rows(out.intra.codyn)
    }}
  
  if (i == 1) {
    out.intra.codyn.cons.all <- out.intra.codyn.all
  } else {
    out.intra.codyn.cons.all <- out.intra.codyn.cons.all %>%
      bind_rows(out.intra.codyn.all)
  }}

saveRDS(out.intra.codyn.cons.all, file="data/beetle_outCodyn-bout-consecutive-all.rds")


########read in all data 
stab.bout.cons.all <- readRDS(file="data/beetle_outStability-bout-consecutive-all.rds")
stab.bout.cons.all$bout_lag <- as.numeric(stab.bout.cons.all$bout_lag)
bat.bout.cons.all <- readRDS(file="data/beetle_outBAT-bout-consecutive-all.rds")
bat.bout.cons.all$bout_lag <- as.numeric(bat.bout.cons.all$bout_lag)
codyn.bout.cons.all <- readRDS(file="data/beetle_outCodyn-bout-consecutive-all.rds")
codyn.bout.cons.all$bout_lag <- as.numeric(codyn.bout.cons.all$bout_lag)

all.bout.cons.all <- cbind(stab.bout.cons.all, bat.bout.cons.all[,6:ncol(bat.bout.cons.all)], codyn.bout.cons.all[,8:ncol(codyn.bout.cons.all)])
saveRDS(all.bout.cons.all, file="data/beetle_outAll-bout-consecutive-all.rds")



################### fish
data.fish.bout <- readRDS(file="data/fish_abund-bout.rds")
sites.id <- unique(data.fish.bout$siteID)

###### Community stability
### Pair-wise community stability at bout-resolution
for (i in 1:length(sites.id)){
  mi <- data.fish.bout %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  yi <- nrow(mi)-1
  yi_ls <- list()
  for (p in 1:yi){
    ps <- seq((p+1),(yi+1),1)
    yi_ls[[p]] <- ps
  }
  yi_vec <- seq(2,(yi+1),1)
  yi_ls2 <- list()
  for (p in 1:yi){
    ps <- seq(1,p,1)
    yi_ls2[[(yi+1)-p]] <- ps
  }
  yi_vec2 <- seq(1,yi,1)
  
  
  for (k in 1:(nrow(mi)-1)){
    mik <- mi[k:nrow(mi),]
    yik <- nrow(mik)-1
    yik_vec <- seq(2,(yik+1),1)
    out.intra.stab <- matrix(NA, yik, 7)
    out.intra.stab <- as.data.frame(out.intra.stab)
    for (j in 1:length(yik_vec)){
      mikj <- mik[1:yik_vec[j],]
      
      cv <- sd(apply(mikj[,8:ncol(mikj)],1,sum))/mean(apply(mikj[,8:ncol(mikj)],1,sum))
      stab <- 1/cv
      out.intra.stab[j,1] <- c(domain)
      out.intra.stab[j,2] <- sites.id[i]
      out.intra.stab[j,3] <- lat
      out.intra.stab[j,4] <- lon
      out.intra.stab[j,5] <- yi_ls2[[k]][j]
      out.intra.stab[j,6] <- cv
      out.intra.stab[j,7] <- stab
    }
    if (k == 1) {
      out.intra.stab.all <- out.intra.stab
    } else {
      out.intra.stab.all <- out.intra.stab.all %>%
        bind_rows(out.intra.stab)
    }} 
  
  if (i == 1) {
    out.intra.stab.consecutive <- out.intra.stab.all
  } else {
    out.intra.stab.consecutive <- out.intra.stab.consecutive %>%
      bind_rows(out.intra.stab.all)
  }}

colnames(out.intra.stab.consecutive) <- c("domainID","siteID","lat", "lon", "bout_lag","bout_cv","bout_stability")
saveRDS(out.intra.stab.consecutive, file="data/fish_outStability-bout-consecutive-all.rds")


###### Community change: Package BAT
### pair-wise comparisons
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
  
  betadiv <- beta(mi[,8:ncol(mi)], abund = TRUE, raref = 0)
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
saveRDS(out.bout.bat.all, file="data/fish_outBAT-bout-consecutive-all.rds")


###### Community change: Package codyn
### pair-wise comparisons
mean.coord <- data.fish.bout %>%
  group_by(domainID, siteID) %>%
  summarize(lat = mean(lat), lon = mean(lon, na.rm=TRUE)) 

data.fish.bout.piv <- data.fish.bout %>%
  pivot_longer(-c(domainID, siteID, lat, lon, year, month, bout), names_to = "scientificName", values_to = "abundance") 


siteids <- unique(data.fish.bout.piv$siteID)

for (i in 1:length(siteids)){
  sel.site <- data.fish.bout.piv %>% 
    filter(siteID == siteids[i]) %>% 
    arrange(year,month)
  
  sel.time <- unique(sel.site$bout)
  sel.time <- as.data.frame(cbind(sel.time, seq(1,length(sel.time),1)))
  colnames(sel.time) <- c("bout","bout_no")
  sel.time$bout_no <- as.numeric(as.character(sel.time$bout_no))
  sel.site <- sel.site %>% left_join(sel.time, by="bout")
  yi <- nrow(sel.time)
  yi_ls <- list()
  for (p in 1:yi){
    ps <- seq(p,yi,1)
    yi_ls[[p]] <- ps
  }
  
  for (j in 1:(yi-1)){
    sel.site.j <- sel.site %>% filter(bout_no %in% yi_ls[[j]])
    ref.time <- as.numeric(yi_ls[[j]][1])
    out.intra.codyn <- RAC_change(sel.site.j, 
                                  time.var = "bout_no", reference.time = ref.time, species.var = "scientificName",
                                  abundance.var = "abundance", replicate.var = "siteID") %>%
      #remove bout comparisons that cross years
      mutate(bout_lag = bout_no2-bout_no) %>%
      left_join(mean.coord, by="siteID") %>%
      select(domainID, siteID, lat, lon, bout_no, bout_no2, bout_lag, richness_change, evenness_change, rank_change, gains, losses)
    
    if (j == 1) {
      out.intra.codyn.all <- out.intra.codyn
    } else {
      out.intra.codyn.all <- out.intra.codyn.all %>%
        bind_rows(out.intra.codyn)
    }}
  
  if (i == 1) {
    out.intra.codyn.cons.all <- out.intra.codyn.all
  } else {
    out.intra.codyn.cons.all <- out.intra.codyn.cons.all %>%
      bind_rows(out.intra.codyn.all)
  }}

saveRDS(out.intra.codyn.cons.all, file="data/fish_outCodyn-bout-consecutive-all.rds")


########read in all data 
stab.bout.cons.all <- readRDS(file="data/fish_outStability-bout-consecutive-all.rds")
stab.bout.cons.all$bout_lag <- as.numeric(stab.bout.cons.all$bout_lag)
bat.bout.cons.all <- readRDS(file="data/fish_outBAT-bout-consecutive-all.rds")
bat.bout.cons.all$bout_lag <- as.numeric(bat.bout.cons.all$bout_lag)
codyn.bout.cons.all <- readRDS(file="data/fish_outCodyn-bout-consecutive-all.rds")
codyn.bout.cons.all$bout_lag <- as.numeric(codyn.bout.cons.all$bout_lag)

all.bout.cons.all <- cbind(stab.bout.cons.all, bat.bout.cons.all[,6:ncol(bat.bout.cons.all)], codyn.bout.cons.all[,8:ncol(codyn.bout.cons.all)])
saveRDS(all.bout.cons.all, file="data/fish_outAll-bout-consecutive-all.rds")



################### aquatic macroinvertebrates
data.macroinv.bout <- readRDS(file="data/macroinv_abund-bout.rds")
sites.id <- unique(data.macroinv.bout$siteID)

###### Community stability
### Pair-wise community stability at bout-resolution
for (i in 1:length(sites.id)){
  mi <- data.macroinv.bout %>%
    filter(siteID == sites.id[i])
  domain <- unique(mi$domainID)
  lat <- mean(mi$lat)
  lon <- mean(mi$lon)
  yi <- nrow(mi)-1
  yi_ls <- list()
  for (p in 1:yi){
    ps <- seq((p+1),(yi+1),1)
    yi_ls[[p]] <- ps
  }
  yi_vec <- seq(2,(yi+1),1)
  yi_ls2 <- list()
  for (p in 1:yi){
    ps <- seq(1,p,1)
    yi_ls2[[(yi+1)-p]] <- ps
  }
  yi_vec2 <- seq(1,yi,1)
  
  
  for (k in 1:(nrow(mi)-1)){
    mik <- mi[k:nrow(mi),]
    yik <- nrow(mik)-1
    yik_vec <- seq(2,(yik+1),1)
    out.intra.stab <- matrix(NA, yik, 7)
    out.intra.stab <- as.data.frame(out.intra.stab)
    for (j in 1:length(yik_vec)){
      mikj <- mik[1:yik_vec[j],]
      
      cv <- sd(apply(mikj[,8:ncol(mikj)],1,sum))/mean(apply(mikj[,8:ncol(mikj)],1,sum))
      stab <- 1/cv
      out.intra.stab[j,1] <- c(domain)
      out.intra.stab[j,2] <- sites.id[i]
      out.intra.stab[j,3] <- lat
      out.intra.stab[j,4] <- lon
      out.intra.stab[j,5] <- yi_ls2[[k]][j]
      out.intra.stab[j,6] <- cv
      out.intra.stab[j,7] <- stab
    }
    if (k == 1) {
      out.intra.stab.all <- out.intra.stab
    } else {
      out.intra.stab.all <- out.intra.stab.all %>%
        bind_rows(out.intra.stab)
    }} 
  
  if (i == 1) {
    out.intra.stab.consecutive <- out.intra.stab.all
  } else {
    out.intra.stab.consecutive <- out.intra.stab.consecutive %>%
      bind_rows(out.intra.stab.all)
  }}

colnames(out.intra.stab.consecutive) <- c("domainID","siteID","lat", "lon", "bout_lag","bout_cv","bout_stability")
saveRDS(out.intra.stab.consecutive, file="data/macroinv_outStability-bout-consecutive-all.rds")


###### Community change: Package BAT
### Pair-wise comparisons
for (i in 1:length(sites.id)){
  mi <- data.macroinv.bout %>%
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
  
  betadiv <- beta(mi[,8:ncol(mi)], abund = TRUE, raref = 0)
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
saveRDS(out.bout.bat.all, file="data/macroinv_outBAT-bout-consecutive-all.rds")

###### Community change: Package codyn
###pair-wise comparisons
mean.coord <- data.macroinv.bout %>%
  group_by(domainID, siteID) %>%
  summarize(lat = mean(lat), lon = mean(lon, na.rm=TRUE)) 

data.macroinv.bout.piv <- data.macroinv.bout %>%
  pivot_longer(-c(domainID, siteID, lat, lon, year, month, bout), names_to = "scientificName", values_to = "abundance") 

siteids <- unique(data.macroinv.bout.piv$siteID)

for (i in 1:length(siteids)){
  sel.site <- data.macroinv.bout.piv %>% 
    filter(siteID == siteids[i]) %>% 
    arrange(year,month)
  
  sel.time <- unique(sel.site$bout)
  sel.time <- as.data.frame(cbind(sel.time, seq(1,length(sel.time),1)))
  colnames(sel.time) <- c("bout","bout_no")
  sel.time$bout_no <- as.numeric(as.character(sel.time$bout_no))
  sel.site <- sel.site %>% left_join(sel.time, by="bout")
  yi <- nrow(sel.time)
  yi_ls <- list()
  for (p in 1:yi){
    ps <- seq(p,yi,1)
    yi_ls[[p]] <- ps
  }
  
  for (j in 1:(yi-1)){
    sel.site.j <- sel.site %>% filter(bout_no %in% yi_ls[[j]])
    ref.time <- as.numeric(yi_ls[[j]][1])
    out.intra.codyn <- RAC_change(sel.site.j, 
                                  time.var = "bout_no", reference.time = ref.time, species.var = "scientificName",
                                  abundance.var = "abundance", replicate.var = "siteID") %>%
      #remove bout comparisons that cross years
      mutate(bout_lag = bout_no2-bout_no) %>%
      left_join(mean.coord, by="siteID") %>%
      select(domainID, siteID, lat, lon, bout_no, bout_no2, bout_lag, richness_change, evenness_change, rank_change, gains, losses)
    
    if (j == 1) {
      out.intra.codyn.all <- out.intra.codyn
    } else {
      out.intra.codyn.all <- out.intra.codyn.all %>%
        bind_rows(out.intra.codyn)
    }}
  
  if (i == 1) {
    out.intra.codyn.cons.all <- out.intra.codyn.all
  } else {
    out.intra.codyn.cons.all <- out.intra.codyn.cons.all %>%
      bind_rows(out.intra.codyn.all)
  }}

saveRDS(out.intra.codyn.cons.all, file="data/macroinv_outCodyn-bout-consecutive-all.rds")

########read in all data 
stab.bout.cons.all <- readRDS(file="data/macroinv_outStability-bout-consecutive-all.rds")
stab.bout.cons.all$bout_lag <- as.numeric(stab.bout.cons.all$bout_lag)
bat.bout.cons.all <- readRDS(file="data/macroinv_outBAT-bout-consecutive-all.rds")
bat.bout.cons.all$bout_lag <- as.numeric(bat.bout.cons.all$bout_lag)
codyn.bout.cons.all <- readRDS(file="data/macroinv_outCodyn-bout-consecutive-all.rds")
codyn.bout.cons.all$bout_lag <- as.numeric(codyn.bout.cons.all$bout_lag)

all.bout.cons.all <- cbind(stab.bout.cons.all, bat.bout.cons.all[,6:ncol(bat.bout.cons.all)], codyn.bout.cons.all[,8:ncol(codyn.bout.cons.all)])
saveRDS(all.bout.cons.all, file="data/macroinv_outAll-bout-consecutive-all.rds")





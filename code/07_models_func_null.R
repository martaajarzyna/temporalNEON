#Load packages
require(rhdf5)
require(here)
require(tidyverse)
select <- dplyr::select
summarize <- dplyr::summarize
require(INLA)

###################MODELS EXPLAINING COMMUNITY STABILITY AS A FUNCTION OF CHANGE IN COMMUNITY FUNCTIONAL COMPOSITION

################### small mammals
### Read in data on stability and community compositional change
mammal.bout <- readRDS(file="data/mammal_outAll-FD-SES-bout_null-consecutive-all.rds")

mammal.bout <- mammal.bout %>% 
  filter(!(bout_lag == 1)) %>% # remove lag times = 1 
  filter_all(all_vars(!is.infinite(.))) %>% # remove ecosystem stability = Inf
  filter(bout_stability < 10) %>% # remove stability outliers
  mutate(bout_stability_sc = scale(bout_stability)) # standardize stability


##### MODELS
### 1. STABILITY vs FUNCTIONAL COMMUNITY CHANGE (SES)
## 1st set of models: Bout-level (x, x^2, x * time, x^2 * time; all with siteid as random effect)
siteids <- unique(mammal.bout$siteID)
siteids <- as.data.frame(cbind(siteids, seq(1,length(siteids),1)))
colnames(siteids) <- c("siteID","siteID1")
siteids$siteID1 <- as.numeric(as.character(siteids$siteID1))
mammal.bout <- mammal.bout %>%
  left_join(siteids,by="siteID")
mammal.bout$siteID2 <- mammal.bout$siteID1 
mammal.bout$siteID3 <- mammal.bout$siteID1

# Jaccard dissimilarity
form <- bout_stability_sc ~ bout_consec_jaccses + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccses, copy="siteID1")
m1 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-FD-SES-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_jaccses + I(bout_consec_jaccses^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccses, copy="siteID1") + f(siteID3, I(bout_consec_jaccses^2), copy="siteID1")
m2 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/mammal_stabsc-jacc-FD-SES-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_jaccses + bout_lag + bout_consec_jaccses:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccses, copy="siteID1")
m3 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/mammal_stabsc-jacc-FD-SES-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_jaccses + I(bout_consec_jaccses^2) + bout_lag + bout_consec_jaccses:bout_lag + I(bout_consec_jaccses^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccses, copy="siteID1") + f(siteID3, I(bout_consec_jaccses^2), copy="siteID1")
m4 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/mammal_stabsc-jacc-FD-SES-bout-m4.rds")


# Jaccard dissimilarity due to replacement
form <- bout_stability_sc ~ bout_consec_jaccreplses + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccreplses, copy="siteID1")
m1 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-FD-SES-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_jaccreplses + I(bout_consec_jaccreplses^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccreplses, copy="siteID1") + f(siteID3, I(bout_consec_jaccreplses^2), copy="siteID1")
m2 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/mammal_stabsc-jaccrepl-FD-SES-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_jaccreplses + bout_lag + bout_consec_jaccreplses:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccreplses, copy="siteID1")
m3 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/mammal_stabsc-jaccrepl-FD-SES-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_jaccreplses + I(bout_consec_jaccreplses^2) + bout_lag + bout_consec_jaccreplses:bout_lag + I(bout_consec_jaccreplses^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccreplses, copy="siteID1") + f(siteID3, I(bout_consec_jaccreplses^2), copy="siteID1")
m4 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/mammal_stabsc-jaccrepl-FD-SES-bout-m4.rds")


# Jaccard dissimilarity due to richness
form <- bout_stability_sc ~ bout_consec_jaccrichses + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrichses, copy="siteID1")
m1 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-FD-SES-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_jaccrichses + I(bout_consec_jaccrichses^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrichses, copy="siteID1") + f(siteID3, I(bout_consec_jaccrichses^2), copy="siteID1")
m2 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/mammal_stabsc-jaccrich-FD-SES-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_jaccrichses + bout_lag + bout_consec_jaccrichses:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrichses, copy="siteID1")
m3 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/mammal_stabsc-jaccrich-FD-SES-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_jaccrichses + I(bout_consec_jaccrichses^2) + bout_lag + bout_consec_jaccrichses:bout_lag + I(bout_consec_jaccrichses^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrichses, copy="siteID1") + f(siteID3, I(bout_consec_jaccrichses^2), copy="siteID1")
m4 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/mammal_stabsc-jaccrich-FD-SES-bout-m4.rds")


# Contribution of replacament to Jaccard dissimilarity
form <- bout_stability_sc ~ bout_consec_contr_jaccreplses + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccreplses, copy="siteID1")
m1 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccreplcontr-FD-SES-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccreplses + I(bout_consec_contr_jaccreplses^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccreplses, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccreplses^2), copy="siteID1")
m2 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/mammal_stabsc-jaccreplcontr-FD-SES-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccreplses + bout_lag + bout_consec_contr_jaccreplses:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccreplses, copy="siteID1")
m3 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/mammal_stabsc-jaccreplcontr-FD-SES-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccreplses + I(bout_consec_contr_jaccreplses^2) + bout_lag + bout_consec_contr_jaccreplses:bout_lag + I(bout_consec_contr_jaccreplses^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccreplses, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccreplses^2), copy="siteID1")
m4 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/mammal_stabsc-jaccreplcontr-FD-SES-bout-m4.rds")


# Contribution of richness to Jaccard dissimilarity
form <- bout_stability_sc ~ bout_consec_contr_jaccrichses + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrichses, copy="siteID1")
m1 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrichcontr-FD-SES-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrichses + I(bout_consec_contr_jaccrichses^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrichses, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccrichses^2), copy="siteID1")
m2 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/mammal_stabsc-jaccrichcontr-FD-SES-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrichses + bout_lag + bout_consec_contr_jaccrichses:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrichses, copy="siteID1")
m3 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/mammal_stabsc-jaccrichcontr-FD-SES-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrichses + I(bout_consec_contr_jaccrichses^2) + bout_lag + bout_consec_contr_jaccrichses:bout_lag + I(bout_consec_contr_jaccrichses^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrichses, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccrichses^2), copy="siteID1")
m4 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/mammal_stabsc-jaccrichcontr-FD-SES-bout-m4.rds")



################### fish
### Read in data on stability and community compositional change
fish.bout <- readRDS(file="data/fish_outAll-FD-SES-bout_null-consecutive-all.rds")

fish.bout <- fish.bout %>% 
  filter(!(bout_lag == 1)) %>% # remove lag times = 1 
  filter_all(all_vars(!is.infinite(.))) %>% # remove ecosystem stability = Inf
  filter(bout_stability < 10) %>% # remove stability outliers
  mutate(bout_stability_sc = scale(bout_stability)) # standardize stability


##### MODELS
### 1. STABILITY vs FUNCTIONAL COMMUNITY CHANGE (SES)
## 1st set of models: Bout-level (x, x^2, x * time, x^2 * time; all with siteid as random effect)
siteids <- unique(fish.bout$siteID)
siteids <- as.data.frame(cbind(siteids, seq(1,length(siteids),1)))
colnames(siteids) <- c("siteID","siteID1")
siteids$siteID1 <- as.numeric(as.character(siteids$siteID1))
fish.bout <- fish.bout %>%
  left_join(siteids,by="siteID")
fish.bout$siteID2 <- fish.bout$siteID1 
fish.bout$siteID3 <- fish.bout$siteID1

# Jaccard dissimilarity
form <- bout_stability_sc ~ bout_consec_jaccses + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccses, copy="siteID1")
m1 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-FD-SES-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_jaccses + I(bout_consec_jaccses^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccses, copy="siteID1") + f(siteID3, I(bout_consec_jaccses^2), copy="siteID1")
m2 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/fish_stabsc-jacc-FD-SES-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_jaccses + bout_lag + bout_consec_jaccses:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccses, copy="siteID1")
m3 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/fish_stabsc-jacc-FD-SES-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_jaccses + I(bout_consec_jaccses^2) + bout_lag + bout_consec_jaccses:bout_lag + I(bout_consec_jaccses^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccses, copy="siteID1") + f(siteID3, I(bout_consec_jaccses^2), copy="siteID1")
m4 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/fish_stabsc-jacc-FD-SES-bout-m4.rds")


# Jaccard dissimilarity due to replacement
form <- bout_stability_sc ~ bout_consec_jaccreplses + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccreplses, copy="siteID1")
m1 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-FD-SES-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_jaccreplses + I(bout_consec_jaccreplses^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccreplses, copy="siteID1") + f(siteID3, I(bout_consec_jaccreplses^2), copy="siteID1")
m2 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/fish_stabsc-jaccrepl-FD-SES-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_jaccreplses + bout_lag + bout_consec_jaccreplses:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccreplses, copy="siteID1")
m3 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/fish_stabsc-jaccrepl-FD-SES-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_jaccreplses + I(bout_consec_jaccreplses^2) + bout_lag + bout_consec_jaccreplses:bout_lag + I(bout_consec_jaccreplses^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccreplses, copy="siteID1") + f(siteID3, I(bout_consec_jaccreplses^2), copy="siteID1")
m4 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/fish_stabsc-jaccrepl-FD-SES-bout-m4.rds")


# Jaccard dissimilarity due to richness
form <- bout_stability_sc ~ bout_consec_jaccrichses + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrichses, copy="siteID1")
m1 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-FD-SES-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_jaccrichses + I(bout_consec_jaccrichses^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrichses, copy="siteID1") + f(siteID3, I(bout_consec_jaccrichses^2), copy="siteID1")
m2 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/fish_stabsc-jaccrich-FD-SES-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_jaccrichses + bout_lag + bout_consec_jaccrichses:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrichses, copy="siteID1")
m3 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/fish_stabsc-jaccrich-FD-SES-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_jaccrichses + I(bout_consec_jaccrichses^2) + bout_lag + bout_consec_jaccrichses:bout_lag + I(bout_consec_jaccrichses^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrichses, copy="siteID1") + f(siteID3, I(bout_consec_jaccrichses^2), copy="siteID1")
m4 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/fish_stabsc-jaccrich-FD-SES-bout-m4.rds")


# Contribution of replacament to Jaccard dissimilarity
form <- bout_stability_sc ~ bout_consec_contr_jaccreplses + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccreplses, copy="siteID1")
m1 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccreplcontr-FD-SES-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccreplses + I(bout_consec_contr_jaccreplses^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccreplses, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccreplses^2), copy="siteID1")
m2 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/fish_stabsc-jaccreplcontr-FD-SES-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccreplses + bout_lag + bout_consec_contr_jaccreplses:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccreplses, copy="siteID1")
m3 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/fish_stabsc-jaccreplcontr-FD-SES-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccreplses + I(bout_consec_contr_jaccreplses^2) + bout_lag + bout_consec_contr_jaccreplses:bout_lag + I(bout_consec_contr_jaccreplses^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccreplses, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccreplses^2), copy="siteID1")
m4 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/fish_stabsc-jaccreplcontr-FD-SES-bout-m4.rds")


# Contribution of richness to Jaccard dissimilarity
form <- bout_stability_sc ~ bout_consec_contr_jaccrichses + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrichses, copy="siteID1")
m1 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrichcontr-FD-SES-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrichses + I(bout_consec_contr_jaccrichses^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrichses, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccrichses^2), copy="siteID1")
m2 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/fish_stabsc-jaccrichcontr-FD-SES-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrichses + bout_lag + bout_consec_contr_jaccrichses:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrichses, copy="siteID1")
m3 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/fish_stabsc-jaccrichcontr-FD-SES-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrichses + I(bout_consec_contr_jaccrichses^2) + bout_lag + bout_consec_contr_jaccrichses:bout_lag + I(bout_consec_contr_jaccrichses^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrichses, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccrichses^2), copy="siteID1")
m4 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/fish_stabsc-jaccrichcontr-FD-SES-bout-m4.rds")


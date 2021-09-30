#Load packages
require(rhdf5)
require(here)
require(tidyverse)
select <- dplyr::select
summarize <- dplyr::summarize
require(INLA)

###################MODELS EXPLAINING COMMUNITY STABILITY AS A FUNCTION OF CHANGE IN COMMUNITY TAXONOMIC COMPOSITION

################### small mammals
### Read in data on stability and community compositional change
mammal.bout <- readRDS(file="data/mammal_outAll-bout-consecutive-all.rds")

mammal.bout <- mammal.bout %>% 
  filter(!(bout_lag == 1)) %>% # remove lag times = 1 
  filter_all(all_vars(!is.infinite(.))) %>% # remove community stability = Inf
  filter(bout_stability < 10) %>% # remove stability outliers
  mutate(bout_stability_sc = scale(bout_stability)) # standardize stability

##### MODELS
### 1. STABILITY vs TAXONOMIC COMMUNITY CHANGE
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
form <- bout_stability_sc ~ bout_consec_jacc + f(siteID1, model = "iid") + f(siteID2, bout_consec_jacc, copy="siteID1")
m1 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_jacc + I(bout_consec_jacc^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_jacc, copy="siteID1") + f(siteID3, I(bout_consec_jacc^2), copy="siteID1")
m2 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/mammal_stabsc-jacc-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_jacc + bout_lag + bout_consec_jacc:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jacc, copy="siteID1")
m3 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/mammal_stabsc-jacc-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_jacc + I(bout_consec_jacc^2) + bout_lag + bout_consec_jacc:bout_lag + I(bout_consec_jacc^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jacc, copy="siteID1") + f(siteID3, I(bout_consec_jacc^2), copy="siteID1")
m4 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/mammal_stabsc-jacc-bout-m4.rds")


# Jaccard dissimilarity due to replacement
form <- bout_stability_sc ~ bout_consec_jaccrepl + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrepl, copy="siteID1")
m1 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_jaccrepl + I(bout_consec_jaccrepl^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrepl, copy="siteID1") + f(siteID3, I(bout_consec_jaccrepl^2), copy="siteID1")
m2 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/mammal_stabsc-jaccrepl-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_jaccrepl + bout_lag + bout_consec_jaccrepl:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrepl, copy="siteID1")
m3 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/mammal_stabsc-jaccrepl-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_jaccrepl + I(bout_consec_jaccrepl^2) + bout_lag + bout_consec_jaccrepl:bout_lag + I(bout_consec_jaccrepl^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrepl, copy="siteID1") + f(siteID3, I(bout_consec_jaccrepl^2), copy="siteID1")
m4 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/mammal_stabsc-jaccrepl-bout-m4.rds")


# Jaccard dissimilarity due to richness
form <- bout_stability_sc ~ bout_consec_jaccrich + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrich, copy="siteID1")
m1 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_jaccrich + I(bout_consec_jaccrich^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrich, copy="siteID1") + f(siteID3, I(bout_consec_jaccrich^2), copy="siteID1")
m2 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/mammal_stabsc-jaccrich-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_jaccrich + bout_lag + bout_consec_jaccrich:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrich, copy="siteID1")
m3 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/mammal_stabsc-jaccrich-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_jaccrich + I(bout_consec_jaccrich^2) + bout_lag + bout_consec_jaccrich:bout_lag + I(bout_consec_jaccrich^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrich, copy="siteID1") + f(siteID3, I(bout_consec_jaccrich^2), copy="siteID1")
m4 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/mammal_stabsc-jaccrich-bout-m4.rds")


# Contribution of replacament to Jaccard dissimilarity
form <- bout_stability_sc ~ bout_consec_contr_jaccrepl + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrepl, copy="siteID1")
m1 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccreplcontr-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrepl + I(bout_consec_contr_jaccrepl^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrepl, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccrepl^2), copy="siteID1")
m2 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/mammal_stabsc-jaccreplcontr-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrepl + bout_lag + bout_consec_contr_jaccrepl:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrepl, copy="siteID1")
m3 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/mammal_stabsc-jaccreplcontr-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrepl + I(bout_consec_contr_jaccrepl^2) + bout_lag + bout_consec_contr_jaccrepl:bout_lag + I(bout_consec_contr_jaccrepl^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrepl, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccrepl^2), copy="siteID1")
m4 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/mammal_stabsc-jaccreplcontr-bout-m4.rds")


# Contribution of richness to Jaccard dissimilarity
form <- bout_stability_sc ~ bout_consec_contr_jaccrich + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrich, copy="siteID1")
m1 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrichcontr-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrich + I(bout_consec_contr_jaccrich^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrich, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccrich^2), copy="siteID1")
m2 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/mammal_stabsc-jaccrichcontr-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrich + bout_lag + bout_consec_contr_jaccrich:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrich, copy="siteID1")
m3 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/mammal_stabsc-jaccrichcontr-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrich + I(bout_consec_contr_jaccrich^2) + bout_lag + bout_consec_contr_jaccrich:bout_lag + I(bout_consec_contr_jaccrich^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrich, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccrich^2), copy="siteID1")
m4 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/mammal_stabsc-jaccrichcontr-bout-m4.rds")


# codyn: richness change
form <- bout_stability_sc ~ richness_change + f(siteID1, model = "iid") + f(siteID2, richness_change, copy="siteID1")
m1 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-richchange-bout-m1.rds")

form <- bout_stability_sc ~ richness_change + I(richness_change^2) + f(siteID1, model = "iid") + f(siteID2, richness_change, copy="siteID1") + f(siteID3, I(richness_change^2), copy="siteID1")
m2 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/mammal_stabsc-richchange-bout-m2.rds")

form <- bout_stability_sc ~ richness_change + bout_lag + richness_change:bout_lag + f(siteID1, model = "iid") + f(siteID2, richness_change, copy="siteID1")
m3 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/mammal_stabsc-richchange-bout-m3.rds")

form <- bout_stability_sc ~ richness_change + I(richness_change^2) + bout_lag + richness_change:bout_lag + I(richness_change^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, richness_change, copy="siteID1") + f(siteID3, I(richness_change^2), copy="siteID1")
m4 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/mammal_stabsc-richchange-bout-m4.rds")


# codyn: evenness change
form <- bout_stability_sc ~ evenness_change + f(siteID1, model = "iid") + f(siteID2, evenness_change, copy="siteID1")
m1 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-evenchange-bout-m1.rds")

form <- bout_stability_sc ~ evenness_change + I(evenness_change^2) + f(siteID1, model = "iid") + f(siteID2, evenness_change, copy="siteID1") + f(siteID3, I(evenness_change^2), copy="siteID1")
m2 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/mammal_stabsc-evenchange-bout-m2.rds")

form <- bout_stability_sc ~ evenness_change + bout_lag + evenness_change:bout_lag + f(siteID1, model = "iid") + f(siteID2, evenness_change, copy="siteID1")
m3 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/mammal_stabsc-evenchange-bout-m3.rds")

form <- bout_stability_sc ~ evenness_change + I(evenness_change^2) + bout_lag + evenness_change:bout_lag + I(evenness_change^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, evenness_change, copy="siteID1") + f(siteID3, I(evenness_change^2), copy="siteID1")
m4 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/mammal_stabsc-evenchange-bout-m4.rds")


# codyn: rank change
form <- bout_stability_sc ~ rank_change + f(siteID1, model = "iid") + f(siteID2, rank_change, copy="siteID1")
m1 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-rankchange-bout-m1.rds")

form <- bout_stability_sc ~ rank_change + I(rank_change^2) + f(siteID1, model = "iid") + f(siteID2, rank_change, copy="siteID1") + f(siteID3, I(rank_change^2), copy="siteID1")
m2 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/mammal_stabsc-rankchange-bout-m2.rds")

form <- bout_stability_sc ~ rank_change + bout_lag + rank_change:bout_lag + f(siteID1, model = "iid") + f(siteID2, rank_change, copy="siteID1")
m3 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/mammal_stabsc-rankchange-bout-m3.rds")

form <- bout_stability_sc ~ rank_change + I(rank_change^2) + bout_lag + rank_change:bout_lag + I(rank_change^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, rank_change, copy="siteID1") + f(siteID3, I(rank_change^2), copy="siteID1")
m4 <- inla(form, data=mammal.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/mammal_stabsc-rankchange-bout-m4.rds")



################### beetles
### Read in data on stability and community compositional change
beetle.bout <- readRDS(file="data/beetle_outAll-bout-consecutive-all.rds")

beetle.bout <- beetle.bout %>% 
  filter(!(bout_lag == 1)) %>% # remove lag times = 1 
  filter_all(all_vars(!is.infinite(.))) %>% # remove community stability = Inf
  filter(bout_stability < 10) %>% # remove stability outliers
  mutate(bout_stability_sc = scale(bout_stability)) # standardize stability


##### MODELS
### 1. STABILITY vs TAXONOMIC COMMUNITY CHANGE 
## 1st set of models: Bout-level (x, x^2, x * time, x^2 * time; all with siteid as random effect)
siteids <- unique(beetle.bout$siteID)
siteids <- as.data.frame(cbind(siteids, seq(1,length(siteids),1)))
colnames(siteids) <- c("siteID","siteID1")
siteids$siteID1 <- as.numeric(as.character(siteids$siteID1))
beetle.bout <- beetle.bout %>%
  left_join(siteids,by="siteID")
beetle.bout$siteID2 <- beetle.bout$siteID1 
beetle.bout$siteID3 <- beetle.bout$siteID1

# Jaccard dissimilarity
form <- bout_stability_sc ~ bout_consec_jacc + f(siteID1, model = "iid") + f(siteID2, bout_consec_jacc, copy="siteID1")
m1 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jacc-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_jacc + I(bout_consec_jacc^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_jacc, copy="siteID1") + f(siteID3, I(bout_consec_jacc^2), copy="siteID1")
m2 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/beetle_stabsc-jacc-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_jacc + bout_lag + bout_consec_jacc:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jacc, copy="siteID1")
m3 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/beetle_stabsc-jacc-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_jacc + I(bout_consec_jacc^2) + bout_lag + bout_consec_jacc:bout_lag + I(bout_consec_jacc^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jacc, copy="siteID1") + f(siteID3, I(bout_consec_jacc^2), copy="siteID1")
m4 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/beetle_stabsc-jacc-bout-m4.rds")


# Jaccard dissimilarity due to replacement
form <- bout_stability_sc ~ bout_consec_jaccrepl + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrepl, copy="siteID1")
m1 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jaccrepl-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_jaccrepl + I(bout_consec_jaccrepl^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrepl, copy="siteID1") + f(siteID3, I(bout_consec_jaccrepl^2), copy="siteID1")
m2 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/beetle_stabsc-jaccrepl-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_jaccrepl + bout_lag + bout_consec_jaccrepl:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrepl, copy="siteID1")
m3 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/beetle_stabsc-jaccrepl-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_jaccrepl + I(bout_consec_jaccrepl^2) + bout_lag + bout_consec_jaccrepl:bout_lag + I(bout_consec_jaccrepl^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrepl, copy="siteID1") + f(siteID3, I(bout_consec_jaccrepl^2), copy="siteID1")
m4 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/beetle_stabsc-jaccrepl-bout-m4.rds")


# Jaccard dissimilarity due to richness
form <- bout_stability_sc ~ bout_consec_jaccrich + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrich, copy="siteID1")
m1 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jaccrich-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_jaccrich + I(bout_consec_jaccrich^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrich, copy="siteID1") + f(siteID3, I(bout_consec_jaccrich^2), copy="siteID1")
m2 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/beetle_stabsc-jaccrich-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_jaccrich + bout_lag + bout_consec_jaccrich:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrich, copy="siteID1")
m3 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/beetle_stabsc-jaccrich-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_jaccrich + I(bout_consec_jaccrich^2) + bout_lag + bout_consec_jaccrich:bout_lag + I(bout_consec_jaccrich^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrich, copy="siteID1") + f(siteID3, I(bout_consec_jaccrich^2), copy="siteID1")
m4 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/beetle_stabsc-jaccrich-bout-m4.rds")


# Contribution of replacament to Jaccard dissimilarity
form <- bout_stability_sc ~ bout_consec_contr_jaccrepl + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrepl, copy="siteID1")
m1 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jaccreplcontr-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrepl + I(bout_consec_contr_jaccrepl^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrepl, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccrepl^2), copy="siteID1")
m2 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/beetle_stabsc-jaccreplcontr-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrepl + bout_lag + bout_consec_contr_jaccrepl:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrepl, copy="siteID1")
m3 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/beetle_stabsc-jaccreplcontr-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrepl + I(bout_consec_contr_jaccrepl^2) + bout_lag + bout_consec_contr_jaccrepl:bout_lag + I(bout_consec_contr_jaccrepl^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrepl, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccrepl^2), copy="siteID1")
m4 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/beetle_stabsc-jaccreplcontr-bout-m4.rds")


# Contribution of richness to Jaccard dissimilarity
form <- bout_stability_sc ~ bout_consec_contr_jaccrich + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrich, copy="siteID1")
m1 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jaccrichcontr-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrich + I(bout_consec_contr_jaccrich^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrich, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccrich^2), copy="siteID1")
m2 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/beetle_stabsc-jaccrichcontr-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrich + bout_lag + bout_consec_contr_jaccrich:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrich, copy="siteID1")
m3 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/beetle_stabsc-jaccrichcontr-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrich + I(bout_consec_contr_jaccrich^2) + bout_lag + bout_consec_contr_jaccrich:bout_lag + I(bout_consec_contr_jaccrich^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrich, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccrich^2), copy="siteID1")
m4 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/beetle_stabsc-jaccrichcontr-bout-m4.rds")


# codyn: richness change
form <- bout_stability_sc ~ richness_change + f(siteID1, model = "iid") + f(siteID2, richness_change, copy="siteID1")
m1 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-richchange-bout-m1.rds")

form <- bout_stability_sc ~ richness_change + I(richness_change^2) + f(siteID1, model = "iid") + f(siteID2, richness_change, copy="siteID1") + f(siteID3, I(richness_change^2), copy="siteID1")
m2 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/beetle_stabsc-richchange-bout-m2.rds")

form <- bout_stability_sc ~ richness_change + bout_lag + richness_change:bout_lag + f(siteID1, model = "iid") + f(siteID2, richness_change, copy="siteID1")
m3 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/beetle_stabsc-richchange-bout-m3.rds")

form <- bout_stability_sc ~ richness_change + I(richness_change^2) + bout_lag + richness_change:bout_lag + I(richness_change^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, richness_change, copy="siteID1") + f(siteID3, I(richness_change^2), copy="siteID1")
m4 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/beetle_stabsc-richchange-bout-m4.rds")


# codyn: evenness change
form <- bout_stability_sc ~ evenness_change + f(siteID1, model = "iid") + f(siteID2, evenness_change, copy="siteID1")
m1 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-evenchange-bout-m1.rds")

form <- bout_stability_sc ~ evenness_change + I(evenness_change^2) + f(siteID1, model = "iid") + f(siteID2, evenness_change, copy="siteID1") + f(siteID3, I(evenness_change^2), copy="siteID1")
m2 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/beetle_stabsc-evenchange-bout-m2.rds")

form <- bout_stability_sc ~ evenness_change + bout_lag + evenness_change:bout_lag + f(siteID1, model = "iid") + f(siteID2, evenness_change, copy="siteID1")
m3 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/beetle_stabsc-evenchange-bout-m3.rds")

form <- bout_stability_sc ~ evenness_change + I(evenness_change^2) + bout_lag + evenness_change:bout_lag + I(evenness_change^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, evenness_change, copy="siteID1") + f(siteID3, I(evenness_change^2), copy="siteID1")
m4 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/beetle_stabsc-evenchange-bout-m4.rds")


# codyn: rank change
form <- bout_stability_sc ~ rank_change + f(siteID1, model = "iid") + f(siteID2, rank_change, copy="siteID1")
m1 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-rankchange-bout-m1.rds")

form <- bout_stability_sc ~ rank_change + I(rank_change^2) + f(siteID1, model = "iid") + f(siteID2, rank_change, copy="siteID1") + f(siteID3, I(rank_change^2), copy="siteID1")
m2 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/beetle_stabsc-rankchange-bout-m2.rds")

form <- bout_stability_sc ~ rank_change + bout_lag + rank_change:bout_lag + f(siteID1, model = "iid") + f(siteID2, rank_change, copy="siteID1")
m3 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/beetle_stabsc-rankchange-bout-m3.rds")

form <- bout_stability_sc ~ rank_change + I(rank_change^2) + bout_lag + rank_change:bout_lag + I(rank_change^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, rank_change, copy="siteID1") + f(siteID3, I(rank_change^2), copy="siteID1")
m4 <- inla(form, data=beetle.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/beetle_stabsc-rankchange-bout-m4.rds")




################### fish
### Read in data on stability and community compositional change
fish.bout <- readRDS(file="data/fish_outAll-bout-consecutive-all.rds")

fish.bout <- fish.bout %>% 
  filter(!(bout_lag == 1)) %>% # remove lag times = 1 
  filter_all(all_vars(!is.infinite(.))) %>% # remove community stability = Inf
  filter(bout_stability < 10) %>% # remove stability outliers
  mutate(bout_stability_sc = scale(bout_stability)) # standardize stability


##### MODELS
### 1. STABILITY vs TAXONOMIC COMMUNITY CHANGE 
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
form <- bout_stability_sc ~ bout_consec_jacc + f(siteID1, model = "iid") + f(siteID2, bout_consec_jacc, copy="siteID1")
m1 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_jacc + I(bout_consec_jacc^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_jacc, copy="siteID1") + f(siteID3, I(bout_consec_jacc^2), copy="siteID1")
m2 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/fish_stabsc-jacc-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_jacc + bout_lag + bout_consec_jacc:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jacc, copy="siteID1")
m3 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/fish_stabsc-jacc-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_jacc + I(bout_consec_jacc^2) + bout_lag + bout_consec_jacc:bout_lag + I(bout_consec_jacc^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jacc, copy="siteID1") + f(siteID3, I(bout_consec_jacc^2), copy="siteID1")
m4 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/fish_stabsc-jacc-bout-m4.rds")


# Jaccard dissimilarity due to replacement
form <- bout_stability_sc ~ bout_consec_jaccrepl + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrepl, copy="siteID1")
m1 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_jaccrepl + I(bout_consec_jaccrepl^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrepl, copy="siteID1") + f(siteID3, I(bout_consec_jaccrepl^2), copy="siteID1")
m2 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/fish_stabsc-jaccrepl-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_jaccrepl + bout_lag + bout_consec_jaccrepl:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrepl, copy="siteID1")
m3 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/fish_stabsc-jaccrepl-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_jaccrepl + I(bout_consec_jaccrepl^2) + bout_lag + bout_consec_jaccrepl:bout_lag + I(bout_consec_jaccrepl^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrepl, copy="siteID1") + f(siteID3, I(bout_consec_jaccrepl^2), copy="siteID1")
m4 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/fish_stabsc-jaccrepl-bout-m4.rds")


# Jaccard dissimilarity due to richness
form <- bout_stability_sc ~ bout_consec_jaccrich + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrich, copy="siteID1")
m1 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_jaccrich + I(bout_consec_jaccrich^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrich, copy="siteID1") + f(siteID3, I(bout_consec_jaccrich^2), copy="siteID1")
m2 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/fish_stabsc-jaccrich-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_jaccrich + bout_lag + bout_consec_jaccrich:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrich, copy="siteID1")
m3 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/fish_stabsc-jaccrich-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_jaccrich + I(bout_consec_jaccrich^2) + bout_lag + bout_consec_jaccrich:bout_lag + I(bout_consec_jaccrich^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrich, copy="siteID1") + f(siteID3, I(bout_consec_jaccrich^2), copy="siteID1")
m4 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/fish_stabsc-jaccrich-bout-m4.rds")


# Contribution of replacament to Jaccard dissimilarity
form <- bout_stability_sc ~ bout_consec_contr_jaccrepl + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrepl, copy="siteID1")
m1 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccreplcontr-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrepl + I(bout_consec_contr_jaccrepl^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrepl, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccrepl^2), copy="siteID1")
m2 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/fish_stabsc-jaccreplcontr-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrepl + bout_lag + bout_consec_contr_jaccrepl:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrepl, copy="siteID1")
m3 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/fish_stabsc-jaccreplcontr-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrepl + I(bout_consec_contr_jaccrepl^2) + bout_lag + bout_consec_contr_jaccrepl:bout_lag + I(bout_consec_contr_jaccrepl^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrepl, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccrepl^2), copy="siteID1")
m4 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/fish_stabsc-jaccreplcontr-bout-m4.rds")


# Contribution of richness to Jaccard dissimilarity
form <- bout_stability_sc ~ bout_consec_contr_jaccrich + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrich, copy="siteID1")
m1 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrichcontr-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrich + I(bout_consec_contr_jaccrich^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrich, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccrich^2), copy="siteID1")
m2 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/fish_stabsc-jaccrichcontr-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrich + bout_lag + bout_consec_contr_jaccrich:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrich, copy="siteID1")
m3 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/fish_stabsc-jaccrichcontr-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrich + I(bout_consec_contr_jaccrich^2) + bout_lag + bout_consec_contr_jaccrich:bout_lag + I(bout_consec_contr_jaccrich^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrich, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccrich^2), copy="siteID1")
m4 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/fish_stabsc-jaccrichcontr-bout-m4.rds")

# codyn: richness change
form <- bout_stability_sc ~ richness_change + f(siteID1, model = "iid") + f(siteID2, richness_change, copy="siteID1")
m1 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-richchange-bout-m1.rds")

form <- bout_stability_sc ~ richness_change + I(richness_change^2) + f(siteID1, model = "iid") + f(siteID2, richness_change, copy="siteID1") + f(siteID3, I(richness_change^2), copy="siteID1")
m2 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/fish_stabsc-richchange-bout-m2.rds")

form <- bout_stability_sc ~ richness_change + bout_lag + richness_change:bout_lag + f(siteID1, model = "iid") + f(siteID2, richness_change, copy="siteID1")
m3 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/fish_stabsc-richchange-bout-m3.rds")

form <- bout_stability_sc ~ richness_change + I(richness_change^2) + bout_lag + richness_change:bout_lag + I(richness_change^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, richness_change, copy="siteID1") + f(siteID3, I(richness_change^2), copy="siteID1")
m4 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/fish_stabsc-richchange-bout-m4.rds")


# codyn: evenness change
form <- bout_stability_sc ~ evenness_change + f(siteID1, model = "iid") + f(siteID2, evenness_change, copy="siteID1")
m1 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-evenchange-bout-m1.rds")

form <- bout_stability_sc ~ evenness_change + I(evenness_change^2) + f(siteID1, model = "iid") + f(siteID2, evenness_change, copy="siteID1") + f(siteID3, I(evenness_change^2), copy="siteID1")
m2 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/fish_stabsc-evenchange-bout-m2.rds")

form <- bout_stability_sc ~ evenness_change + bout_lag + evenness_change:bout_lag + f(siteID1, model = "iid") + f(siteID2, evenness_change, copy="siteID1")
m3 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/fish_stabsc-evenchange-bout-m3.rds")

form <- bout_stability_sc ~ evenness_change + I(evenness_change^2) + bout_lag + evenness_change:bout_lag + I(evenness_change^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, evenness_change, copy="siteID1") + f(siteID3, I(evenness_change^2), copy="siteID1")
m4 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/fish_stabsc-evenchange-bout-m4.rds")


# codyn: rank change
form <- bout_stability_sc ~ rank_change + f(siteID1, model = "iid") + f(siteID2, rank_change, copy="siteID1")
m1 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-rankchange-bout-m1.rds")

form <- bout_stability_sc ~ rank_change + I(rank_change^2) + f(siteID1, model = "iid") + f(siteID2, rank_change, copy="siteID1") + f(siteID3, I(rank_change^2), copy="siteID1")
m2 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/fish_stabsc-rankchange-bout-m2.rds")

form <- bout_stability_sc ~ rank_change + bout_lag + rank_change:bout_lag + f(siteID1, model = "iid") + f(siteID2, rank_change, copy="siteID1")
m3 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/fish_stabsc-rankchange-bout-m3.rds")

form <- bout_stability_sc ~ rank_change + I(rank_change^2) + bout_lag + rank_change:bout_lag + I(rank_change^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, rank_change, copy="siteID1") + f(siteID3, I(rank_change^2), copy="siteID1")
m4 <- inla(form, data=fish.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/fish_stabsc-rankchange-bout-m4.rds")


################### macroinv
### Read in data on stability and community compositional change
macroinv.bout <- readRDS(file="Data/macroinv_outAll-bout-consecutive-all.rds")


macroinv.bout <- macroinv.bout %>% 
  filter(!(bout_lag == 1)) %>% # remove lag times = 1 
  filter_all(all_vars(!is.infinite(.))) %>% # remove community stability = Inf
  filter(bout_stability < 10) %>% # remove stability outliers
  mutate(bout_stability_sc = scale(bout_stability)) # standardize stability

##### MODELS
### 1. STABILITY vs TAXONOMIC COMMUNITY CHANGE 
## 1st set of models: Bout-level (x, x^2, x * time, x^2 * time; all with siteid as random effect)
siteids <- unique(macroinv.bout$siteID)
siteids <- as.data.frame(cbind(siteids, seq(1,length(siteids),1)))
colnames(siteids) <- c("siteID","siteID1")
siteids$siteID1 <- as.numeric(as.character(siteids$siteID1))
macroinv.bout <- macroinv.bout %>%
  left_join(siteids,by="siteID")
macroinv.bout$siteID2 <- macroinv.bout$siteID1 
macroinv.bout$siteID3 <- macroinv.bout$siteID1

# Jaccard dissimilarity
form <- bout_stability_sc ~ bout_consec_jacc + f(siteID1, model = "iid") + f(siteID2, bout_consec_jacc, copy="siteID1")
m1 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jacc-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_jacc + I(bout_consec_jacc^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_jacc, copy="siteID1") + f(siteID3, I(bout_consec_jacc^2), copy="siteID1")
m2 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/macroinv_stabsc-jacc-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_jacc + bout_lag + bout_consec_jacc:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jacc, copy="siteID1")
m3 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/macroinv_stabsc-jacc-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_jacc + I(bout_consec_jacc^2) + bout_lag + bout_consec_jacc:bout_lag + I(bout_consec_jacc^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jacc, copy="siteID1") + f(siteID3, I(bout_consec_jacc^2), copy="siteID1")
m4 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/macroinv_stabsc-jacc-bout-m4.rds")


# Jaccard dissimilarity due to replacement
form <- bout_stability_sc ~ bout_consec_jaccrepl + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrepl, copy="siteID1")
m1 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jaccrepl-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_jaccrepl + I(bout_consec_jaccrepl^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrepl, copy="siteID1") + f(siteID3, I(bout_consec_jaccrepl^2), copy="siteID1")
m2 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/macroinv_stabsc-jaccrepl-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_jaccrepl + bout_lag + bout_consec_jaccrepl:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrepl, copy="siteID1")
m3 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/macroinv_stabsc-jaccrepl-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_jaccrepl + I(bout_consec_jaccrepl^2) + bout_lag + bout_consec_jaccrepl:bout_lag + I(bout_consec_jaccrepl^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrepl, copy="siteID1") + f(siteID3, I(bout_consec_jaccrepl^2), copy="siteID1")
m4 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/macroinv_stabsc-jaccrepl-bout-m4.rds")


# Jaccard dissimilarity due to richness
form <- bout_stability_sc ~ bout_consec_jaccrich + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrich, copy="siteID1")
m1 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jaccrich-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_jaccrich + I(bout_consec_jaccrich^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrich, copy="siteID1") + f(siteID3, I(bout_consec_jaccrich^2), copy="siteID1")
m2 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/macroinv_stabsc-jaccrich-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_jaccrich + bout_lag + bout_consec_jaccrich:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrich, copy="siteID1")
m3 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/macroinv_stabsc-jaccrich-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_jaccrich + I(bout_consec_jaccrich^2) + bout_lag + bout_consec_jaccrich:bout_lag + I(bout_consec_jaccrich^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_jaccrich, copy="siteID1") + f(siteID3, I(bout_consec_jaccrich^2), copy="siteID1")
m4 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/macroinv_stabsc-jaccrich-bout-m4.rds")


# Contribution of replacament to Jaccard dissimilarity
form <- bout_stability_sc ~ bout_consec_contr_jaccrepl + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrepl, copy="siteID1")
m1 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jaccreplcontr-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrepl + I(bout_consec_contr_jaccrepl^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrepl, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccrepl^2), copy="siteID1")
m2 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/macroinv_stabsc-jaccreplcontr-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrepl + bout_lag + bout_consec_contr_jaccrepl:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrepl, copy="siteID1")
m3 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/macroinv_stabsc-jaccreplcontr-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrepl + I(bout_consec_contr_jaccrepl^2) + bout_lag + bout_consec_contr_jaccrepl:bout_lag + I(bout_consec_contr_jaccrepl^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrepl, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccrepl^2), copy="siteID1")
m4 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/macroinv_stabsc-jaccreplcontr-bout-m4.rds")


# Contribution of richness to Jaccard dissimilarity
form <- bout_stability_sc ~ bout_consec_contr_jaccrich + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrich, copy="siteID1")
m1 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jaccrichcontr-bout-m1.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrich + I(bout_consec_contr_jaccrich^2) + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrich, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccrich^2), copy="siteID1")
m2 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/macroinv_stabsc-jaccrichcontr-bout-m2.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrich + bout_lag + bout_consec_contr_jaccrich:bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrich, copy="siteID1")
m3 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/macroinv_stabsc-jaccrichcontr-bout-m3.rds")

form <- bout_stability_sc ~ bout_consec_contr_jaccrich + I(bout_consec_contr_jaccrich^2) + bout_lag + bout_consec_contr_jaccrich:bout_lag + I(bout_consec_contr_jaccrich^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, bout_consec_contr_jaccrich, copy="siteID1") + f(siteID3, I(bout_consec_contr_jaccrich^2), copy="siteID1")
m4 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/macroinv_stabsc-jaccrichcontr-bout-m4.rds")


# codyn: richness change
form <- bout_stability_sc ~ richness_change + f(siteID1, model = "iid") + f(siteID2, richness_change, copy="siteID1")
m1 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-richchange-bout-m1.rds")

form <- bout_stability_sc ~ richness_change + I(richness_change^2) + f(siteID1, model = "iid") + f(siteID2, richness_change, copy="siteID1") + f(siteID3, I(richness_change^2), copy="siteID1")
m2 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/macroinv_stabsc-richchange-bout-m2.rds")

form <- bout_stability_sc ~ richness_change + bout_lag + richness_change:bout_lag + f(siteID1, model = "iid") + f(siteID2, richness_change, copy="siteID1")
m3 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/macroinv_stabsc-richchange-bout-m3.rds")

form <- bout_stability_sc ~ richness_change + I(richness_change^2) + bout_lag + richness_change:bout_lag + I(richness_change^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, richness_change, copy="siteID1") + f(siteID3, I(richness_change^2), copy="siteID1")
m4 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/macroinv_stabsc-richchange-bout-m4.rds")


# codyn: evenness change
form <- bout_stability_sc ~ evenness_change + f(siteID1, model = "iid") + f(siteID2, evenness_change, copy="siteID1")
m1 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-evenchange-bout-m1.rds")

form <- bout_stability_sc ~ evenness_change + I(evenness_change^2) + f(siteID1, model = "iid") + f(siteID2, evenness_change, copy="siteID1") + f(siteID3, I(evenness_change^2), copy="siteID1")
m2 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/macroinv_stabsc-evenchange-bout-m2.rds")

form <- bout_stability_sc ~ evenness_change + bout_lag + evenness_change:bout_lag + f(siteID1, model = "iid") + f(siteID2, evenness_change, copy="siteID1")
m3 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/macroinv_stabsc-evenchange-bout-m3.rds")

form <- bout_stability_sc ~ evenness_change + I(evenness_change^2) + bout_lag + evenness_change:bout_lag + I(evenness_change^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, evenness_change, copy="siteID1") + f(siteID3, I(evenness_change^2), copy="siteID1")
m4 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/macroinv_stabsc-evenchange-bout-m4.rds")


# codyn: rank change
form <- bout_stability_sc ~ rank_change + f(siteID1, model = "iid") + f(siteID2, rank_change, copy="siteID1")
m1 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-rankchange-bout-m1.rds")

form <- bout_stability_sc ~ rank_change + I(rank_change^2) + f(siteID1, model = "iid") + f(siteID2, rank_change, copy="siteID1") + f(siteID3, I(rank_change^2), copy="siteID1")
m2 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m2, file="Models_out/macroinv_stabsc-rankchange-bout-m2.rds")

form <- bout_stability_sc ~ rank_change + bout_lag + rank_change:bout_lag + f(siteID1, model = "iid") + f(siteID2, rank_change, copy="siteID1")
m3 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m3, file="Models_out/macroinv_stabsc-rankchange-bout-m3.rds")

form <- bout_stability_sc ~ rank_change + I(rank_change^2) + bout_lag + rank_change:bout_lag + I(rank_change^2):bout_lag + f(siteID1, model = "iid") + f(siteID2, rank_change, copy="siteID1") + f(siteID3, I(rank_change^2), copy="siteID1")
m4 <- inla(form, data=macroinv.bout, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m4, file="Models_out/macroinv_stabsc-rankchange-bout-m4.rds")



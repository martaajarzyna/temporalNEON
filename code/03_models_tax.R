#Load packages
require(rhdf5)
require(here)
require(tidyverse)
select <- dplyr::select
summarize <- dplyr::summarize
require(INLA)

###################MODELS EXPLAINING ECOSYSTEM STABILITY AS A FUNCTION OF CHANGE IN COMMUNITY TAXONOMIC COMPOSITION

################### small mammals
### Read in data on stability and community compositional change
mammal.bout.stab <- readRDS(file="data/mammal_outStability-bout.rds") %>%
  select(siteID, year, bout_cv, bout_stability)
mammal.year.stab <- readRDS(file="data/mammal_outStability-year.rds") %>%
  select(siteID, year_cv, year_stability)
mammal.bout.bat <- readRDS(file="data/mammal_outBAT-bout.rds")
mammal.year.bat <- readRDS(file="data/mammal_outBAT-year.rds")
mammal.bout.codyn <- readRDS(file="data/mammal_outCodyn-bout.rds")
mammal.year.codyn <- readRDS(file="data/mammal_outCodyn-year.rds")

# quantify contributions of replacement and richness to dissimilarity 
mammal.bout.bat <- mammal.bout.bat %>% 
  mutate(bout_contr_jaccrepl = bout_jaccrepl/bout_jacc, bout_contr_jaccrich = bout_jaccrich/bout_jacc)
mammal.year.bat <- mammal.year.bat %>% 
  mutate(year_contr_jaccrepl = year_jaccrepl/year_jacc, year_contr_jaccrich = year_jaccrich/year_jacc)

# join with stability
mammal.bout.bat <- mammal.bout.bat %>%
  left_join(mammal.bout.stab, by = c("siteID","year"))
mammal.year.bat <- mammal.year.bat %>%
  left_join(mammal.year.stab, by = c("siteID"))
mammal.bout.codyn <- mammal.bout.codyn %>%
  left_join(mammal.bout.stab, by = c("siteID","year"))
mammal.year.codyn <- mammal.year.codyn %>%
  left_join(mammal.year.stab, by = c("siteID"))

# remove stability outliers
mammal.bout.bat <- mammal.bout.bat %>%
  filter(bout_stability < 10)
mammal.bout.codyn <- mammal.bout.codyn %>%
  filter(bout_stability < 10)
mammal.year.bat <- mammal.year.bat %>%
  filter(year_stability < 5)
mammal.year.codyn <- mammal.year.codyn %>%
  filter(year_stability < 5)

# obtain site-level means
mammal.bout.bat.m <- mammal.bout.bat %>%
  dplyr::group_by(domainID, siteID, lat, lon, year, bout_cv, bout_stability) %>%
  dplyr::summarize(bout_jacc = mean(bout_jacc, na.rm=TRUE), bout_jaccrepl = mean(bout_jaccrepl, na.rm=TRUE), bout_jaccrich = mean(bout_jaccrich, na.rm=TRUE),
            bout_contr_jaccrepl = mean(bout_contr_jaccrepl, na.rm=TRUE), bout_contr_jaccrich = mean(bout_contr_jaccrich, na.rm=TRUE)) %>%
  dplyr::select(domainID, siteID, lat, lon, year, bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich, bout_cv, bout_stability) %>%
  ungroup

mammal.year.bat.m <- mammal.year.bat %>%
  group_by(domainID, siteID, lat, lon, year_cv, year_stability) %>%
  summarize(year_jacc = mean(year_jacc, na.rm=TRUE), year_jaccrepl = mean(year_jaccrepl, na.rm=TRUE), year_jaccrich = mean(year_jaccrich, na.rm=TRUE),
            year_contr_jaccrepl = mean(year_contr_jaccrepl, na.rm=TRUE), year_contr_jaccrich = mean(year_contr_jaccrich, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich, year_cv, year_stability) %>%
  ungroup

mammal.bout.codyn.m <- mammal.bout.codyn %>%
  group_by(domainID, siteID, lat, lon, year, bout_cv, bout_stability) %>%
  summarize(richness_change = mean(richness_change, na.rm=TRUE), evenness_change = mean(evenness_change, na.rm=TRUE), rank_change = mean(rank_change, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, year, richness_change, evenness_change, rank_change, bout_cv, bout_stability) %>%
  ungroup

mammal.year.codyn.m <- mammal.year.codyn %>%
  group_by(domainID, siteID, lat, lon, year_cv, year_stability) %>%
  summarize(richness_change = mean(richness_change, na.rm=TRUE), evenness_change = mean(evenness_change, na.rm=TRUE), rank_change = mean(rank_change, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, richness_change, evenness_change, rank_change, year_cv, year_stability) %>%
  ungroup


# scale stability, latitude, and longitude, remove NAs, replace 0 and 1 values for BAT metrics and rank change with 0.0001 and 0.999 so that models with beta distribution can be fit
mammal.bout.bat <- mammal.bout.bat %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(bout_stability)) %>%
  mutate_at(vars(bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich), ~replace(.,. == 1, 0.9999))

mammal.year.bat <- mammal.year.bat %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(year_stability)) %>%
  mutate_at(vars(year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich), ~replace(.,. == 1, 0.9999)) 

mammal.bout.codyn <- mammal.bout.codyn %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(bout_stability)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 1, 0.9999))

mammal.year.codyn <- mammal.year.codyn %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(year_stability)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 1, 0.9999))

mammal.bout.bat.m <- mammal.bout.bat.m %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(bout_stability)) %>%
  mutate_at(vars(bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich), ~replace(.,. == 1, 0.9999))

mammal.year.bat.m <- mammal.year.bat.m %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(year_stability)) %>%
  mutate_at(vars(year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich), ~replace(.,. == 1, 0.9999)) 

mammal.bout.codyn.m <- mammal.bout.codyn.m %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(bout_stability)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 1, 0.9999))

mammal.year.codyn.m <- mammal.year.codyn.m %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(year_stability)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 1, 0.9999))



##### MODELS
### 1. STABILITY vs TAXONOMIC COMMUNITY CHANGE (total of 32 models)
## 1st set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, BAT
# Jaccard dissimilarity
form <- stabsc ~ bout_jacc
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-bout-lin.rds")

form <- stabsc ~ bout_jacc + I(bout_jacc^2)
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-bout-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ bout_jaccrepl
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-bout-lin.rds")

form <- stabsc ~ bout_jaccrepl + I(bout_jaccrepl^2)
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-bout-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ bout_jaccrich
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-bout-lin.rds")

form <- stabsc ~ bout_jaccrich + I(bout_jaccrich^2)
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-bout-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrepl
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-bout-lin.rds")

form <- stabsc ~ bout_contr_jaccrepl + I(bout_contr_jaccrepl^2)
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-bout-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrich
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-bout-lin.rds")

form <- stabsc ~ bout_contr_jaccrich + I(bout_contr_jaccrich^2)
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-bout-quad.rds")


## 2nd set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, Codyn
# Richness change
form <- stabsc ~ richness_change
m1 <- inla(form, data=mammal.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-rich-bout-lin.rds")

form <- stabsc ~ richness_change + I(richness_change^2)
m1 <- inla(form, data=mammal.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-rich-bout-quad.rds")

# Evenness change
form <- stabsc ~ evenness_change
m1 <- inla(form, data=mammal.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-even-bout-lin.rds")

form <- stabsc ~ evenness_change + I(evenness_change^2)
m1 <- inla(form, data=mammal.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-even-bout-quad.rds")

# Rank change
form <- stabsc ~ rank_change
m1 <- inla(form, data=mammal.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-rank-bout-lin.rds")

form <- stabsc ~ rank_change + I(rank_change^2)
m1 <- inla(form, data=mammal.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-rank-bout-quad.rds")


## 3rd set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, BAT: site-level means
# Jaccard dissimilarity
form <- stabsc ~ bout_jacc
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-bout-lin-S.rds")

form <- stabsc ~ bout_jacc + I(bout_jacc^2)
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-bout-quad-S.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ bout_jaccrepl
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-bout-lin-S.rds")

form <- stabsc ~ bout_jaccrepl + I(bout_jaccrepl^2)
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-bout-quad-S.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ bout_jaccrich
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-bout-lin-S.rds")

form <- stabsc ~ bout_jaccrich + I(bout_jaccrich^2)
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-bout-quad-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrepl
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-bout-lin-S.rds")

form <- stabsc ~ bout_contr_jaccrepl + I(bout_contr_jaccrepl^2)
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-bout-quad-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrich
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-bout-lin-S.rds")

form <- stabsc ~ bout_contr_jaccrich + I(bout_contr_jaccrich^2)
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-bout-quad-S.rds")


## 4th set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, Codyn: site-level means
# Richness change
form <- stabsc ~ richness_change
m1 <- inla(form, data=mammal.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-rich-bout-lin-S.rds")

form <- stabsc ~ richness_change + I(richness_change^2)
m1 <- inla(form, data=mammal.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-rich-bout-quad-S.rds")

# Evenness change
form <- stabsc ~ evenness_change
m1 <- inla(form, data=mammal.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-even-bout-lin-S.rds")

form <- stabsc ~ evenness_change + I(evenness_change^2)
m1 <- inla(form, data=mammal.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-even-bout-quad-S.rds")

# Rank change
form <- stabsc ~ rank_change
m1 <- inla(form, data=mammal.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-rank-bout-lin-S.rds")

form <- stabsc ~ rank_change + I(rank_change^2)
m1 <- inla(form, data=mammal.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-rank-bout-quad-S.rds")


## 5th set of models: Year-level (inter-annual) stability ~ Year-level (inter-annual) community change, BAT
# Jaccard dissimilarity
form <- stabsc ~ year_jacc
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-year-lin.rds")

form <- stabsc ~ year_jacc + I(year_jacc^2)
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-year-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ year_jaccrepl
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-year-lin.rds")

form <- stabsc ~ year_jaccrepl + I(year_jaccrepl^2)
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-year-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ year_jaccrich
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-year-lin.rds")

form <- stabsc ~ year_jaccrich + I(year_jaccrich^2)
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-year-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrepl
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-year-lin.rds")

form <- stabsc ~ year_contr_jaccrepl + I(year_contr_jaccrepl^2)
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-year-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrich
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-year-lin.rds")

form <- stabsc ~ year_contr_jaccrich + I(year_contr_jaccrich^2)
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-year-quad.rds")


## 6th set of models: Year-level (intra-annual) stability ~ Year-level (intra-annual) community change, Codyn
# Richness change
form <- stabsc ~ richness_change
m1 <- inla(form, data=mammal.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-rich-year-lin.rds")

form <- stabsc ~ richness_change + I(richness_change^2)
m1 <- inla(form, data=mammal.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-rich-year-quad.rds")

# Evenness change
form <- stabsc ~ evenness_change
m1 <- inla(form, data=mammal.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-even-year-lin.rds")

form <- stabsc ~ evenness_change + I(evenness_change^2)
m1 <- inla(form, data=mammal.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-even-year-quad.rds")

# Rank change
form <- stabsc ~ rank_change
m1 <- inla(form, data=mammal.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-rank-year-lin.rds")

form <- stabsc ~ rank_change + I(rank_change^2)
m1 <- inla(form, data=mammal.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-rank-year-quad.rds")


## 7th set of models: Year-level (intra-annual) stability ~ Year-level (intra-annual) community change, BAT: site-level means
# Jaccard dissimilarity
form <- stabsc ~ year_jacc
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-year-lin-S.rds")

form <- stabsc ~ year_jacc + I(year_jacc^2)
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-year-quad-S.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ year_jaccrepl
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-year-lin-S.rds")

form <- stabsc ~ year_jaccrepl + I(year_jaccrepl^2)
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-year-quad-S.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ year_jaccrich
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-year-lin-S.rds")

form <- stabsc ~ year_jaccrich + I(year_jaccrich^2)
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-year-quad-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrepl
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-year-lin-S.rds")

form <- stabsc ~ year_contr_jaccrepl + I(year_contr_jaccrepl^2)
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-year-quad-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrich
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-year-lin-S.rds")

form <- stabsc ~ year_contr_jaccrich + I(year_contr_jaccrich^2)
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-year-quad-S.rds")


## 8th set of models: Year-level (intra-annual) stability ~ Year-level (intra-annual) community change, Codyn: site-level means
# Richness change
form <- stabsc ~ richness_change
m1 <- inla(form, data=mammal.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-rich-year-lin-S.rds")

form <- stabsc ~ richness_change + I(richness_change^2)
m1 <- inla(form, data=mammal.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-rich-year-quad-S.rds")

# Evenness change
form <- stabsc ~ evenness_change
m1 <- inla(form, data=mammal.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-even-year-lin-S.rds")

form <- stabsc ~ evenness_change + I(evenness_change^2)
m1 <- inla(form, data=mammal.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-even-year-quad-S.rds")

# Rank change
form <- stabsc ~ rank_change
m1 <- inla(form, data=mammal.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-rank-year-lin-S.rds")

form <- stabsc ~ rank_change + I(rank_change^2)
m1 <- inla(form, data=mammal.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-rank-year-quad-S.rds")


### 2. STABILITY vs LATITUDE (total of 2 models)
## Bout-level (intra-annual) stability vs latitude, Site-level means
form <- stabsc ~ latsc
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-latsc-bout-S.rds")

## Year-level (intra-annual) stability vs latitude, Site-level means
form <- stabsc ~ latsc
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-latsc-year-S.rds")



### 3. TAXONOMIC COMMUNITY CHANGE vs LATITUDE (total of 32 models)
## 1st set of models: Bout-level (intra-annual) community change vs latitude, BAT
# Jaccard dissimilarity
form <- bout_jacc ~ latsc
m1 <- inla(form, data=mammal.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jacc-latsc-bout.rds")

# Jaccard dissimilarity due to replacement
form <- bout_jaccrepl ~ latsc
m1 <- inla(form, data=mammal.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccrepl-latsc-bout.rds")

# Jaccard dissimilarity due to richness
form <- bout_jaccrich ~ latsc
m1 <- inla(form, data=mammal.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccrich-latsc-bout.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- bout_contr_jaccrepl ~ latsc
m1 <- inla(form, data=mammal.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrepl-latsc-bout.rds")

# Contribution of richness to Jaccard dissimilarity
form <- bout_contr_jaccrich ~ latsc
m1 <- inla(form, data=mammal.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrich-latsc-bout.rds")


## 2nd set of models: Bout-level (intra-annual) community change vs latitude, Codyn
# Richness change
form <- richness_change ~ latsc
m1 <- inla(form, data=mammal.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_rich-latsc-bout.rds")

# Evenness change
form <- evenness_change ~ latsc
m1 <- inla(form, data=mammal.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_even-latsc-bout.rds")

# Rank change
form <- rank_change ~ latsc
m1 <- inla(form, data=mammal.bout.codyn, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_rank-latsc-bout.rds")


## 3rd set of models: Bout-level (intra-annual) community change vs latitude, BAT: site-level means
# Jaccard dissimilarity
form <- bout_jacc ~ latsc
m1 <- inla(form, data=mammal.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jacc-latsc-bout-S.rds")

# Jaccard dissimilarity due to replacement
form <- bout_jaccrepl ~ latsc
m1 <- inla(form, data=mammal.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccrepl-latsc-bout-S.rds")

# Jaccard dissimilarity due to richness
form <- bout_jaccrich ~ latsc
m1 <- inla(form, data=mammal.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccrich-latsc-bout-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- bout_contr_jaccrepl ~ latsc
m1 <- inla(form, data=mammal.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrepl-latsc-bout-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- bout_contr_jaccrich ~ latsc
m1 <- inla(form, data=mammal.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrich-latsc-bout-S.rds")


## 4th set of models: Bout-level (intra-annual) community change vs latitude, Codyn: site-level means
# Richness change
form <- richness_change ~ latsc
m1 <- inla(form, data=mammal.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_rich-latsc-bout-S.rds")

# Evenness change
form <- evenness_change ~ latsc
m1 <- inla(form, data=mammal.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_even-latsc-bout-S.rds")

# Rank change
form <- rank_change ~ latsc
m1 <- inla(form, data=mammal.bout.codyn.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_rank-latsc-bout-S.rds")


## 5th set of models: Year-level (inter-annual) community change vs latitude, BAT
# Jaccard dissimilarity
form <- year_jacc ~ latsc
m1 <- inla(form, data=mammal.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jacc-latsc-year.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccrepl ~ latsc
m1 <- inla(form, data=mammal.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccrepl-latsc-year.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ latsc
m1 <- inla(form, data=mammal.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccrich-latsc-year.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- year_contr_jaccrepl ~ latsc
m1 <- inla(form, data=mammal.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrepl-latsc-year.rds")

# Contribution of richness to Jaccard dissimilarity
form <- year_contr_jaccrich ~ latsc
m1 <- inla(form, data=mammal.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrich-latsc-year.rds")


## 6th set of models: Year-level (inter-annual) community change vs latitude, Codyn
# Richness change
form <- richness_change ~ latsc
m1 <- inla(form, data=mammal.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_rich-latsc-year.rds")

# Evenness change
form <- evenness_change ~ latsc
m1 <- inla(form, data=mammal.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_even-latsc-year.rds")

# Rank change
form <- rank_change ~ latsc
m1 <- inla(form, data=mammal.year.codyn, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_rank-latsc-year.rds")


## 7th set of models: Year-level (inter-annual) community change vs latitude, BAT: site-level means
# Jaccard dissimilarity
form <- year_jacc ~ latsc
m1 <- inla(form, data=mammal.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jacc-latsc-year-S.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccrepl ~ latsc
m1 <- inla(form, data=mammal.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccrepl-latsc-year-S.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ latsc
m1 <- inla(form, data=mammal.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccrich-latsc-year-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- year_contr_jaccrepl ~ latsc
m1 <- inla(form, data=mammal.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrepl-latsc-year-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- year_contr_jaccrich ~ latsc
m1 <- inla(form, data=mammal.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrich-latsc-year-S.rds")


## 8th set of models: Year-level (inter-annual) community change vs latitude, Codyn: site-level means
# Richness change
form <- richness_change ~ latsc
m1 <- inla(form, data=mammal.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_rich-latsc-year-S.rds")

# Evenness change
form <- evenness_change ~ latsc
m1 <- inla(form, data=mammal.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_even-latsc-year-S.rds")

# Rank change
form <- rank_change ~ latsc
m1 <- inla(form, data=mammal.year.codyn.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_rank-latsc-year-S.rds")


### 4.YEAR-LEVEL STABILITY vs BOUT-LEVEL STABILITY (total of 2 models)
mammal.year.bat.sel <- mammal.year.bat %>%
  mutate(year_stabsc = stabsc) %>%
  select(siteID, year_stabsc) %>%
  group_by(siteID) %>%
  summarize(year_stabsc = mean(year_stabsc, na.rm=TRUE)) 

mammal.bout.year.bat <- mammal.bout.bat %>%
  left_join(mammal.year.bat.sel, by = c("siteID")) %>%
  mutate(bout_stabsc = stabsc)
mammal.bout.year.bat.m <- mammal.bout.year.bat %>%
  group_by(siteID) %>%
  summarize(bout_stabsc = mean(bout_stabsc, na.rm=TRUE), year_stabsc = mean(year_stabsc, na.rm=TRUE)) 

# All values
form <- year_stabsc ~ bout_stabsc
m1 <- inla(form, data=mammal.bout.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearstabsc-boutstabsc.rds")

# Site-level means: this is the model we want, otherwise replication of the same observations
form <-  year_stabsc ~ bout_stabsc
m1 <- inla(form, data=mammal.bout.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearstabsc-boutstabsc-S.rds")


### 5. YEAR-LEVEL TAXONOMIC COMMUNITY CHANGE vs BOUT-LEVEL TAXONOMIC COMMUNITY CHANGE (total of 12 models)
mammal.year.bat.sel <- mammal.year.bat.m %>%
  select(siteID, year_jacc, year_jaccrepl, year_jaccrich)
mammal.bout.year.bat <- mammal.bout.bat %>%
  left_join(mammal.year.bat.sel, by = c("siteID"))

mammal.bout.year.bat.m <- mammal.bout.bat.m %>%
  left_join(mammal.year.bat.sel, by = c("siteID"))

mammal.year.codyn.sel <- mammal.year.codyn.m %>%
  transmute(siteID = siteID, year_richness_change = richness_change, year_evenness_change = evenness_change, year_rank_change = rank_change) 
mammal.bout.year.codyn <- mammal.bout.codyn %>%
  left_join(mammal.year.codyn.sel, by = c("siteID")) %>% 
  mutate(bout_richness_change = richness_change, bout_evenness_change = evenness_change, bout_rank_change = rank_change) %>%
  select(-c(richness_change, evenness_change, rank_change))

mammal.bout.year.codyn.m <- mammal.bout.codyn.m %>%
  left_join(mammal.year.codyn.sel, by = c("siteID")) %>% 
  mutate(bout_richness_change = richness_change, bout_evenness_change = evenness_change, bout_rank_change = rank_change) %>%
  select(-c(richness_change, evenness_change, rank_change))

## 1st set of models: BAT
# Jaccard dissimilarity
form <- year_jacc ~ bout_jacc
m1 <- inla(form, data=mammal.bout.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearjacc-boutjacc.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccrepl ~ bout_jaccrepl
m1 <- inla(form, data=mammal.bout.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearjaccrepl-boutjaccrepl.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ bout_jaccrich
m1 <- inla(form, data=mammal.bout.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearjaccrich-boutjaccrich.rds")


## 2nd set of models: Codyn
# Richness change
form <- year_richness_change ~ bout_richness_change
m1 <- inla(form, data=mammal.bout.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearrich-boutrich.rds")

# Evenness change
form <- year_evenness_change ~ bout_evenness_change
m1 <- inla(form, data=mammal.bout.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yeareven-bouteven.rds")

# Rank change
form <- year_rank_change ~ bout_rank_change
m1 <- inla(form, data=mammal.bout.year.codyn, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearrank-boutrank.rds")


## 3rd set of models: BAT, Site-level means
# Jaccard dissimilarity
form <- year_jacc ~ bout_jacc
m1 <- inla(form, data=mammal.bout.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearjacc-boutjacc-S.rds")

# Jaccard dissimilarity due to replacement
form <-  year_jaccrepl ~ bout_jaccrepl
m1 <- inla(form, data=mammal.bout.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearjaccrepl-boutjaccrepl-S.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ bout_jaccrich
m1 <- inla(form, data=mammal.bout.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearjaccrich-boutjaccrich-S.rds")


## 4th set of models: Codyn, Site-level means
# Richness change
form <- year_richness_change ~ bout_richness_change
m1 <- inla(form, data=mammal.bout.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearrich-boutrich-S.rds")

# Evenness change
form <- year_evenness_change ~ bout_evenness_change
m1 <- inla(form, data=mammal.bout.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yeareven-bouteven-S.rds")

# Rank change
form <- year_rank_change ~ bout_rank_change
m1 <- inla(form, data=mammal.bout.year.codyn.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearrank-boutrank-S.rds")




################### ground beetles
### Read in data on stability and community compositional change
beetle.bout.stab <- readRDS(file="data/beetle_outStability-bout.rds") %>%
  select(siteID, year, bout_cv, bout_stability)
beetle.year.stab <- readRDS(file="data/beetle_outStability-year.rds") %>%
  select(siteID, year_cv, year_stability)
beetle.bout.bat <- readRDS(file="data/beetle_outBAT-bout.rds")
beetle.year.bat <- readRDS(file="data/beetle_outBAT-year.rds")
beetle.bout.codyn <- readRDS(file="data/beetle_outCodyn-bout.rds")
beetle.year.codyn <- readRDS(file="data/beetle_outCodyn-year.rds")

# quantify contributions of replacement and richness to dissimilarity 
beetle.bout.bat <- beetle.bout.bat %>% 
  mutate(bout_contr_jaccrepl = bout_jaccrepl/bout_jacc, bout_contr_jaccrich = bout_jaccrich/bout_jacc)
beetle.year.bat <- beetle.year.bat %>% 
  mutate(year_contr_jaccrepl = year_jaccrepl/year_jacc, year_contr_jaccrich = year_jaccrich/year_jacc)

# join with stability
beetle.bout.bat <- beetle.bout.bat %>%
  left_join(beetle.bout.stab, by = c("siteID","year"))
beetle.year.bat <- beetle.year.bat %>%
  left_join(beetle.year.stab, by = c("siteID"))
beetle.bout.codyn <- beetle.bout.codyn %>%
  left_join(beetle.bout.stab, by = c("siteID","year"))
beetle.year.codyn <- beetle.year.codyn %>%
  left_join(beetle.year.stab, by = c("siteID"))

# remove stability outliers
# (there are no outliers)

# obtain site-level means
beetle.bout.bat.m <- beetle.bout.bat %>%
  group_by(domainID, siteID, lat, lon, year, bout_cv, bout_stability) %>%
  summarize(bout_jacc = mean(bout_jacc, na.rm=TRUE), bout_jaccrepl = mean(bout_jaccrepl, na.rm=TRUE), bout_jaccrich = mean(bout_jaccrich, na.rm=TRUE),
            bout_contr_jaccrepl = mean(bout_contr_jaccrepl, na.rm=TRUE), bout_contr_jaccrich = mean(bout_contr_jaccrich, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, year, bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich, bout_cv, bout_stability) %>%
  ungroup

beetle.year.bat.m <- beetle.year.bat %>%
  group_by(domainID, siteID, lat, lon, year_cv, year_stability) %>%
  summarize(year_jacc = mean(year_jacc, na.rm=TRUE), year_jaccrepl = mean(year_jaccrepl, na.rm=TRUE), year_jaccrich = mean(year_jaccrich, na.rm=TRUE),
            year_contr_jaccrepl = mean(year_contr_jaccrepl, na.rm=TRUE), year_contr_jaccrich = mean(year_contr_jaccrich, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich, year_cv, year_stability) %>%
  ungroup

beetle.bout.codyn.m <- beetle.bout.codyn %>%
  group_by(domainID, siteID, lat, lon, year, bout_cv, bout_stability) %>%
  summarize(richness_change = mean(richness_change, na.rm=TRUE), evenness_change = mean(evenness_change, na.rm=TRUE), rank_change = mean(rank_change, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, year, richness_change, evenness_change, rank_change, bout_cv, bout_stability) %>%
  ungroup

beetle.year.codyn.m <- beetle.year.codyn %>%
  group_by(domainID, siteID, lat, lon, year_cv, year_stability) %>%
  summarize(richness_change = mean(richness_change, na.rm=TRUE), evenness_change = mean(evenness_change, na.rm=TRUE), rank_change = mean(rank_change, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, richness_change, evenness_change, rank_change, year_cv, year_stability) %>%
  ungroup


# scale stability, latitude, and longitude, remove NAs, replace 0 and 1 values for BAT metrics and rank change with 0.0001 and 0.999 so that models with beta distribution can be fit
beetle.bout.bat <- beetle.bout.bat %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(bout_stability)) %>%
  mutate_at(vars(bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich), ~replace(.,. == 1, 0.9999))

beetle.year.bat <- beetle.year.bat %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(year_stability)) %>%
  mutate_at(vars(year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich), ~replace(.,. == 1, 0.9999)) 

beetle.bout.codyn <- beetle.bout.codyn %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(bout_stability)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 1, 0.9999))

beetle.year.codyn <- beetle.year.codyn %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(year_stability)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 1, 0.9999))

beetle.bout.bat.m <- beetle.bout.bat.m %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(bout_stability)) %>%
  mutate_at(vars(bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich), ~replace(.,. == 1, 0.9999))

beetle.year.bat.m <- beetle.year.bat.m %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(year_stability)) %>%
  mutate_at(vars(year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich), ~replace(.,. == 1, 0.9999)) 

beetle.bout.codyn.m <- beetle.bout.codyn.m %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(bout_stability)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 1, 0.9999))

beetle.year.codyn.m <- beetle.year.codyn.m %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(year_stability)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 1, 0.9999))


##### MODELS
### 1. STABILITY vs TAXONOMIC COMMUNITY CHANGE (total of 32 models)
## 1st set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, BAT
# Jaccard dissimilarity
form <- stabsc ~ bout_jacc
m1 <- inla(form, data=beetle.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jacc-bout-lin.rds")

form <- stabsc ~ bout_jacc + I(bout_jacc^2)
m1 <- inla(form, data=beetle.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jacc-bout-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ bout_jaccrepl
m1 <- inla(form, data=beetle.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jaccrepl-bout-lin.rds")

form <- stabsc ~ bout_jaccrepl + I(bout_jaccrepl^2)
m1 <- inla(form, data=beetle.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jaccrepl-bout-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ bout_jaccrich
m1 <- inla(form, data=beetle.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jaccrich-bout-lin.rds")

form <- stabsc ~ bout_jaccrich + I(bout_jaccrich^2)
m1 <- inla(form, data=beetle.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jaccrich-bout-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrepl
m1 <- inla(form, data=beetle.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-contrjaccrepl-bout-lin.rds")

form <- stabsc ~ bout_contr_jaccrepl + I(bout_contr_jaccrepl^2)
m1 <- inla(form, data=beetle.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-contrjaccrepl-bout-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrich
m1 <- inla(form, data=beetle.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-contrjaccrich-bout-lin.rds")

form <- stabsc ~ bout_contr_jaccrich + I(bout_contr_jaccrich^2)
m1 <- inla(form, data=beetle.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-contrjaccrich-bout-quad.rds")


## 2nd set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, Codyn
# Richness change
form <- stabsc ~ richness_change
m1 <- inla(form, data=beetle.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-rich-bout-lin.rds")

form <- stabsc ~ richness_change + I(richness_change^2)
m1 <- inla(form, data=beetle.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-rich-bout-quad.rds")

# Evenness change
form <- stabsc ~ evenness_change
m1 <- inla(form, data=beetle.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-even-bout-lin.rds")

form <- stabsc ~ evenness_change + I(evenness_change^2)
m1 <- inla(form, data=beetle.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-even-bout-quad.rds")

# Rank change
form <- stabsc ~ rank_change
m1 <- inla(form, data=beetle.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-rank-bout-lin.rds")

form <- stabsc ~ rank_change + I(rank_change^2)
m1 <- inla(form, data=beetle.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-rank-bout-quad.rds")


## 3rd set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, BAT: site-level means
# Jaccard dissimilarity
form <- stabsc ~ bout_jacc
m1 <- inla(form, data=beetle.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jacc-bout-lin-S.rds")

form <- stabsc ~ bout_jacc + I(bout_jacc^2)
m1 <- inla(form, data=beetle.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jacc-bout-quad-S.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ bout_jaccrepl
m1 <- inla(form, data=beetle.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jaccrepl-bout-lin-S.rds")

form <- stabsc ~ bout_jaccrepl + I(bout_jaccrepl^2)
m1 <- inla(form, data=beetle.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jaccrepl-bout-quad-S.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ bout_jaccrich
m1 <- inla(form, data=beetle.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jaccrich-bout-lin-S.rds")

form <- stabsc ~ bout_jaccrich + I(bout_jaccrich^2)
m1 <- inla(form, data=beetle.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jaccrich-bout-quad-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrepl
m1 <- inla(form, data=beetle.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-contrjaccrepl-bout-lin-S.rds")

form <- stabsc ~ bout_contr_jaccrepl + I(bout_contr_jaccrepl^2)
m1 <- inla(form, data=beetle.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-contrjaccrepl-bout-quad-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrich
m1 <- inla(form, data=beetle.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-contrjaccrich-bout-lin-S.rds")

form <- stabsc ~ bout_contr_jaccrich + I(bout_contr_jaccrich^2)
m1 <- inla(form, data=beetle.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-contrjaccrich-bout-quad-S.rds")


## 4th set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, Codyn: site-level means
# Richness change
form <- stabsc ~ richness_change
m1 <- inla(form, data=beetle.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-rich-bout-lin-S.rds")

form <- stabsc ~ richness_change + I(richness_change^2)
m1 <- inla(form, data=beetle.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-rich-bout-quad-S.rds")

# Evenness change
form <- stabsc ~ evenness_change
m1 <- inla(form, data=beetle.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-even-bout-lin-S.rds")

form <- stabsc ~ evenness_change + I(evenness_change^2)
m1 <- inla(form, data=beetle.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-even-bout-quad-S.rds")

# Rank change
form <- stabsc ~ rank_change
m1 <- inla(form, data=beetle.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-rank-bout-lin-S.rds")

form <- stabsc ~ rank_change + I(rank_change^2)
m1 <- inla(form, data=beetle.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-rank-bout-quad-S.rds")


## 5th set of models: Year-level (inter-annual) stability ~ Year-level (inter-annual) community change, BAT
# Jaccard dissimilarity
form <- stabsc ~ year_jacc
m1 <- inla(form, data=beetle.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jacc-year-lin.rds")

form <- stabsc ~ year_jacc + I(year_jacc^2)
m1 <- inla(form, data=beetle.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jacc-year-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ year_jaccrepl
m1 <- inla(form, data=beetle.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jaccrepl-year-lin.rds")

form <- stabsc ~ year_jaccrepl + I(year_jaccrepl^2)
m1 <- inla(form, data=beetle.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jaccrepl-year-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ year_jaccrich
m1 <- inla(form, data=beetle.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jaccrich-year-lin.rds")

form <- stabsc ~ year_jaccrich + I(year_jaccrich^2)
m1 <- inla(form, data=beetle.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jaccrich-year-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrepl
m1 <- inla(form, data=beetle.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-contrjaccrepl-year-lin.rds")

form <- stabsc ~ year_contr_jaccrepl + I(year_contr_jaccrepl^2)
m1 <- inla(form, data=beetle.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-contrjaccrepl-year-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrich
m1 <- inla(form, data=beetle.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-contrjaccrich-year-lin.rds")

form <- stabsc ~ year_contr_jaccrich + I(year_contr_jaccrich^2)
m1 <- inla(form, data=beetle.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-contrjaccrich-year-quad.rds")


## 6th set of models: Year-level (intra-annual) stability ~ Year-level (intra-annual) community change, Codyn
# Richness change
form <- stabsc ~ richness_change
m1 <- inla(form, data=beetle.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-rich-year-lin.rds")

form <- stabsc ~ richness_change + I(richness_change^2)
m1 <- inla(form, data=beetle.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-rich-year-quad.rds")

# Evenness change
form <- stabsc ~ evenness_change
m1 <- inla(form, data=beetle.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-even-year-lin.rds")

form <- stabsc ~ evenness_change + I(evenness_change^2)
m1 <- inla(form, data=beetle.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-even-year-quad.rds")

# Rank change
form <- stabsc ~ rank_change
m1 <- inla(form, data=beetle.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-rank-year-lin.rds")

form <- stabsc ~ rank_change + I(rank_change^2)
m1 <- inla(form, data=beetle.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-rank-year-quad.rds")


## 7th set of models: Year-level (intra-annual) stability ~ Year-level (intra-annual) community change, BAT: site-level means
# Jaccard dissimilarity
form <- stabsc ~ year_jacc
m1 <- inla(form, data=beetle.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jacc-year-lin-S.rds")

form <- stabsc ~ year_jacc + I(year_jacc^2)
m1 <- inla(form, data=beetle.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jacc-year-quad-S.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ year_jaccrepl
m1 <- inla(form, data=beetle.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jaccrepl-year-lin-S.rds")

form <- stabsc ~ year_jaccrepl + I(year_jaccrepl^2)
m1 <- inla(form, data=beetle.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jaccrepl-year-quad-S.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ year_jaccrich
m1 <- inla(form, data=beetle.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jaccrich-year-lin-S.rds")

form <- stabsc ~ year_jaccrich + I(year_jaccrich^2)
m1 <- inla(form, data=beetle.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-jaccrich-year-quad-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrepl
m1 <- inla(form, data=beetle.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-contrjaccrepl-year-lin-S.rds")

form <- stabsc ~ year_contr_jaccrepl + I(year_contr_jaccrepl^2)
m1 <- inla(form, data=beetle.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-contrjaccrepl-year-quad-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrich
m1 <- inla(form, data=beetle.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-contrjaccrich-year-lin-S.rds")

form <- stabsc ~ year_contr_jaccrich + I(year_contr_jaccrich^2)
m1 <- inla(form, data=beetle.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-contrjaccrich-year-quad-S.rds")


## 8th set of models: Year-level (intra-annual) stability ~ Year-level (intra-annual) community change, Codyn: site-level means
# Richness change
form <- stabsc ~ richness_change
m1 <- inla(form, data=beetle.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-rich-year-lin-S.rds")

form <- stabsc ~ richness_change + I(richness_change^2)
m1 <- inla(form, data=beetle.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-rich-year-quad-S.rds")

# Evenness change
form <- stabsc ~ evenness_change
m1 <- inla(form, data=beetle.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-even-year-lin-S.rds")

form <- stabsc ~ evenness_change + I(evenness_change^2)
m1 <- inla(form, data=beetle.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-even-year-quad-S.rds")

# Rank change
form <- stabsc ~ rank_change
m1 <- inla(form, data=beetle.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-rank-year-lin-S.rds")

form <- stabsc ~ rank_change + I(rank_change^2)
m1 <- inla(form, data=beetle.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-rank-year-quad-S.rds")


### 2. STABILITY vs LATITUDE (total of 2 models)
## Bout-level (intra-annual) stability vs latitude, Site-level means
form <- stabsc ~ latsc
m1 <- inla(form, data=beetle.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-latsc-bout-S.rds")

## Year-level (intra-annual) stability vs latitude, Site-level means
form <- stabsc ~ latsc
m1 <- inla(form, data=beetle.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_stabsc-latsc-year-S.rds")



### 3. TAXONOMIC COMMUNITY CHANGE vs LATITUDE (total of 32 models)
## 1st set of models: Bout-level (intra-annual) community change vs latitude, BAT
# Jaccard dissimilarity
form <- bout_jacc ~ latsc
m1 <- inla(form, data=beetle.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_jacc-latsc-bout.rds")

# Jaccard dissimilarity due to replacement
form <- bout_jaccrepl ~ latsc
m1 <- inla(form, data=beetle.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_jaccrepl-latsc-bout.rds")

# Jaccard dissimilarity due to richness
form <- bout_jaccrich ~ latsc
m1 <- inla(form, data=beetle.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_jaccrich-latsc-bout.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- bout_contr_jaccrepl ~ latsc
m1 <- inla(form, data=beetle.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_contrjaccrepl-latsc-bout.rds")

# Contribution of richness to Jaccard dissimilarity
form <- bout_contr_jaccrich ~ latsc
m1 <- inla(form, data=beetle.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_contrjaccrich-latsc-bout.rds")


## 2nd set of models: Bout-level (intra-annual) community change vs latitude, Codyn
# Richness change
form <- richness_change ~ latsc
m1 <- inla(form, data=beetle.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_rich-latsc-bout.rds")

# Evenness change
form <- evenness_change ~ latsc
m1 <- inla(form, data=beetle.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_even-latsc-bout.rds")

# Rank change
form <- rank_change ~ latsc
m1 <- inla(form, data=beetle.bout.codyn, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_rank-latsc-bout.rds")


## 3rd set of models: Bout-level (intra-annual) community change vs latitude, BAT: site-level means
# Jaccard dissimilarity
form <- bout_jacc ~ latsc
m1 <- inla(form, data=beetle.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_jacc-latsc-bout-S.rds")

# Jaccard dissimilarity due to replacement
form <- bout_jaccrepl ~ latsc
m1 <- inla(form, data=beetle.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_jaccrepl-latsc-bout-S.rds")

# Jaccard dissimilarity due to richness
form <- bout_jaccrich ~ latsc
m1 <- inla(form, data=beetle.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_jaccrich-latsc-bout-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- bout_contr_jaccrepl ~ latsc
m1 <- inla(form, data=beetle.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_contrjaccrepl-latsc-bout-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- bout_contr_jaccrich ~ latsc
m1 <- inla(form, data=beetle.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_contrjaccrich-latsc-bout-S.rds")


## 4th set of models: Bout-level (intra-annual) community change vs latitude, Codyn: site-level means
# Richness change
form <- richness_change ~ latsc
m1 <- inla(form, data=beetle.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_rich-latsc-bout-S.rds")

# Evenness change
form <- evenness_change ~ latsc
m1 <- inla(form, data=beetle.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_even-latsc-bout-S.rds")

# Rank change
form <- rank_change ~ latsc
m1 <- inla(form, data=beetle.bout.codyn.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_rank-latsc-bout-S.rds")


## 5th set of models: Year-level (inter-annual) community change vs latitude, BAT
# Jaccard dissimilarity
form <- year_jacc ~ latsc
m1 <- inla(form, data=beetle.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_jacc-latsc-year.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccrepl ~ latsc
m1 <- inla(form, data=beetle.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_jaccrepl-latsc-year.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ latsc
m1 <- inla(form, data=beetle.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_jaccrich-latsc-year.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- year_contr_jaccrepl ~ latsc
m1 <- inla(form, data=beetle.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_contrjaccrepl-latsc-year.rds")

# Contribution of richness to Jaccard dissimilarity
form <- year_contr_jaccrich ~ latsc
m1 <- inla(form, data=beetle.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_contrjaccrich-latsc-year.rds")


## 6th set of models: Year-level (inter-annual) community change vs latitude, Codyn
# Richness change
form <- richness_change ~ latsc
m1 <- inla(form, data=beetle.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_rich-latsc-year.rds")

# Evenness change
form <- evenness_change ~ latsc
m1 <- inla(form, data=beetle.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_even-latsc-year.rds")

# Rank change
form <- rank_change ~ latsc
m1 <- inla(form, data=beetle.year.codyn, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_rank-latsc-year.rds")


## 7th set of models: Year-level (inter-annual) community change vs latitude, BAT: site-level means
# Jaccard dissimilarity
form <- year_jacc ~ latsc
m1 <- inla(form, data=beetle.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_jacc-latsc-year-S.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccrepl ~ latsc
m1 <- inla(form, data=beetle.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_jaccrepl-latsc-year-S.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ latsc
m1 <- inla(form, data=beetle.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_jaccrich-latsc-year-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- year_contr_jaccrepl ~ latsc
m1 <- inla(form, data=beetle.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_contrjaccrepl-latsc-year-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- year_contr_jaccrich ~ latsc
m1 <- inla(form, data=beetle.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_contrjaccrich-latsc-year-S.rds")


## 8th set of models: Year-level (inter-annual) community change vs latitude, Codyn: site-level means
# Richness change
form <- richness_change ~ latsc
m1 <- inla(form, data=beetle.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_rich-latsc-year-S.rds")

# Evenness change
form <- evenness_change ~ latsc
m1 <- inla(form, data=beetle.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_even-latsc-year-S.rds")

# Rank change
form <- rank_change ~ latsc
m1 <- inla(form, data=beetle.year.codyn.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_rank-latsc-year-S.rds")


### 4.YEAR-LEVEL STABILITY vs BOUT-LEVEL STABILITY (total of 2 models)
beetle.year.bat.sel <- beetle.year.bat %>%
  mutate(year_stabsc = stabsc) %>%
  select(siteID, year_stabsc) %>%
  group_by(siteID) %>%
  summarize(year_stabsc = mean(year_stabsc, na.rm=TRUE)) 

beetle.bout.year.bat <- beetle.bout.bat %>%
  left_join(beetle.year.bat.sel, by = c("siteID")) %>%
  mutate(bout_stabsc = stabsc)
beetle.bout.year.bat.m <- beetle.bout.year.bat %>%
  group_by(siteID) %>%
  summarize(bout_stabsc = mean(bout_stabsc, na.rm=TRUE), year_stabsc = mean(year_stabsc, na.rm=TRUE)) 

# All values
form <- year_stabsc ~ bout_stabsc
m1 <- inla(form, data=beetle.bout.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_yearstabsc-boutstabsc.rds")

# Site-level means: this is the model we want, otherwise replication of the same observations
form <-  year_stabsc ~ bout_stabsc
m1 <- inla(form, data=beetle.bout.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_yearstabsc-boutstabsc-S.rds")


### 5. YEAR-LEVEL TAXONOMIC COMMUNITY CHANGE vs BOUT-LEVEL TAXONOMIC COMMUNITY CHANGE (total of 12 models)
beetle.year.bat.sel <- beetle.year.bat.m %>%
  select(siteID, year_jacc, year_jaccrepl, year_jaccrich)
beetle.bout.year.bat <- beetle.bout.bat %>%
  left_join(beetle.year.bat.sel, by = c("siteID"))

beetle.bout.year.bat.m <- beetle.bout.bat.m %>%
  left_join(beetle.year.bat.sel, by = c("siteID"))

beetle.year.codyn.sel <- beetle.year.codyn.m %>%
  transmute(siteID = siteID, year_richness_change = richness_change, year_evenness_change = evenness_change, year_rank_change = rank_change) 
beetle.bout.year.codyn <- beetle.bout.codyn %>%
  left_join(beetle.year.codyn.sel, by = c("siteID")) %>% 
  mutate(bout_richness_change = richness_change, bout_evenness_change = evenness_change, bout_rank_change = rank_change) %>%
  select(-c(richness_change, evenness_change, rank_change))

beetle.bout.year.codyn.m <- beetle.bout.codyn.m %>%
  left_join(beetle.year.codyn.sel, by = c("siteID")) %>% 
  mutate(bout_richness_change = richness_change, bout_evenness_change = evenness_change, bout_rank_change = rank_change) %>%
  select(-c(richness_change, evenness_change, rank_change))

## 1st set of models: BAT
# Jaccard dissimilarity
form <- year_jacc ~ bout_jacc
m1 <- inla(form, data=beetle.bout.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_yearjacc-boutjacc.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccrepl ~ bout_jaccrepl
m1 <- inla(form, data=beetle.bout.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_yearjaccrepl-boutjaccrepl.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ bout_jaccrich
m1 <- inla(form, data=beetle.bout.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_yearjaccrich-boutjaccrich.rds")


## 2nd set of models: Codyn
# Richness change
form <- year_richness_change ~ bout_richness_change
m1 <- inla(form, data=beetle.bout.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_yearrich-boutrich.rds")

# Evenness change
form <- year_evenness_change ~ bout_evenness_change
m1 <- inla(form, data=beetle.bout.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_yeareven-bouteven.rds")

# Rank change
form <- year_rank_change ~ bout_rank_change
m1 <- inla(form, data=beetle.bout.year.codyn, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_yearrank-boutrank.rds")


## 3rd set of models: BAT, Site-level means
# Jaccard dissimilarity
form <- year_jacc ~ bout_jacc
m1 <- inla(form, data=beetle.bout.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_yearjacc-boutjacc-S.rds")

# Jaccard dissimilarity due to replacement
form <-  year_jaccrepl ~ bout_jaccrepl
m1 <- inla(form, data=beetle.bout.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_yearjaccrepl-boutjaccrepl-S.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ bout_jaccrich
m1 <- inla(form, data=beetle.bout.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_yearjaccrich-boutjaccrich-S.rds")


## 4th set of models: Codyn, Site-level means
# Richness change
form <- year_richness_change ~ bout_richness_change
m1 <- inla(form, data=beetle.bout.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_yearrich-boutrich-S.rds")

# Evenness change
form <- year_evenness_change ~ bout_evenness_change
m1 <- inla(form, data=beetle.bout.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_yeareven-bouteven-S.rds")

# Rank change
form <- year_rank_change ~ bout_rank_change
m1 <- inla(form, data=beetle.bout.year.codyn.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/beetle_yearrank-boutrank-S.rds")








################### fish
### Read in data on stability and community compositional change
fish.bout.stab <- readRDS(file="data/fish_outStability-bout.rds") %>%
  select(siteID, year, bout_cv, bout_stability)
fish.year.stab <- readRDS(file="data/fish_outStability-year.rds") %>%
  select(siteID, year_cv, year_stability)
fish.bout.bat <- readRDS(file="data/fish_outBAT-bout.rds")
fish.year.bat <- readRDS(file="data/fish_outBAT-year.rds")
fish.bout.codyn <- readRDS(file="data/fish_outCodyn-bout.rds")
fish.year.codyn <- readRDS(file="data/fish_outCodyn-year.rds")

# quantify contributions of replacement and richness to dissimilarity 
fish.bout.bat <- fish.bout.bat %>% 
  mutate(bout_contr_jaccrepl = bout_jaccrepl/bout_jacc, bout_contr_jaccrich = bout_jaccrich/bout_jacc)
fish.year.bat <- fish.year.bat %>% 
  mutate(year_contr_jaccrepl = year_jaccrepl/year_jacc, year_contr_jaccrich = year_jaccrich/year_jacc)

# join with stability
fish.bout.bat <- fish.bout.bat %>%
  left_join(fish.bout.stab, by = c("siteID","year"))
fish.year.bat <- fish.year.bat %>%
  left_join(fish.year.stab, by = c("siteID"))
fish.bout.codyn <- fish.bout.codyn %>%
  left_join(fish.bout.stab, by = c("siteID","year"))
fish.year.codyn <- fish.year.codyn %>%
  left_join(fish.year.stab, by = c("siteID"))

# remove stability outliers
fish.bout.bat <- fish.bout.bat %>%
  filter(bout_stability < 30)
fish.bout.codyn <- fish.bout.codyn %>%
  filter(bout_stability < 30)

# obtain site-level means
fish.bout.bat.m <- fish.bout.bat %>%
  group_by(domainID, siteID, lat, lon, year, bout_cv, bout_stability) %>%
  summarize(bout_jacc = mean(bout_jacc, na.rm=TRUE), bout_jaccrepl = mean(bout_jaccrepl, na.rm=TRUE), bout_jaccrich = mean(bout_jaccrich, na.rm=TRUE),
            bout_contr_jaccrepl = mean(bout_contr_jaccrepl, na.rm=TRUE), bout_contr_jaccrich = mean(bout_contr_jaccrich, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, year, bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich, bout_cv, bout_stability) %>%
  ungroup

fish.year.bat.m <- fish.year.bat %>%
  group_by(domainID, siteID, lat, lon, year_cv, year_stability) %>%
  summarize(year_jacc = mean(year_jacc, na.rm=TRUE), year_jaccrepl = mean(year_jaccrepl, na.rm=TRUE), year_jaccrich = mean(year_jaccrich, na.rm=TRUE),
            year_contr_jaccrepl = mean(year_contr_jaccrepl, na.rm=TRUE), year_contr_jaccrich = mean(year_contr_jaccrich, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich, year_cv, year_stability) %>%
  ungroup

fish.bout.codyn.m <- fish.bout.codyn %>%
  group_by(domainID, siteID, lat, lon, year, bout_cv, bout_stability) %>%
  summarize(richness_change = mean(richness_change, na.rm=TRUE), evenness_change = mean(evenness_change, na.rm=TRUE), rank_change = mean(rank_change, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, year, richness_change, evenness_change, rank_change, bout_cv, bout_stability) %>%
  ungroup

fish.year.codyn.m <- fish.year.codyn %>%
  group_by(domainID, siteID, lat, lon, year_cv, year_stability) %>%
  summarize(richness_change = mean(richness_change, na.rm=TRUE), evenness_change = mean(evenness_change, na.rm=TRUE), rank_change = mean(rank_change, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, richness_change, evenness_change, rank_change, year_cv, year_stability) %>%
  ungroup


# scale stability, latitude, and longitude, remove NAs, replace 0 and 1 values for BAT metrics and rank change with 0.0001 and 0.999 so that models with beta distribution can be fit
fish.bout.bat <- fish.bout.bat %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(bout_stability)) %>%
  mutate_at(vars(bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich), ~replace(.,. == 1, 0.9999))

fish.year.bat <- fish.year.bat %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(year_stability)) %>%
  mutate_at(vars(year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich), ~replace(.,. == 1, 0.9999)) 

fish.bout.codyn <- fish.bout.codyn %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(bout_stability)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 1, 0.9999))

fish.year.codyn <- fish.year.codyn %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(year_stability)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 1, 0.9999))

fish.bout.bat.m <- fish.bout.bat.m %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(bout_stability)) %>%
  mutate_at(vars(bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich), ~replace(.,. == 1, 0.9999))

fish.year.bat.m <- fish.year.bat.m %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(year_stability)) %>%
  mutate_at(vars(year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich), ~replace(.,. == 1, 0.9999)) 

fish.bout.codyn.m <- fish.bout.codyn.m %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(bout_stability)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 1, 0.9999))

fish.year.codyn.m <- fish.year.codyn.m %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(year_stability)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 1, 0.9999))


##### MODELS
### 1. STABILITY vs TAXONOMIC COMMUNITY CHANGE (total of 32 models)
## 1st set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, BAT
# Jaccard dissimilarity
form <- stabsc ~ bout_jacc
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-bout-lin.rds")

form <- stabsc ~ bout_jacc + I(bout_jacc^2)
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-bout-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ bout_jaccrepl
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-bout-lin.rds")

form <- stabsc ~ bout_jaccrepl + I(bout_jaccrepl^2)
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-bout-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ bout_jaccrich
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-bout-lin.rds")

form <- stabsc ~ bout_jaccrich + I(bout_jaccrich^2)
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-bout-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrepl
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-bout-lin.rds")

form <- stabsc ~ bout_contr_jaccrepl + I(bout_contr_jaccrepl^2)
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-bout-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrich
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-bout-lin.rds")

form <- stabsc ~ bout_contr_jaccrich + I(bout_contr_jaccrich^2)
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-bout-quad.rds")


## 2nd set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, Codyn
# Richness change
form <- stabsc ~ richness_change
m1 <- inla(form, data=fish.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-rich-bout-lin.rds")

form <- stabsc ~ richness_change + I(richness_change^2)
m1 <- inla(form, data=fish.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-rich-bout-quad.rds")

# Evenness change
form <- stabsc ~ evenness_change
m1 <- inla(form, data=fish.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-even-bout-lin.rds")

form <- stabsc ~ evenness_change + I(evenness_change^2)
m1 <- inla(form, data=fish.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-even-bout-quad.rds")

# Rank change
form <- stabsc ~ rank_change
m1 <- inla(form, data=fish.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-rank-bout-lin.rds")

form <- stabsc ~ rank_change + I(rank_change^2)
m1 <- inla(form, data=fish.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-rank-bout-quad.rds")


## 3rd set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, BAT: site-level means
# Jaccard dissimilarity
form <- stabsc ~ bout_jacc
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-bout-lin-S.rds")

form <- stabsc ~ bout_jacc + I(bout_jacc^2)
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-bout-quad-S.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ bout_jaccrepl
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-bout-lin-S.rds")

form <- stabsc ~ bout_jaccrepl + I(bout_jaccrepl^2)
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-bout-quad-S.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ bout_jaccrich
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-bout-lin-S.rds")

form <- stabsc ~ bout_jaccrich + I(bout_jaccrich^2)
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-bout-quad-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrepl
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-bout-lin-S.rds")

form <- stabsc ~ bout_contr_jaccrepl + I(bout_contr_jaccrepl^2)
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-bout-quad-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrich
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-bout-lin-S.rds")

form <- stabsc ~ bout_contr_jaccrich + I(bout_contr_jaccrich^2)
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-bout-quad-S.rds")


## 4th set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, Codyn: site-level means
# Richness change
form <- stabsc ~ richness_change
m1 <- inla(form, data=fish.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-rich-bout-lin-S.rds")

form <- stabsc ~ richness_change + I(richness_change^2)
m1 <- inla(form, data=fish.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-rich-bout-quad-S.rds")

# Evenness change
form <- stabsc ~ evenness_change
m1 <- inla(form, data=fish.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-even-bout-lin-S.rds")

form <- stabsc ~ evenness_change + I(evenness_change^2)
m1 <- inla(form, data=fish.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-even-bout-quad-S.rds")

# Rank change
form <- stabsc ~ rank_change
m1 <- inla(form, data=fish.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-rank-bout-lin-S.rds")

form <- stabsc ~ rank_change + I(rank_change^2)
m1 <- inla(form, data=fish.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-rank-bout-quad-S.rds")


## 5th set of models: Year-level (inter-annual) stability ~ Year-level (inter-annual) community change, BAT
# Jaccard dissimilarity
form <- stabsc ~ year_jacc
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-year-lin.rds")

form <- stabsc ~ year_jacc + I(year_jacc^2)
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-year-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ year_jaccrepl
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-year-lin.rds")

form <- stabsc ~ year_jaccrepl + I(year_jaccrepl^2)
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-year-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ year_jaccrich
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-year-lin.rds")

form <- stabsc ~ year_jaccrich + I(year_jaccrich^2)
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-year-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrepl
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-year-lin.rds")

form <- stabsc ~ year_contr_jaccrepl + I(year_contr_jaccrepl^2)
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-year-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrich
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-year-lin.rds")

form <- stabsc ~ year_contr_jaccrich + I(year_contr_jaccrich^2)
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-year-quad.rds")


## 6th set of models: Year-level (intra-annual) stability ~ Year-level (intra-annual) community change, Codyn
# Richness change
form <- stabsc ~ richness_change
m1 <- inla(form, data=fish.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-rich-year-lin.rds")

form <- stabsc ~ richness_change + I(richness_change^2)
m1 <- inla(form, data=fish.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-rich-year-quad.rds")

# Evenness change
form <- stabsc ~ evenness_change
m1 <- inla(form, data=fish.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-even-year-lin.rds")

form <- stabsc ~ evenness_change + I(evenness_change^2)
m1 <- inla(form, data=fish.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-even-year-quad.rds")

# Rank change
form <- stabsc ~ rank_change
m1 <- inla(form, data=fish.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-rank-year-lin.rds")

form <- stabsc ~ rank_change + I(rank_change^2)
m1 <- inla(form, data=fish.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-rank-year-quad.rds")


## 7th set of models: Year-level (intra-annual) stability ~ Year-level (intra-annual) community change, BAT: site-level means
# Jaccard dissimilarity
form <- stabsc ~ year_jacc
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-year-lin-S.rds")

form <- stabsc ~ year_jacc + I(year_jacc^2)
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-year-quad-S.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ year_jaccrepl
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-year-lin-S.rds")

form <- stabsc ~ year_jaccrepl + I(year_jaccrepl^2)
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-year-quad-S.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ year_jaccrich
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-year-lin-S.rds")

form <- stabsc ~ year_jaccrich + I(year_jaccrich^2)
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-year-quad-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrepl
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-year-lin-S.rds")

form <- stabsc ~ year_contr_jaccrepl + I(year_contr_jaccrepl^2)
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-year-quad-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrich
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-year-lin-S.rds")

form <- stabsc ~ year_contr_jaccrich + I(year_contr_jaccrich^2)
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-year-quad-S.rds")


## 8th set of models: Year-level (intra-annual) stability ~ Year-level (intra-annual) community change, Codyn: site-level means
# Richness change
form <- stabsc ~ richness_change
m1 <- inla(form, data=fish.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-rich-year-lin-S.rds")

form <- stabsc ~ richness_change + I(richness_change^2)
m1 <- inla(form, data=fish.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-rich-year-quad-S.rds")

# Evenness change
form <- stabsc ~ evenness_change
m1 <- inla(form, data=fish.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-even-year-lin-S.rds")

form <- stabsc ~ evenness_change + I(evenness_change^2)
m1 <- inla(form, data=fish.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-even-year-quad-S.rds")

# Rank change
form <- stabsc ~ rank_change
m1 <- inla(form, data=fish.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-rank-year-lin-S.rds")

form <- stabsc ~ rank_change + I(rank_change^2)
m1 <- inla(form, data=fish.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-rank-year-quad-S.rds")


### 2. STABILITY vs LATITUDE (total of 2 models)
## Bout-level (intra-annual) stability vs latitude, Site-level means
form <- stabsc ~ latsc
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-latsc-bout-S.rds")

## Year-level (intra-annual) stability vs latitude, Site-level means
form <- stabsc ~ latsc
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-latsc-year-S.rds")



### 3. TAXONOMIC COMMUNITY CHANGE vs LATITUDE (total of 32 models)
## 1st set of models: Bout-level (intra-annual) community change vs latitude, BAT
# Jaccard dissimilarity
form <- bout_jacc ~ latsc
m1 <- inla(form, data=fish.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jacc-latsc-bout.rds")

# Jaccard dissimilarity due to replacement
form <- bout_jaccrepl ~ latsc
m1 <- inla(form, data=fish.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccrepl-latsc-bout.rds")

# Jaccard dissimilarity due to richness
form <- bout_jaccrich ~ latsc
m1 <- inla(form, data=fish.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccrich-latsc-bout.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- bout_contr_jaccrepl ~ latsc
m1 <- inla(form, data=fish.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrepl-latsc-bout.rds")

# Contribution of richness to Jaccard dissimilarity
form <- bout_contr_jaccrich ~ latsc
m1 <- inla(form, data=fish.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrich-latsc-bout.rds")


## 2nd set of models: Bout-level (intra-annual) community change vs latitude, Codyn
# Richness change
form <- richness_change ~ latsc
m1 <- inla(form, data=fish.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_rich-latsc-bout.rds")

# Evenness change
form <- evenness_change ~ latsc
m1 <- inla(form, data=fish.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_even-latsc-bout.rds")

# Rank change
form <- rank_change ~ latsc
m1 <- inla(form, data=fish.bout.codyn, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_rank-latsc-bout.rds")


## 3rd set of models: Bout-level (intra-annual) community change vs latitude, BAT: site-level means
# Jaccard dissimilarity
form <- bout_jacc ~ latsc
m1 <- inla(form, data=fish.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jacc-latsc-bout-S.rds")

# Jaccard dissimilarity due to replacement
form <- bout_jaccrepl ~ latsc
m1 <- inla(form, data=fish.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccrepl-latsc-bout-S.rds")

# Jaccard dissimilarity due to richness
form <- bout_jaccrich ~ latsc
m1 <- inla(form, data=fish.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccrich-latsc-bout-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- bout_contr_jaccrepl ~ latsc
m1 <- inla(form, data=fish.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrepl-latsc-bout-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- bout_contr_jaccrich ~ latsc
m1 <- inla(form, data=fish.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrich-latsc-bout-S.rds")


## 4th set of models: Bout-level (intra-annual) community change vs latitude, Codyn: site-level means
# Richness change
form <- richness_change ~ latsc
m1 <- inla(form, data=fish.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_rich-latsc-bout-S.rds")

# Evenness change
form <- evenness_change ~ latsc
m1 <- inla(form, data=fish.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_even-latsc-bout-S.rds")

# Rank change
form <- rank_change ~ latsc
m1 <- inla(form, data=fish.bout.codyn.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_rank-latsc-bout-S.rds")


## 5th set of models: Year-level (inter-annual) community change vs latitude, BAT
# Jaccard dissimilarity
form <- year_jacc ~ latsc
m1 <- inla(form, data=fish.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jacc-latsc-year.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccrepl ~ latsc
m1 <- inla(form, data=fish.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccrepl-latsc-year.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ latsc
m1 <- inla(form, data=fish.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccrich-latsc-year.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- year_contr_jaccrepl ~ latsc
m1 <- inla(form, data=fish.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrepl-latsc-year.rds")

# Contribution of richness to Jaccard dissimilarity
form <- year_contr_jaccrich ~ latsc
m1 <- inla(form, data=fish.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrich-latsc-year.rds")


## 6th set of models: Year-level (inter-annual) community change vs latitude, Codyn
# Richness change
form <- richness_change ~ latsc
m1 <- inla(form, data=fish.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_rich-latsc-year.rds")

# Evenness change
form <- evenness_change ~ latsc
m1 <- inla(form, data=fish.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_even-latsc-year.rds")

# Rank change
form <- rank_change ~ latsc
m1 <- inla(form, data=fish.year.codyn, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_rank-latsc-year.rds")


## 7th set of models: Year-level (inter-annual) community change vs latitude, BAT: site-level means
# Jaccard dissimilarity
form <- year_jacc ~ latsc
m1 <- inla(form, data=fish.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jacc-latsc-year-S.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccrepl ~ latsc
m1 <- inla(form, data=fish.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccrepl-latsc-year-S.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ latsc
m1 <- inla(form, data=fish.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccrich-latsc-year-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- year_contr_jaccrepl ~ latsc
m1 <- inla(form, data=fish.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrepl-latsc-year-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- year_contr_jaccrich ~ latsc
m1 <- inla(form, data=fish.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrich-latsc-year-S.rds")


## 8th set of models: Year-level (inter-annual) community change vs latitude, Codyn: site-level means
# Richness change
form <- richness_change ~ latsc
m1 <- inla(form, data=fish.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_rich-latsc-year-S.rds")

# Evenness change
form <- evenness_change ~ latsc
m1 <- inla(form, data=fish.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_even-latsc-year-S.rds")

# Rank change
form <- rank_change ~ latsc
m1 <- inla(form, data=fish.year.codyn.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_rank-latsc-year-S.rds")


### 4.YEAR-LEVEL STABILITY vs BOUT-LEVEL STABILITY (total of 2 models)
fish.year.bat.sel <- fish.year.bat %>%
  mutate(year_stabsc = stabsc) %>%
  select(siteID, year_stabsc) %>%
  group_by(siteID) %>%
  summarize(year_stabsc = mean(year_stabsc, na.rm=TRUE)) 

fish.bout.year.bat <- fish.bout.bat %>%
  left_join(fish.year.bat.sel, by = c("siteID")) %>%
  mutate(bout_stabsc = stabsc)
fish.bout.year.bat.m <- fish.bout.year.bat %>%
  group_by(siteID) %>%
  summarize(bout_stabsc = mean(bout_stabsc, na.rm=TRUE), year_stabsc = mean(year_stabsc, na.rm=TRUE)) 

# All values
form <- year_stabsc ~ bout_stabsc
m1 <- inla(form, data=fish.bout.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearstabsc-boutstabsc.rds")

# Site-level means: this is the model we want, otherwise replication of the same observations
form <-  year_stabsc ~ bout_stabsc
m1 <- inla(form, data=fish.bout.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearstabsc-boutstabsc-S.rds")


### 5. YEAR-LEVEL TAXONOMIC COMMUNITY CHANGE vs BOUT-LEVEL TAXONOMIC COMMUNITY CHANGE (total of 12 models)
fish.year.bat.sel <- fish.year.bat.m %>%
  select(siteID, year_jacc, year_jaccrepl, year_jaccrich)
fish.bout.year.bat <- fish.bout.bat %>%
  left_join(fish.year.bat.sel, by = c("siteID"))

fish.bout.year.bat.m <- fish.bout.bat.m %>%
  left_join(fish.year.bat.sel, by = c("siteID"))

fish.year.codyn.sel <- fish.year.codyn.m %>%
  transmute(siteID = siteID, year_richness_change = richness_change, year_evenness_change = evenness_change, year_rank_change = rank_change) 
fish.bout.year.codyn <- fish.bout.codyn %>%
  left_join(fish.year.codyn.sel, by = c("siteID")) %>% 
  mutate(bout_richness_change = richness_change, bout_evenness_change = evenness_change, bout_rank_change = rank_change) %>%
  select(-c(richness_change, evenness_change, rank_change))

fish.bout.year.codyn.m <- fish.bout.codyn.m %>%
  left_join(fish.year.codyn.sel, by = c("siteID")) %>% 
  mutate(bout_richness_change = richness_change, bout_evenness_change = evenness_change, bout_rank_change = rank_change) %>%
  select(-c(richness_change, evenness_change, rank_change))

## 1st set of models: BAT
# Jaccard dissimilarity
form <- year_jacc ~ bout_jacc
m1 <- inla(form, data=fish.bout.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearjacc-boutjacc.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccrepl ~ bout_jaccrepl
m1 <- inla(form, data=fish.bout.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearjaccrepl-boutjaccrepl.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ bout_jaccrich
m1 <- inla(form, data=fish.bout.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearjaccrich-boutjaccrich.rds")


## 2nd set of models: Codyn
# Richness change
form <- year_richness_change ~ bout_richness_change
m1 <- inla(form, data=fish.bout.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearrich-boutrich.rds")

# Evenness change
form <- year_evenness_change ~ bout_evenness_change
m1 <- inla(form, data=fish.bout.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yeareven-bouteven.rds")

# Rank change
form <- year_rank_change ~ bout_rank_change
m1 <- inla(form, data=fish.bout.year.codyn, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearrank-boutrank.rds")


## 3rd set of models: BAT, Site-level means
# Jaccard dissimilarity
form <- year_jacc ~ bout_jacc
m1 <- inla(form, data=fish.bout.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearjacc-boutjacc-S.rds")

# Jaccard dissimilarity due to replacement
form <-  year_jaccrepl ~ bout_jaccrepl
m1 <- inla(form, data=fish.bout.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearjaccrepl-boutjaccrepl-S.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ bout_jaccrich
m1 <- inla(form, data=fish.bout.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearjaccrich-boutjaccrich-S.rds")


## 4th set of models: Codyn, Site-level means
# Richness change
form <- year_richness_change ~ bout_richness_change
m1 <- inla(form, data=fish.bout.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearrich-boutrich-S.rds")

# Evenness change
form <- year_evenness_change ~ bout_evenness_change
m1 <- inla(form, data=fish.bout.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yeareven-bouteven-S.rds")

# Rank change
form <- year_rank_change ~ bout_rank_change
m1 <- inla(form, data=fish.bout.year.codyn.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearrank-boutrank-S.rds")







################### aquatic macroinvertebrates
### Read in data on stability and community compositional change
macroinv.bout.stab <- readRDS(file="data/macroinv_outStability-bout.rds") %>%
  select(siteID, year, bout_cv, bout_stability)
macroinv.year.stab <- readRDS(file="data/macroinv_outStability-year.rds") %>%
  select(siteID, year_cv, year_stability)
macroinv.bout.bat <- readRDS(file="data/macroinv_outBAT-bout.rds")
macroinv.year.bat <- readRDS(file="data/macroinv_outBAT-year.rds")
macroinv.bout.codyn <- readRDS(file="data/macroinv_outCodyn-bout.rds")
macroinv.year.codyn <- readRDS(file="data/macroinv_outCodyn-year.rds")

# quantify contributions of replacement and richness to dissimilarity 
macroinv.bout.bat <- macroinv.bout.bat %>% 
  mutate(bout_contr_jaccrepl = bout_jaccrepl/bout_jacc, bout_contr_jaccrich = bout_jaccrich/bout_jacc)
macroinv.year.bat <- macroinv.year.bat %>% 
  mutate(year_contr_jaccrepl = year_jaccrepl/year_jacc, year_contr_jaccrich = year_jaccrich/year_jacc)

# join with stability
macroinv.bout.bat <- macroinv.bout.bat %>%
  left_join(macroinv.bout.stab, by = c("siteID","year"))
macroinv.year.bat <- macroinv.year.bat %>%
  left_join(macroinv.year.stab, by = c("siteID"))
macroinv.bout.codyn <- macroinv.bout.codyn %>%
  left_join(macroinv.bout.stab, by = c("siteID","year"))
macroinv.year.codyn <- macroinv.year.codyn %>%
  left_join(macroinv.year.stab, by = c("siteID"))

# remove stability outliers
macroinv.bout.bat <- macroinv.bout.bat %>%
  filter(bout_stability < 50)
macroinv.bout.codyn <- macroinv.bout.codyn %>%
  filter(bout_stability < 50)
macroinv.year.bat <- macroinv.year.bat %>%
  filter(year_stability < 15)
macroinv.year.codyn <- macroinv.year.codyn %>%
  filter(year_stability < 15)

# obtain site-level means
macroinv.bout.bat.m <- macroinv.bout.bat %>%
  group_by(domainID, siteID, lat, lon, year, bout_cv, bout_stability) %>%
  summarize(bout_jacc = mean(bout_jacc, na.rm=TRUE), bout_jaccrepl = mean(bout_jaccrepl, na.rm=TRUE), bout_jaccrich = mean(bout_jaccrich, na.rm=TRUE),
            bout_contr_jaccrepl = mean(bout_contr_jaccrepl, na.rm=TRUE), bout_contr_jaccrich = mean(bout_contr_jaccrich, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, year, bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich, bout_cv, bout_stability) %>%
  ungroup

macroinv.year.bat.m <- macroinv.year.bat %>%
  group_by(domainID, siteID, lat, lon, year_cv, year_stability) %>%
  summarize(year_jacc = mean(year_jacc, na.rm=TRUE), year_jaccrepl = mean(year_jaccrepl, na.rm=TRUE), year_jaccrich = mean(year_jaccrich, na.rm=TRUE),
            year_contr_jaccrepl = mean(year_contr_jaccrepl, na.rm=TRUE), year_contr_jaccrich = mean(year_contr_jaccrich, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich, year_cv, year_stability) %>%
  ungroup

macroinv.bout.codyn.m <- macroinv.bout.codyn %>%
  group_by(domainID, siteID, lat, lon, year, bout_cv, bout_stability) %>%
  summarize(richness_change = mean(richness_change, na.rm=TRUE), evenness_change = mean(evenness_change, na.rm=TRUE), rank_change = mean(rank_change, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, year, richness_change, evenness_change, rank_change, bout_cv, bout_stability) %>%
  ungroup

macroinv.year.codyn.m <- macroinv.year.codyn %>%
  group_by(domainID, siteID, lat, lon, year_cv, year_stability) %>%
  summarize(richness_change = mean(richness_change, na.rm=TRUE), evenness_change = mean(evenness_change, na.rm=TRUE), rank_change = mean(rank_change, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, richness_change, evenness_change, rank_change, year_cv, year_stability) %>%
  ungroup


# scale stability, latitude, and longitude, remove NAs, replace 0 and 1 values for BAT metrics and rank change with 0.0001 and 0.999 so that models with beta distribution can be fit
macroinv.bout.bat <- macroinv.bout.bat %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(bout_stability)) %>%
  mutate_at(vars(bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich), ~replace(.,. == 1, 0.9999))

macroinv.year.bat <- macroinv.year.bat %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(year_stability)) %>%
  mutate_at(vars(year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich), ~replace(.,. == 1, 0.9999)) 

macroinv.bout.codyn <- macroinv.bout.codyn %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(bout_stability)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 1, 0.9999))

macroinv.year.codyn <- macroinv.year.codyn %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(year_stability)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 1, 0.9999))

macroinv.bout.bat.m <- macroinv.bout.bat.m %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(bout_stability)) %>%
  mutate_at(vars(bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich), ~replace(.,. == 1, 0.9999))

macroinv.year.bat.m <- macroinv.year.bat.m %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(year_stability)) %>%
  mutate_at(vars(year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich), ~replace(.,. == 1, 0.9999)) 

macroinv.bout.codyn.m <- macroinv.bout.codyn.m %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(bout_stability)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 1, 0.9999))

macroinv.year.codyn.m <- macroinv.year.codyn.m %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(year_stability)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 0, 0.0001)) %>%
  mutate_at(vars(rank_change), ~replace(.,. == 1, 0.9999))

##### MODELS
### 1. STABILITY vs TAXONOMIC COMMUNITY CHANGE (total of 32 models)
## 1st set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, BAT
# Jaccard dissimilarity
form <- stabsc ~ bout_jacc
m1 <- inla(form, data=macroinv.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jacc-bout-lin.rds")

form <- stabsc ~ bout_jacc + I(bout_jacc^2)
m1 <- inla(form, data=macroinv.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jacc-bout-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ bout_jaccrepl
m1 <- inla(form, data=macroinv.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jaccrepl-bout-lin.rds")

form <- stabsc ~ bout_jaccrepl + I(bout_jaccrepl^2)
m1 <- inla(form, data=macroinv.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jaccrepl-bout-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ bout_jaccrich
m1 <- inla(form, data=macroinv.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jaccrich-bout-lin.rds")

form <- stabsc ~ bout_jaccrich + I(bout_jaccrich^2)
m1 <- inla(form, data=macroinv.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jaccrich-bout-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrepl
m1 <- inla(form, data=macroinv.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-contrjaccrepl-bout-lin.rds")

form <- stabsc ~ bout_contr_jaccrepl + I(bout_contr_jaccrepl^2)
m1 <- inla(form, data=macroinv.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-contrjaccrepl-bout-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrich
m1 <- inla(form, data=macroinv.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-contrjaccrich-bout-lin.rds")

form <- stabsc ~ bout_contr_jaccrich + I(bout_contr_jaccrich^2)
m1 <- inla(form, data=macroinv.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-contrjaccrich-bout-quad.rds")


## 2nd set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, Codyn
# Richness change
form <- stabsc ~ richness_change
m1 <- inla(form, data=macroinv.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-rich-bout-lin.rds")

form <- stabsc ~ richness_change + I(richness_change^2)
m1 <- inla(form, data=macroinv.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-rich-bout-quad.rds")

# Evenness change
form <- stabsc ~ evenness_change
m1 <- inla(form, data=macroinv.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-even-bout-lin.rds")

form <- stabsc ~ evenness_change + I(evenness_change^2)
m1 <- inla(form, data=macroinv.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-even-bout-quad.rds")

# Rank change
form <- stabsc ~ rank_change
m1 <- inla(form, data=macroinv.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-rank-bout-lin.rds")

form <- stabsc ~ rank_change + I(rank_change^2)
m1 <- inla(form, data=macroinv.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-rank-bout-quad.rds")


## 3rd set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, BAT: site-level means
# Jaccard dissimilarity
form <- stabsc ~ bout_jacc
m1 <- inla(form, data=macroinv.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jacc-bout-lin-S.rds")

form <- stabsc ~ bout_jacc + I(bout_jacc^2)
m1 <- inla(form, data=macroinv.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jacc-bout-quad-S.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ bout_jaccrepl
m1 <- inla(form, data=macroinv.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jaccrepl-bout-lin-S.rds")

form <- stabsc ~ bout_jaccrepl + I(bout_jaccrepl^2)
m1 <- inla(form, data=macroinv.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jaccrepl-bout-quad-S.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ bout_jaccrich
m1 <- inla(form, data=macroinv.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jaccrich-bout-lin-S.rds")

form <- stabsc ~ bout_jaccrich + I(bout_jaccrich^2)
m1 <- inla(form, data=macroinv.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jaccrich-bout-quad-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrepl
m1 <- inla(form, data=macroinv.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-contrjaccrepl-bout-lin-S.rds")

form <- stabsc ~ bout_contr_jaccrepl + I(bout_contr_jaccrepl^2)
m1 <- inla(form, data=macroinv.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-contrjaccrepl-bout-quad-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrich
m1 <- inla(form, data=macroinv.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-contrjaccrich-bout-lin-S.rds")

form <- stabsc ~ bout_contr_jaccrich + I(bout_contr_jaccrich^2)
m1 <- inla(form, data=macroinv.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-contrjaccrich-bout-quad-S.rds")


## 4th set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, Codyn: site-level means
# Richness change
form <- stabsc ~ richness_change
m1 <- inla(form, data=macroinv.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-rich-bout-lin-S.rds")

form <- stabsc ~ richness_change + I(richness_change^2)
m1 <- inla(form, data=macroinv.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-rich-bout-quad-S.rds")

# Evenness change
form <- stabsc ~ evenness_change
m1 <- inla(form, data=macroinv.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-even-bout-lin-S.rds")

form <- stabsc ~ evenness_change + I(evenness_change^2)
m1 <- inla(form, data=macroinv.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-even-bout-quad-S.rds")

# Rank change
form <- stabsc ~ rank_change
m1 <- inla(form, data=macroinv.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-rank-bout-lin-S.rds")

form <- stabsc ~ rank_change + I(rank_change^2)
m1 <- inla(form, data=macroinv.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-rank-bout-quad-S.rds")


## 5th set of models: Year-level (inter-annual) stability ~ Year-level (inter-annual) community change, BAT
# Jaccard dissimilarity
form <- stabsc ~ year_jacc
m1 <- inla(form, data=macroinv.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jacc-year-lin.rds")

form <- stabsc ~ year_jacc + I(year_jacc^2)
m1 <- inla(form, data=macroinv.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jacc-year-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ year_jaccrepl
m1 <- inla(form, data=macroinv.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jaccrepl-year-lin.rds")

form <- stabsc ~ year_jaccrepl + I(year_jaccrepl^2)
m1 <- inla(form, data=macroinv.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jaccrepl-year-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ year_jaccrich
m1 <- inla(form, data=macroinv.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jaccrich-year-lin.rds")

form <- stabsc ~ year_jaccrich + I(year_jaccrich^2)
m1 <- inla(form, data=macroinv.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jaccrich-year-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrepl
m1 <- inla(form, data=macroinv.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-contrjaccrepl-year-lin.rds")

form <- stabsc ~ year_contr_jaccrepl + I(year_contr_jaccrepl^2)
m1 <- inla(form, data=macroinv.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-contrjaccrepl-year-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrich
m1 <- inla(form, data=macroinv.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-contrjaccrich-year-lin.rds")

form <- stabsc ~ year_contr_jaccrich + I(year_contr_jaccrich^2)
m1 <- inla(form, data=macroinv.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-contrjaccrich-year-quad.rds")


## 6th set of models: Year-level (intra-annual) stability ~ Year-level (intra-annual) community change, Codyn
# Richness change
form <- stabsc ~ richness_change
m1 <- inla(form, data=macroinv.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-rich-year-lin.rds")

form <- stabsc ~ richness_change + I(richness_change^2)
m1 <- inla(form, data=macroinv.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-rich-year-quad.rds")

# Evenness change
form <- stabsc ~ evenness_change
m1 <- inla(form, data=macroinv.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-even-year-lin.rds")

form <- stabsc ~ evenness_change + I(evenness_change^2)
m1 <- inla(form, data=macroinv.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-even-year-quad.rds")

# Rank change
form <- stabsc ~ rank_change
m1 <- inla(form, data=macroinv.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-rank-year-lin.rds")

form <- stabsc ~ rank_change + I(rank_change^2)
m1 <- inla(form, data=macroinv.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-rank-year-quad.rds")


## 7th set of models: Year-level (intra-annual) stability ~ Year-level (intra-annual) community change, BAT: site-level means
# Jaccard dissimilarity
form <- stabsc ~ year_jacc
m1 <- inla(form, data=macroinv.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jacc-year-lin-S.rds")

form <- stabsc ~ year_jacc + I(year_jacc^2)
m1 <- inla(form, data=macroinv.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jacc-year-quad-S.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ year_jaccrepl
m1 <- inla(form, data=macroinv.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jaccrepl-year-lin-S.rds")

form <- stabsc ~ year_jaccrepl + I(year_jaccrepl^2)
m1 <- inla(form, data=macroinv.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jaccrepl-year-quad-S.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ year_jaccrich
m1 <- inla(form, data=macroinv.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jaccrich-year-lin-S.rds")

form <- stabsc ~ year_jaccrich + I(year_jaccrich^2)
m1 <- inla(form, data=macroinv.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-jaccrich-year-quad-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrepl
m1 <- inla(form, data=macroinv.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-contrjaccrepl-year-lin-S.rds")

form <- stabsc ~ year_contr_jaccrepl + I(year_contr_jaccrepl^2)
m1 <- inla(form, data=macroinv.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-contrjaccrepl-year-quad-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrich
m1 <- inla(form, data=macroinv.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-contrjaccrich-year-lin-S.rds")

form <- stabsc ~ year_contr_jaccrich + I(year_contr_jaccrich^2)
m1 <- inla(form, data=macroinv.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-contrjaccrich-year-quad-S.rds")


## 8th set of models: Year-level (intra-annual) stability ~ Year-level (intra-annual) community change, Codyn: site-level means
# Richness change
form <- stabsc ~ richness_change
m1 <- inla(form, data=macroinv.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-rich-year-lin-S.rds")

form <- stabsc ~ richness_change + I(richness_change^2)
m1 <- inla(form, data=macroinv.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-rich-year-quad-S.rds")

# Evenness change
form <- stabsc ~ evenness_change
m1 <- inla(form, data=macroinv.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-even-year-lin-S.rds")

form <- stabsc ~ evenness_change + I(evenness_change^2)
m1 <- inla(form, data=macroinv.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-even-year-quad-S.rds")

# Rank change
form <- stabsc ~ rank_change
m1 <- inla(form, data=macroinv.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-rank-year-lin-S.rds")

form <- stabsc ~ rank_change + I(rank_change^2)
m1 <- inla(form, data=macroinv.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-rank-year-quad-S.rds")


### 2. STABILITY vs LATITUDE (total of 2 models)
## Bout-level (intra-annual) stability vs latitude, Site-level means
form <- stabsc ~ latsc
m1 <- inla(form, data=macroinv.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-latsc-bout-S.rds")

## Year-level (intra-annual) stability vs latitude, Site-level means
form <- stabsc ~ latsc
m1 <- inla(form, data=macroinv.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_stabsc-latsc-year-S.rds")



### 3. TAXONOMIC COMMUNITY CHANGE vs LATITUDE (total of 32 models)
## 1st set of models: Bout-level (intra-annual) community change vs latitude, BAT
# Jaccard dissimilarity
form <- bout_jacc ~ latsc
m1 <- inla(form, data=macroinv.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_jacc-latsc-bout.rds")

# Jaccard dissimilarity due to replacement
form <- bout_jaccrepl ~ latsc
m1 <- inla(form, data=macroinv.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_jaccrepl-latsc-bout.rds")

# Jaccard dissimilarity due to richness
form <- bout_jaccrich ~ latsc
m1 <- inla(form, data=macroinv.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_jaccrich-latsc-bout.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- bout_contr_jaccrepl ~ latsc
m1 <- inla(form, data=macroinv.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_contrjaccrepl-latsc-bout.rds")

# Contribution of richness to Jaccard dissimilarity
form <- bout_contr_jaccrich ~ latsc
m1 <- inla(form, data=macroinv.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_contrjaccrich-latsc-bout.rds")


## 2nd set of models: Bout-level (intra-annual) community change vs latitude, Codyn
# Richness change
form <- richness_change ~ latsc
m1 <- inla(form, data=macroinv.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_rich-latsc-bout.rds")

# Evenness change
form <- evenness_change ~ latsc
m1 <- inla(form, data=macroinv.bout.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_even-latsc-bout.rds")

# Rank change
form <- rank_change ~ latsc
m1 <- inla(form, data=macroinv.bout.codyn, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_rank-latsc-bout.rds")


## 3rd set of models: Bout-level (intra-annual) community change vs latitude, BAT: site-level means
# Jaccard dissimilarity
form <- bout_jacc ~ latsc
m1 <- inla(form, data=macroinv.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_jacc-latsc-bout-S.rds")

# Jaccard dissimilarity due to replacement
form <- bout_jaccrepl ~ latsc
m1 <- inla(form, data=macroinv.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_jaccrepl-latsc-bout-S.rds")

# Jaccard dissimilarity due to richness
form <- bout_jaccrich ~ latsc
m1 <- inla(form, data=macroinv.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_jaccrich-latsc-bout-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- bout_contr_jaccrepl ~ latsc
m1 <- inla(form, data=macroinv.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_contrjaccrepl-latsc-bout-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- bout_contr_jaccrich ~ latsc
m1 <- inla(form, data=macroinv.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_contrjaccrich-latsc-bout-S.rds")


## 4th set of models: Bout-level (intra-annual) community change vs latitude, Codyn: site-level means
# Richness change
form <- richness_change ~ latsc
m1 <- inla(form, data=macroinv.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_rich-latsc-bout-S.rds")

# Evenness change
form <- evenness_change ~ latsc
m1 <- inla(form, data=macroinv.bout.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_even-latsc-bout-S.rds")

# Rank change
form <- rank_change ~ latsc
m1 <- inla(form, data=macroinv.bout.codyn.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_rank-latsc-bout-S.rds")


## 5th set of models: Year-level (inter-annual) community change vs latitude, BAT
# Jaccard dissimilarity
form <- year_jacc ~ latsc
m1 <- inla(form, data=macroinv.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_jacc-latsc-year.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccrepl ~ latsc
m1 <- inla(form, data=macroinv.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_jaccrepl-latsc-year.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ latsc
m1 <- inla(form, data=macroinv.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_jaccrich-latsc-year.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- year_contr_jaccrepl ~ latsc
m1 <- inla(form, data=macroinv.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_contrjaccrepl-latsc-year.rds")

# Contribution of richness to Jaccard dissimilarity
form <- year_contr_jaccrich ~ latsc
m1 <- inla(form, data=macroinv.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_contrjaccrich-latsc-year.rds")


## 6th set of models: Year-level (inter-annual) community change vs latitude, Codyn
# Richness change
form <- richness_change ~ latsc
m1 <- inla(form, data=macroinv.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_rich-latsc-year.rds")

# Evenness change
form <- evenness_change ~ latsc
m1 <- inla(form, data=macroinv.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_even-latsc-year.rds")

# Rank change
form <- rank_change ~ latsc
m1 <- inla(form, data=macroinv.year.codyn, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_rank-latsc-year.rds")


## 7th set of models: Year-level (inter-annual) community change vs latitude, BAT: site-level means
# Jaccard dissimilarity
form <- year_jacc ~ latsc
m1 <- inla(form, data=macroinv.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_jacc-latsc-year-S.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccrepl ~ latsc
m1 <- inla(form, data=macroinv.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_jaccrepl-latsc-year-S.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ latsc
m1 <- inla(form, data=macroinv.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_jaccrich-latsc-year-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- year_contr_jaccrepl ~ latsc
m1 <- inla(form, data=macroinv.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_contrjaccrepl-latsc-year-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- year_contr_jaccrich ~ latsc
m1 <- inla(form, data=macroinv.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_contrjaccrich-latsc-year-S.rds")


## 8th set of models: Year-level (inter-annual) community change vs latitude, Codyn: site-level means
# Richness change
form <- richness_change ~ latsc
m1 <- inla(form, data=macroinv.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_rich-latsc-year-S.rds")

# Evenness change
form <- evenness_change ~ latsc
m1 <- inla(form, data=macroinv.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_even-latsc-year-S.rds")

# Rank change
form <- rank_change ~ latsc
m1 <- inla(form, data=macroinv.year.codyn.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_rank-latsc-year-S.rds")


### 4.YEAR-LEVEL STABILITY vs BOUT-LEVEL STABILITY (total of 2 models)
macroinv.year.bat.sel <- macroinv.year.bat %>%
  mutate(year_stabsc = stabsc) %>%
  select(siteID, year_stabsc) %>%
  group_by(siteID) %>%
  summarize(year_stabsc = mean(year_stabsc, na.rm=TRUE)) 

macroinv.bout.year.bat <- macroinv.bout.bat %>%
  left_join(macroinv.year.bat.sel, by = c("siteID")) %>%
  mutate(bout_stabsc = stabsc)
macroinv.bout.year.bat.m <- macroinv.bout.year.bat %>%
  group_by(siteID) %>%
  summarize(bout_stabsc = mean(bout_stabsc, na.rm=TRUE), year_stabsc = mean(year_stabsc, na.rm=TRUE)) 

# All values
form <- year_stabsc ~ bout_stabsc
m1 <- inla(form, data=macroinv.bout.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_yearstabsc-boutstabsc.rds")

# Site-level means: this is the model we want, otherwise replication of the same observations
form <-  year_stabsc ~ bout_stabsc
m1 <- inla(form, data=macroinv.bout.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_yearstabsc-boutstabsc-S.rds")


### 5. YEAR-LEVEL TAXONOMIC COMMUNITY CHANGE vs BOUT-LEVEL TAXONOMIC COMMUNITY CHANGE (total of 12 models)
macroinv.year.bat.sel <- macroinv.year.bat.m %>%
  select(siteID, year_jacc, year_jaccrepl, year_jaccrich)
macroinv.bout.year.bat <- macroinv.bout.bat %>%
  left_join(macroinv.year.bat.sel, by = c("siteID"))

macroinv.bout.year.bat.m <- macroinv.bout.bat.m %>%
  left_join(macroinv.year.bat.sel, by = c("siteID"))

macroinv.year.codyn.sel <- macroinv.year.codyn.m %>%
  transmute(siteID = siteID, year_richness_change = richness_change, year_evenness_change = evenness_change, year_rank_change = rank_change) 
macroinv.bout.year.codyn <- macroinv.bout.codyn %>%
  left_join(macroinv.year.codyn.sel, by = c("siteID")) %>% 
  mutate(bout_richness_change = richness_change, bout_evenness_change = evenness_change, bout_rank_change = rank_change) %>%
  select(-c(richness_change, evenness_change, rank_change))

macroinv.bout.year.codyn.m <- macroinv.bout.codyn.m %>%
  left_join(macroinv.year.codyn.sel, by = c("siteID")) %>% 
  mutate(bout_richness_change = richness_change, bout_evenness_change = evenness_change, bout_rank_change = rank_change) %>%
  select(-c(richness_change, evenness_change, rank_change))

## 1st set of models: BAT
# Jaccard dissimilarity
form <- year_jacc ~ bout_jacc
m1 <- inla(form, data=macroinv.bout.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_yearjacc-boutjacc.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccrepl ~ bout_jaccrepl
m1 <- inla(form, data=macroinv.bout.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_yearjaccrepl-boutjaccrepl.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ bout_jaccrich
m1 <- inla(form, data=macroinv.bout.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_yearjaccrich-boutjaccrich.rds")


## 2nd set of models: Codyn
# Richness change
form <- year_richness_change ~ bout_richness_change
m1 <- inla(form, data=macroinv.bout.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_yearrich-boutrich.rds")

# Evenness change
form <- year_evenness_change ~ bout_evenness_change
m1 <- inla(form, data=macroinv.bout.year.codyn, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_yeareven-bouteven.rds")

# Rank change
form <- year_rank_change ~ bout_rank_change
m1 <- inla(form, data=macroinv.bout.year.codyn, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_yearrank-boutrank.rds")


## 3rd set of models: BAT, Site-level means
# Jaccard dissimilarity
form <- year_jacc ~ bout_jacc
m1 <- inla(form, data=macroinv.bout.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_yearjacc-boutjacc-S.rds")

# Jaccard dissimilarity due to replacement
form <-  year_jaccrepl ~ bout_jaccrepl
m1 <- inla(form, data=macroinv.bout.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_yearjaccrepl-boutjaccrepl-S.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ bout_jaccrich
m1 <- inla(form, data=macroinv.bout.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_yearjaccrich-boutjaccrich-S.rds")


## 4th set of models: Codyn, Site-level means
# Richness change
form <- year_richness_change ~ bout_richness_change
m1 <- inla(form, data=macroinv.bout.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_yearrich-boutrich-S.rds")

# Evenness change
form <- year_evenness_change ~ bout_evenness_change
m1 <- inla(form, data=macroinv.bout.year.codyn.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_yeareven-bouteven-S.rds")

# Rank change
form <- year_rank_change ~ bout_rank_change
m1 <- inla(form, data=macroinv.bout.year.codyn.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/macroinv_yearrank-boutrank-S.rds")





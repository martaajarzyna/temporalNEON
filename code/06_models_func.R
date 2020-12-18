#Load packages
require(rhdf5)
require(here)
require(tidyverse)
select <- dplyr::select
summarize <- dplyr::summarize
require(INLA)

###################MODELS EXPLAINING ECOSYSTEM STABILITY AS A FUNCTION OF CHANGE IN COMMUNITY FUNCTIONAL COMPOSITION: ABSOLUTE

################### small mammals
### Read in data on stability and community compositional change
mammal.bout.stab <- readRDS(file="data/mammal_outStability-bout.rds") %>%
  select(siteID, year, bout_cv, bout_stability)
mammal.year.stab <- readRDS(file="data/mammal_outStability-year.rds") %>%
  select(siteID, year_cv, year_stability)
mammal.bout.bat <- readRDS(file="data/mammal_outBAT-FD-bout.rds")
mammal.year.bat <- readRDS(file="data/mammal_outBAT-FD-year.rds")

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

# remove stability outliers
mammal.bout.bat <- mammal.bout.bat %>%
  filter(bout_stability < 10)
mammal.year.bat <- mammal.year.bat %>%
  filter(year_stability < 5)

# obtain site-level means
mammal.bout.bat.m <- mammal.bout.bat %>%
  group_by(domainID, siteID, lat, lon, year, bout_cv, bout_stability) %>%
  summarize(bout_jacc = mean(bout_jacc, na.rm=TRUE), bout_jaccrepl = mean(bout_jaccrepl, na.rm=TRUE), bout_jaccrich = mean(bout_jaccrich, na.rm=TRUE),
            bout_contr_jaccrepl = mean(bout_contr_jaccrepl, na.rm=TRUE), bout_contr_jaccrich = mean(bout_contr_jaccrich, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, year, bout_jacc, bout_jaccrepl, bout_jaccrich, bout_contr_jaccrepl, bout_contr_jaccrich, bout_cv, bout_stability) %>%
  ungroup

mammal.year.bat.m <- mammal.year.bat %>%
  group_by(domainID, siteID, lat, lon, year_cv, year_stability) %>%
  summarize(year_jacc = mean(year_jacc, na.rm=TRUE), year_jaccrepl = mean(year_jaccrepl, na.rm=TRUE), year_jaccrich = mean(year_jaccrich, na.rm=TRUE),
            year_contr_jaccrepl = mean(year_contr_jaccrepl, na.rm=TRUE), year_contr_jaccrich = mean(year_contr_jaccrich, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, year_jacc, year_jaccrepl, year_jaccrich, year_contr_jaccrepl, year_contr_jaccrich, year_cv, year_stability) %>%
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



##### MODELS
### 1. STABILITY vs FUNCTIONAL COMMUNITY CHANGE (total of 20 models)
## 1st set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, BAT
# Jaccard dissimilarity
form <- stabsc ~ bout_jacc
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-FD-bout-lin.rds")

form <- stabsc ~ bout_jacc + I(bout_jacc^2)
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-FD-bout-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ bout_jaccrepl
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-FD-bout-lin.rds")

form <- stabsc ~ bout_jaccrepl + I(bout_jaccrepl^2)
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-FD-bout-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ bout_jaccrich
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-FD-bout-lin.rds")

form <- stabsc ~ bout_jaccrich + I(bout_jaccrich^2)
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-FD-bout-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrepl
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-FD-bout-lin.rds")

form <- stabsc ~ bout_contr_jaccrepl + I(bout_contr_jaccrepl^2)
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-FD-bout-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrich
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-FD-bout-lin.rds")

form <- stabsc ~ bout_contr_jaccrich + I(bout_contr_jaccrich^2)
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-FD-bout-quad.rds")


## 2nd set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, BAT: site-level means
# Jaccard dissimilarity
form <- stabsc ~ bout_jacc
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-FD-bout-S-lin.rds")

form <- stabsc ~ bout_jacc + I(bout_jacc^2)
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-FD-bout-S-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ bout_jaccrepl
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-FD-bout-S-lin.rds")

form <- stabsc ~ bout_jaccrepl + I(bout_jaccrepl^2)
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-FD-bout-S-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ bout_jaccrich
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-FD-bout-S-lin.rds")

form <- stabsc ~ bout_jaccrich + I(bout_jaccrich^2)
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-FD-bout-S-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrepl
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-FD-bout-S-lin.rds")

form <- stabsc ~ bout_contr_jaccrepl + I(bout_contr_jaccrepl^2)
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-FD-bout-S-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrich
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-FD-bout-S-lin.rds")

form <- stabsc ~ bout_contr_jaccrich + I(bout_contr_jaccrich^2)
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-FD-bout-S-quad.rds")


## 3rd set of models: Year-level (inter-annual) stability ~ Year-level (inter-annual) community change, BAT
# Jaccard dissimilarity
form <- stabsc ~ year_jacc
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-FD-year-lin.rds")

form <- stabsc ~ year_jacc + I(year_jacc^2)
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-FD-year-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ year_jaccrepl
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-FD-year-lin.rds")

form <- stabsc ~ year_jaccrepl + I(year_jaccrepl^2)
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-FD-year-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ year_jaccrich
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-FD-year-lin.rds")

form <- stabsc ~ year_jaccrich + I(year_jaccrich^2)
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-FD-year-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrepl
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-FD-year-lin.rds")

form <- stabsc ~ year_contr_jaccrepl + I(year_contr_jaccrepl^2)
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-FD-year-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrich
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-FD-year-lin.rds")

form <- stabsc ~ year_contr_jaccrich + I(year_contr_jaccrich^2)
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-FD-year-quad.rds")

## 4th set of models: Year-level (inter-annual) stability ~ Year-level (inter-annual) community change, BAT: site-level means
# Jaccard dissimilarity
form <- stabsc ~ year_jacc
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-FD-year-S-lin.rds")

form <- stabsc ~ year_jacc + I(year_jacc^2)
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-FD-year-S-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ year_jaccrepl
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-FD-year-S-lin.rds")

form <- stabsc ~ year_jaccrepl + I(year_jaccrepl^2)
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-FD-year-S-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ year_jaccrich
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-FD-year-S-lin.rds")

form <- stabsc ~ year_jaccrich + I(year_jaccrich^2)
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-FD-year-S-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrepl
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-FD-year-S-lin.rds")

form <- stabsc ~ year_contr_jaccrepl + I(year_contr_jaccrepl^2)
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-FD-year-S-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrich
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-FD-year-S-lin.rds")

form <- stabsc ~ year_contr_jaccrich + I(year_contr_jaccrich^2)
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-FD-year-S-quad.rds")


### 2. FUNCTIONAL COMMUNITY CHANGE vs LATITUDE (total of 2- models)
## 1st set of models: Bout-level (intra-annual) community change vs latitude, BAT
# Jaccard dissimilarity
form <- bout_jacc ~ latsc
m1 <- inla(form, data=mammal.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jacc-latsc-FD-bout.rds")

# Jaccard dissimilarity due to replacement
form <- bout_jaccrepl ~ latsc
m1 <- inla(form, data=mammal.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccrepl-latsc-FD-bout.rds")

# Jaccard dissimilarity due to richness
form <- bout_jaccrich ~ latsc
m1 <- inla(form, data=mammal.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccrich-latsc-FD-bout.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- bout_contr_jaccrepl ~ latsc
m1 <- inla(form, data=mammal.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrepl-latsc-FD-bout.rds")

# Contribution of richness to Jaccard dissimilarity
form <- bout_contr_jaccrich ~ latsc
m1 <- inla(form, data=mammal.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrich-latsc-FD-bout.rds")


## 2nd set of models: Bout-level (intra-annual) community change vs latitude, BAT: site-level means
# Jaccard dissimilarity
form <- bout_jacc ~ latsc
m1 <- inla(form, data=mammal.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jacc-latsc-FD-bout-S.rds")

# Jaccard dissimilarity due to replacement
form <- bout_jaccrepl ~ latsc
m1 <- inla(form, data=mammal.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccrepl-latsc-FD-bout-S.rds")

# Jaccard dissimilarity due to richness
form <- bout_jaccrich ~ latsc
m1 <- inla(form, data=mammal.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccrich-latsc-FD-bout-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- bout_contr_jaccrepl ~ latsc
m1 <- inla(form, data=mammal.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrepl-latsc-FD-bout-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- bout_contr_jaccrich ~ latsc
m1 <- inla(form, data=mammal.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrich-latsc-FD-bout-S.rds")


## 3rd set of models: Year-level (inter-annual) community change vs latitude, BAT
# Jaccard dissimilarity
form <- year_jacc ~ latsc
m1 <- inla(form, data=mammal.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jacc-latsc-FD-year.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccrepl ~ latsc
m1 <- inla(form, data=mammal.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccrepl-latsc-FD-year.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ latsc
m1 <- inla(form, data=mammal.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccrich-latsc-FD-year.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- year_contr_jaccrepl ~ latsc
m1 <- inla(form, data=mammal.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrepl-latsc-FD-year.rds")

# Contribution of richness to Jaccard dissimilarity
form <- year_contr_jaccrich ~ latsc
m1 <- inla(form, data=mammal.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrich-latsc-FD-year.rds")


## 4th set of models: Year-level (inter-annual) community change vs latitude, BAT: site-level means
# Jaccard dissimilarity
form <- year_jacc ~ latsc
m1 <- inla(form, data=mammal.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jacc-latsc-FD-year-S.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccrepl ~ latsc
m1 <- inla(form, data=mammal.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccrepl-latsc-FD-year-S.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ latsc
m1 <- inla(form, data=mammal.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccrich-latsc-FD-year-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- year_contr_jaccrepl ~ latsc
m1 <- inla(form, data=mammal.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrepl-latsc-FD-year-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- year_contr_jaccrich ~ latsc
m1 <- inla(form, data=mammal.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrich-latsc-FD-year-S.rds")



### 3. YEAR-LEVEL FUNCTIONAL COMMUNITY CHANGE vs BOUT-LEVEL TAXONOMIC COMMUNITY CHANGE (total of 6 models)
mammal.year.bat.sel <- mammal.year.bat.m %>%
  select(siteID, year_jacc, year_jaccrepl, year_jaccrich)
mammal.bout.year.bat <- mammal.bout.bat %>%
  left_join(mammal.year.bat.sel, by = c("siteID"))

mammal.bout.year.bat.m <- mammal.bout.bat.m %>%
  left_join(mammal.year.bat.sel, by = c("siteID"))

## 1st set of models: BAT
# Jaccard dissimilarity
form <- year_jacc ~bout_jacc
m1 <- inla(form, data=mammal.bout.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearjacc-boutjacc-FD.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccrepl ~ bout_jaccrepl
m1 <- inla(form, data=mammal.bout.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearjaccrepl-boutjaccrepl-FD.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ bout_jaccrich
m1 <- inla(form, data=mammal.bout.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearjaccrich-boutjaccrich-FD.rds")


## 2nd set of models: BAT, Site-level means
# Jaccard dissimilarity
form <- year_jacc ~ bout_jacc
m1 <- inla(form, data=mammal.bout.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearjacc-boutjacc-FD-S.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccrepl ~ bout_jaccrepl
m1 <- inla(form, data=mammal.bout.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearjaccrepl-boutjaccrepl-FD-S.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ bout_jaccrich
m1 <- inla(form, data=mammal.bout.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearjaccrich-boutjaccrich-FD-S.rds")






################### fish

### Read in data on stability and community compositional change
fish.bout.stab <- readRDS(file="data/fish_outStability-bout.rds") %>%
  select(siteID, year, bout_cv, bout_stability)
fish.year.stab <- readRDS(file="data/fish_outStability-year.rds") %>%
  select(siteID, year_cv, year_stability)
fish.bout.bat <- readRDS(file="data/fish_outBAT-FD-bout.rds")
fish.year.bat <- readRDS(file="data/fish_outBAT-FD-year.rds")

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

# remove stability outliers
fish.bout.bat <- fish.bout.bat %>%
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




##### MODELS
### 1. STABILITY vs FUNCTIONAL COMMUNITY CHANGE (total of 20 models)
## 1st set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, BAT
# Jaccard dissimilarity
form <- stabsc ~ bout_jacc
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-FD-bout-lin.rds")

form <- stabsc ~ bout_jacc + I(bout_jacc^2)
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-FD-bout-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ bout_jaccrepl
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-FD-bout-lin.rds")

form <- stabsc ~ bout_jaccrepl + I(bout_jaccrepl^2)
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-FD-bout-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ bout_jaccrich
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-FD-bout-lin.rds")

form <- stabsc ~ bout_jaccrich + I(bout_jaccrich^2)
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-FD-bout-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrepl
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-FD-bout-lin.rds")

form <- stabsc ~ bout_contr_jaccrepl + I(bout_contr_jaccrepl^2)
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-FD-bout-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrich
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-FD-bout-lin.rds")

form <- stabsc ~ bout_contr_jaccrich + I(bout_contr_jaccrich^2)
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-FD-bout-quad.rds")


## 2nd set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, BAT: site-level means
# Jaccard dissimilarity
form <- stabsc ~ bout_jacc
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-FD-bout-S-lin.rds")

form <- stabsc ~ bout_jacc + I(bout_jacc^2)
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-FD-bout-S-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ bout_jaccrepl
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-FD-bout-S-lin.rds")

form <- stabsc ~ bout_jaccrepl + I(bout_jaccrepl^2)
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-FD-bout-S-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ bout_jaccrich
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-FD-bout-S-lin.rds")

form <- stabsc ~ bout_jaccrich + I(bout_jaccrich^2)
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-FD-bout-S-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrepl
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-FD-bout-S-lin.rds")

form <- stabsc ~ bout_contr_jaccrepl + I(bout_contr_jaccrepl^2)
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-FD-bout-S-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrich
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-FD-bout-S-lin.rds")

form <- stabsc ~ bout_contr_jaccrich + I(bout_contr_jaccrich^2)
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-FD-bout-S-quad.rds")


## 3rd set of models: Year-level (inter-annual) stability ~ Year-level (inter-annual) community change, BAT
# Jaccard dissimilarity
form <- stabsc ~ year_jacc
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-FD-year-lin.rds")

form <- stabsc ~ year_jacc + I(year_jacc^2)
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-FD-year-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ year_jaccrepl
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-FD-year-lin.rds")

form <- stabsc ~ year_jaccrepl + I(year_jaccrepl^2)
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-FD-year-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ year_jaccrich
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-FD-year-lin.rds")

form <- stabsc ~ year_jaccrich + I(year_jaccrich^2)
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-FD-year-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrepl
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-FD-year-lin.rds")

form <- stabsc ~ year_contr_jaccrepl + I(year_contr_jaccrepl^2)
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-FD-year-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrich
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-FD-year-lin.rds")

form <- stabsc ~ year_contr_jaccrich + I(year_contr_jaccrich^2)
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-FD-year-quad.rds")

## 4th set of models: Year-level (inter-annual) stability ~ Year-level (inter-annual) community change, BAT: site-level means
# Jaccard dissimilarity
form <- stabsc ~ year_jacc
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-FD-year-S-lin.rds")

form <- stabsc ~ year_jacc + I(year_jacc^2)
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-FD-year-S-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ year_jaccrepl
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-FD-year-S-lin.rds")

form <- stabsc ~ year_jaccrepl + I(year_jaccrepl^2)
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-FD-year-S-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ year_jaccrich
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-FD-year-S-lin.rds")

form <- stabsc ~ year_jaccrich + I(year_jaccrich^2)
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-FD-year-S-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrepl
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-FD-year-S-lin.rds")

form <- stabsc ~ year_contr_jaccrepl + I(year_contr_jaccrepl^2)
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-FD-year-S-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrich
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-FD-year-S-lin.rds")

form <- stabsc ~ year_contr_jaccrich + I(year_contr_jaccrich^2)
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-FD-year-S-quad.rds")


### 2. FUNCTIONAL COMMUNITY CHANGE vs LATITUDE (total of 2- models)
## 1st set of models: Bout-level (intra-annual) community change vs latitude, BAT
# Jaccard dissimilarity
form <- bout_jacc ~ latsc
m1 <- inla(form, data=fish.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jacc-latsc-FD-bout.rds")

# Jaccard dissimilarity due to replacement
form <- bout_jaccrepl ~ latsc
m1 <- inla(form, data=fish.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccrepl-latsc-FD-bout.rds")

# Jaccard dissimilarity due to richness
form <- bout_jaccrich ~ latsc
m1 <- inla(form, data=fish.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccrich-latsc-FD-bout.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- bout_contr_jaccrepl ~ latsc
m1 <- inla(form, data=fish.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrepl-latsc-FD-bout.rds")

# Contribution of richness to Jaccard dissimilarity
form <- bout_contr_jaccrich ~ latsc
m1 <- inla(form, data=fish.bout.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrich-latsc-FD-bout.rds")


## 2nd set of models: Bout-level (intra-annual) community change vs latitude, BAT: site-level means
# Jaccard dissimilarity
form <- bout_jacc ~ latsc
m1 <- inla(form, data=fish.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jacc-latsc-FD-bout-S.rds")

# Jaccard dissimilarity due to replacement
form <- bout_jaccrepl ~ latsc
m1 <- inla(form, data=fish.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccrepl-latsc-FD-bout-S.rds")

# Jaccard dissimilarity due to richness
form <- bout_jaccrich ~ latsc
m1 <- inla(form, data=fish.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccrich-latsc-FD-bout-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- bout_contr_jaccrepl ~ latsc
m1 <- inla(form, data=fish.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrepl-latsc-FD-bout-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- bout_contr_jaccrich ~ latsc
m1 <- inla(form, data=fish.bout.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrich-latsc-FD-bout-S.rds")


## 3rd set of models: Year-level (inter-annual) community change vs latitude, BAT
# Jaccard dissimilarity
form <- year_jacc ~ latsc
m1 <- inla(form, data=fish.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jacc-latsc-FD-year.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccrepl ~ latsc
m1 <- inla(form, data=fish.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccrepl-latsc-FD-year.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ latsc
m1 <- inla(form, data=fish.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccrich-latsc-FD-year.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- year_contr_jaccrepl ~ latsc
m1 <- inla(form, data=fish.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrepl-latsc-FD-year.rds")

# Contribution of richness to Jaccard dissimilarity
form <- year_contr_jaccrich ~ latsc
m1 <- inla(form, data=fish.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrich-latsc-FD-year.rds")


## 4th set of models: Year-level (inter-annual) community change vs latitude, BAT: site-level means
# Jaccard dissimilarity
form <- year_jacc ~ latsc
m1 <- inla(form, data=fish.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jacc-latsc-FD-year-S.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccrepl ~ latsc
m1 <- inla(form, data=fish.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccrepl-latsc-FD-year-S.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ latsc
m1 <- inla(form, data=fish.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccrich-latsc-FD-year-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- year_contr_jaccrepl ~ latsc
m1 <- inla(form, data=fish.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrepl-latsc-FD-year-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- year_contr_jaccrich ~ latsc
m1 <- inla(form, data=fish.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrich-latsc-FD-year-S.rds")



### 3. YEAR-LEVEL FUNCTIONAL COMMUNITY CHANGE vs BOUT-LEVEL TAXONOMIC COMMUNITY CHANGE (total of 6 models)
fish.year.bat.sel <- fish.year.bat.m %>%
  select(siteID, year_jacc, year_jaccrepl, year_jaccrich)
fish.bout.year.bat <- fish.bout.bat %>%
  left_join(fish.year.bat.sel, by = c("siteID"))

fish.bout.year.bat.m <- fish.bout.bat.m %>%
  left_join(fish.year.bat.sel, by = c("siteID"))

## 1st set of models: BAT
# Jaccard dissimilarity
form <- year_jacc ~ bout_jacc
m1 <- inla(form, data=fish.bout.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearjacc-boutjacc-FD.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccrepl ~ bout_jaccrepl
m1 <- inla(form, data=fish.bout.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearjaccrepl-boutjaccrepl-FD.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ bout_jaccrich
m1 <- inla(form, data=fish.bout.year.bat, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearjaccrich-boutjaccrich-FD.rds")


## 2nd set of models: BAT, Site-level means
# Jaccard dissimilarity
form <- year_jacc ~ bout_jacc
m1 <- inla(form, data=fish.bout.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearjacc-boutjacc-FD-S.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccrepl ~ bout_jaccrepl
m1 <- inla(form, data=fish.bout.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearjaccrepl-boutjaccrepl-FD-S.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrich ~ bout_jaccrich
m1 <- inla(form, data=fish.bout.year.bat.m, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearjaccrich-boutjaccrich-FD-S.rds")











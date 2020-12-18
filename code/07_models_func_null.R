#Load packages
require(rhdf5)
require(here)
require(tidyverse)
select <- dplyr::select
require(INLA)

###################MODELS EXPLAINING ECOSYSTEM STABILITY AS A FUNCTION OF CHANGE IN COMMUNITY FUNCTIONAL COMPOSITION: NULL

################### small mammals
### Read in data on stability and community compositional change
mammal.bout.stab <- readRDS(file="data/mammal_outStability-bout.rds") %>%
  select(siteID, year, bout_cv, bout_stability)
mammal.year.stab <- readRDS(file="data/mammal_outStability-year.rds") %>%
  select(siteID, year_cv, year_stability)
mammal.bout.bat <- readRDS(file="data/mammal_outBAT-FD-SES-bout_null.rds")
mammal.year.bat <- readRDS(file="data/mammal_outBAT-FD-SES-year_null.rds")

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
  summarize(bout_jaccses = mean(bout_jaccses, na.rm=TRUE), bout_jaccreplses = mean(bout_jaccreplses, na.rm=TRUE), bout_jaccrichses = mean(bout_jaccrichses, na.rm=TRUE),
            bout_contr_jaccreplses = mean(bout_contr_jaccreplses, na.rm=TRUE), bout_contr_jaccrichses = mean(bout_contr_jaccrichses, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, year, bout_jaccses, bout_jaccreplses, bout_jaccrichses, bout_contr_jaccreplses, bout_contr_jaccrichses, bout_cv, bout_stability) %>%
  ungroup

mammal.year.bat.m <- mammal.year.bat %>%
  group_by(domainID, siteID, lat, lon, year_cv, year_stability) %>%
  summarize(year_jaccses = mean(year_jaccses, na.rm=TRUE), year_jaccreplses = mean(year_jaccreplses, na.rm=TRUE), year_jaccrichses = mean(year_jaccrichses, na.rm=TRUE),
            year_contr_jaccreplses = mean(year_contr_jaccreplses, na.rm=TRUE), year_contr_jaccrichses = mean(year_contr_jaccrichses, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, year_jaccses, year_jaccreplses, year_jaccrichses, year_contr_jaccreplses, year_contr_jaccrichses, year_cv, year_stability) %>%
  ungroup

# scale stability, latitude, and longitude, remove NAs, replace 0 and 1 values for BAT metrics and rank change with 0.0001 and 0.999 so that models with beta distribution can be fit
mammal.bout.bat <- mammal.bout.bat %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(bout_stability))

mammal.year.bat <- mammal.year.bat %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(year_stability))

mammal.bout.bat.m <- mammal.bout.bat.m %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(bout_stability))

mammal.year.bat.m <- mammal.year.bat.m %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(year_stability)) 



##### MODELS
### 1. STABILITY vs FUNCTIONAL COMMUNITY CHANGE (total of 20 models)
## 1st set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, BAT
# Jaccard dissimilarity
form <- stabsc ~ bout_jaccses
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-FD-SES-bout-lin.rds")

form <- stabsc ~ bout_jaccses + I(bout_jaccses^2)
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-FD-SES-bout-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ bout_jaccreplses
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-FD-SES-bout-lin.rds")

form <- stabsc ~ bout_jaccreplses + I(bout_jaccreplses^2)
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-FD-SES-bout-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ bout_jaccrichses
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-FD-SES-bout-lin.rds")

form <- stabsc ~ bout_jaccrichses + I(bout_jaccrichses^2)
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-FD-SES-bout-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccreplses
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-FD-SES-bout-lin.rds")

form <- stabsc ~ bout_contr_jaccreplses + I(bout_contr_jaccreplses)
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-FD-SES-bout-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrichses
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-FD-SES-bout-lin.rds")

form <- stabsc ~ bout_contr_jaccrichses + I(bout_contr_jaccrichses)
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-FD-SES-bout-quad.rds")

## 2nd set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, BAT: site-level means
# Jaccard dissimilarity
form <- stabsc ~ bout_jaccses
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-FD-SES-bout-S-lin.rds")

form <- stabsc ~ bout_jaccses + I(bout_jaccses)
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-FD-SES-bout-S-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ bout_jaccreplses
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-FD-SES-bout-S-lin.rds")

form <- stabsc ~ bout_jaccreplses + I(bout_jaccreplses)
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-FD-SES-bout-S-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ bout_jaccrichses
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-FD-SES-bout-S-lin.rds")

form <- stabsc ~ bout_jaccrichses + I(bout_jaccrichses)
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-FD-SES-bout-S-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccreplses
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-FD-SES-bout-S-lin.rds")

form <- stabsc ~ bout_contr_jaccreplses + I(bout_contr_jaccreplses)
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-FD-SES-bout-S-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrichses
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-FD-SES-bout-S-lin.rds")

form <- stabsc ~ bout_contr_jaccrichses + I(bout_contr_jaccrichses)
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-FD-SES-bout-S-quad.rds")

## 3rd set of models: Year-level (inter-annual) stability ~ Year-level (inter-annual) community change, BAT
# Jaccard dissimilarity
form <- stabsc ~ year_jaccses
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-FD-SES-year-lin.rds")

form <- stabsc ~ year_jaccses + I(year_jaccses)
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-FD-SES-year-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ year_jaccreplses
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-FD-SES-year-lin.rds")

form <- stabsc ~ year_jaccreplses + I(year_jaccreplses)
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-FD-SES-year-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ year_jaccrichses
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-FD-SES-year-lin.rds")

form <- stabsc ~ year_jaccrichses + I(year_jaccrichses)
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-FD-SES-year-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccreplses
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-FD-SES-year-lin.rds")

form <- stabsc ~ year_contr_jaccreplses + I(year_contr_jaccreplses)
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-FD-SES-year-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrichses
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-FD-SES-year-lin.rds")

form <- stabsc ~ year_contr_jaccrichses + I(year_contr_jaccrichses)
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-FD-SES-year-quad.rds")

## 4th set of models: Year-level (inter-annual) stability ~ Year-level (inter-annual) community change, BAT: site-level means
# Jaccard dissimilarity
form <- stabsc ~ year_jaccses
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-FD-SES-year-S-lin.rds")

form <- stabsc ~ year_jaccses + I(year_jaccses)
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jacc-FD-SES-year-S-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ year_jaccreplses
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-FD-SES-year-S-lin.rds")

form <- stabsc ~ year_jaccreplses + I(year_jaccreplses)
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrepl-FD-SES-year-S-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ year_jaccrichses
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-FD-SES-year-S-lin.rds")

form <- stabsc ~ year_jaccrichses + I(year_jaccrichses)
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-jaccrich-FD-SES-year-S-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccreplses
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-FD-SES-year-S-lin.rds")

form <- stabsc ~ year_contr_jaccreplses + I(year_contr_jaccreplses)
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrepl-FD-SES-year-S-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrichses
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-FD-SES-year-S-lin.rds")

form <- stabsc ~ year_contr_jaccrichses + I(year_contr_jaccrichses)
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_stabsc-contrjaccrich-FD-SES-year-S-quad.rds")


### 2. FUNCTIONAL COMMUNITY CHANGE vs LATITUDE (total of 2- models)
## 1st set of models: Bout-level (intra-annual) community change vs latitude, BAT
# Jaccard dissimilarity
form <- bout_jaccses ~ latsc
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jacc-latsc-FD-SES-bout.rds")

# Jaccard dissimilarity due to replacement
form <- bout_jaccreplses ~ latsc
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccreplses-latsc-FD-SES-bout.rds")

# Jaccard dissimilarity due to richness
form <- bout_jaccrichses ~ latsc
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccrichses-latsc-FD-SES-bout.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- bout_contr_jaccreplses ~ latsc
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrepl-latsc-FD-SES-bout.rds")

# Contribution of richness to Jaccard dissimilarity
form <- bout_contr_jaccrichses ~ latsc
m1 <- inla(form, data=mammal.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrich-latsc-FD-SES-bout.rds")


## 2nd set of models: Bout-level (intra-annual) community change vs latitude, BAT: site-level means
# Jaccard dissimilarity
form <- bout_jaccses ~ latsc
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jacc-latsc-FD-SES-bout-S.rds")

# Jaccard dissimilarity due to replacement
form <- bout_jaccreplses ~ latsc
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccreplses-latsc-FD-SES-bout-S.rds")

# Jaccard dissimilarity due to richness
form <- bout_jaccrichses ~ latsc
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccrichses-latsc-FD-SES-bout-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- bout_contr_jaccreplses ~ latsc
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrepl-latsc-FD-SES-bout-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- bout_contr_jaccrichses ~ latsc
m1 <- inla(form, data=mammal.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrich-latsc-FD-SES-bout-S.rds")


## 3rd set of models: Year-level (inter-annual) community change vs latitude, BAT
# Jaccard dissimilarity
form <- year_jaccses ~ latsc
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jacc-latsc-FD-SES-year.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccreplses ~ latsc
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccreplses-latsc-FD-SES-year.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrichses ~ latsc
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccrichses-latsc-FD-SES-year.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- year_contr_jaccreplses ~ latsc
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrepl-latsc-FD-SES-year.rds")

# Contribution of richness to Jaccard dissimilarity
form <- year_contr_jaccrichses ~ latsc
m1 <- inla(form, data=mammal.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrich-latsc-FD-SES-year.rds")


## 4th set of models: Year-level (inter-annual) community change vs latitude, BAT: site-level means
# Jaccard dissimilarity
form <- year_jaccses ~ latsc
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jacc-latsc-FD-SES-year-S.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccreplses ~ latsc
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccreplses-latsc-FD-SES-year-S.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrichses ~ latsc
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_jaccrichses-latsc-FD-SES-year-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- year_contr_jaccreplses ~ latsc
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrepl-latsc-FD-SES-year-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- year_contr_jaccrichses ~ latsc
m1 <- inla(form, data=mammal.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_contrjaccrich-latsc-FD-SES-year-S.rds")



### 3. YEAR-LEVEL FUNCTIONAL COMMUNITY CHANGE vs BOUT-LEVEL TAXONOMIC COMMUNITY CHANGE (total of 6 models)
mammal.year.bat.sel <- mammal.year.bat.m %>%
  select(siteID, year_jaccses, year_jaccreplses, year_jaccrichses)
mammal.bout.year.bat <- mammal.bout.bat %>%
  left_join(mammal.year.bat.sel, by = c("siteID"))

mammal.bout.year.bat.m <- mammal.bout.bat.m %>%
  left_join(mammal.year.bat.sel, by = c("siteID"))

## 1st set of models: BAT
# Jaccard dissimilarity
form <- year_jaccses ~ bout_jaccses
m1 <- inla(form, data=mammal.bout.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearjacc-boutjacc-FD-SES.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccreplses ~ bout_jaccreplses
m1 <- inla(form, data=mammal.bout.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearjaccrepl-boutjaccrepl-FD-SES.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrichses ~ bout_jaccrichses
m1 <- inla(form, data=mammal.bout.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearjaccrich-boutjaccrich-FD-SES.rds")


## 2nd set of models: BAT, Site-level means
# Jaccard dissimilarity
form <- year_jaccses ~ bout_jaccses
m1 <- inla(form, data=mammal.bout.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearjacc-boutjacc-FD-SES-S.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccreplses ~ bout_jaccreplses
m1 <- inla(form, data=mammal.bout.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearjaccrepl-boutjaccrepl-FD-SES-S.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrichses ~ bout_jaccrichses
m1 <- inla(form, data=mammal.bout.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/mammal_yearjaccrich-boutjaccrich-FD-SES-S.rds")






################### fish
### Read in data on stability and community compositional change
fish.bout.stab <- readRDS(file="data/fish_outStability-bout.rds") %>%
  select(siteID, year, bout_cv, bout_stability)
fish.year.stab <- readRDS(file="data/fish_outStability-year.rds") %>%
  select(siteID, year_cv, year_stability)
fish.bout.bat <- readRDS(file="data/fish_outBAT-FD-SES-bout_null.rds")
fish.year.bat <- readRDS(file="data/fish_outBAT-FD-SES-year_null.rds")

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
  summarize(bout_jaccses = mean(bout_jaccses, na.rm=TRUE), bout_jaccreplses = mean(bout_jaccreplses, na.rm=TRUE), bout_jaccrichses = mean(bout_jaccrichses, na.rm=TRUE),
            bout_contr_jaccreplses = mean(bout_contr_jaccreplses, na.rm=TRUE), bout_contr_jaccrichses = mean(bout_contr_jaccrichses, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, year, bout_jaccses, bout_jaccreplses, bout_jaccrichses, bout_contr_jaccreplses, bout_contr_jaccrichses, bout_cv, bout_stability) %>%
  ungroup

fish.year.bat.m <- fish.year.bat %>%
  group_by(domainID, siteID, lat, lon, year_cv, year_stability) %>%
  summarize(year_jaccses = mean(year_jaccses, na.rm=TRUE), year_jaccreplses = mean(year_jaccreplses, na.rm=TRUE), year_jaccrichses = mean(year_jaccrichses, na.rm=TRUE),
            year_contr_jaccreplses = mean(year_contr_jaccreplses, na.rm=TRUE), year_contr_jaccrichses = mean(year_contr_jaccrichses, na.rm=TRUE)) %>%
  select(domainID, siteID, lat, lon, year_jaccses, year_jaccreplses, year_jaccrichses, year_contr_jaccreplses, year_contr_jaccrichses, year_cv, year_stability) %>%
  ungroup

# scale stability, latitude, and longitude, remove NAs, replace 0 and 1 values for BAT metrics and rank change with 0.0001 and 0.999 so that models with beta distribution can be fit
fish.bout.bat <- fish.bout.bat %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(bout_stability))

fish.year.bat <- fish.year.bat %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(year_stability))

fish.bout.bat.m <- fish.bout.bat.m %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(bout_stability))

fish.year.bat.m <- fish.year.bat.m %>% 
  filter_all(all_vars(!is.infinite(.))) %>% #remove ecosystem stability = Inf
  mutate(latsc = scale(lat), lonsc = scale(lon), stabsc = scale(year_stability)) 



##### MODELS
### 1. STABILITY vs FUNCTIONAL COMMUNITY CHANGE (total of 20 models)
## 1st set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, BAT
# Jaccard dissimilarity
form <- stabsc ~ bout_jaccses
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-FD-SES-bout-lin.rds")

form <- stabsc ~ bout_jaccses + I(bout_jaccses^2)
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-FD-SES-bout-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ bout_jaccreplses
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-FD-SES-bout-lin.rds")

form <- stabsc ~ bout_jaccreplses + I(bout_jaccreplses^2)
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-FD-SES-bout-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ bout_jaccrichses
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-FD-SES-bout-lin.rds")

form <- stabsc ~ bout_jaccrichses + I(bout_jaccrichses^2)
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-FD-SES-bout-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccreplses
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-FD-SES-bout-lin.rds")

form <- stabsc ~ bout_contr_jaccreplses + I(bout_contr_jaccreplses)
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-FD-SES-bout-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrichses
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-FD-SES-bout-lin.rds")

form <- stabsc ~ bout_contr_jaccrichses + I(bout_contr_jaccrichses)
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-FD-SES-bout-quad.rds")

## 2nd set of models: Bout-level (intra-annual) stability ~ Bout-level (intra-annual) community change, BAT: site-level means
# Jaccard dissimilarity
form <- stabsc ~ bout_jaccses
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-FD-SES-bout-S-lin.rds")

form <- stabsc ~ bout_jaccses + I(bout_jaccses)
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-FD-SES-bout-S-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ bout_jaccreplses
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-FD-SES-bout-S-lin.rds")

form <- stabsc ~ bout_jaccreplses + I(bout_jaccreplses)
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-FD-SES-bout-S-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ bout_jaccrichses
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-FD-SES-bout-S-lin.rds")

form <- stabsc ~ bout_jaccrichses + I(bout_jaccrichses)
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-FD-SES-bout-S-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccreplses
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-FD-SES-bout-S-lin.rds")

form <- stabsc ~ bout_contr_jaccreplses + I(bout_contr_jaccreplses)
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-FD-SES-bout-S-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ bout_contr_jaccrichses
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-FD-SES-bout-S-lin.rds")

form <- stabsc ~ bout_contr_jaccrichses + I(bout_contr_jaccrichses)
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-FD-SES-bout-S-quad.rds")

## 3rd set of models: Year-level (inter-annual) stability ~ Year-level (inter-annual) community change, BAT
# Jaccard dissimilarity
form <- stabsc ~ year_jaccses
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-FD-SES-year-lin.rds")

form <- stabsc ~ year_jaccses + I(year_jaccses)
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-FD-SES-year-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ year_jaccreplses
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-FD-SES-year-lin.rds")

form <- stabsc ~ year_jaccreplses + I(year_jaccreplses)
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-FD-SES-year-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ year_jaccrichses
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-FD-SES-year-lin.rds")

form <- stabsc ~ year_jaccrichses + I(year_jaccrichses)
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-FD-SES-year-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccreplses
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-FD-SES-year-lin.rds")

form <- stabsc ~ year_contr_jaccreplses + I(year_contr_jaccreplses)
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-FD-SES-year-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrichses
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-FD-SES-year-lin.rds")

form <- stabsc ~ year_contr_jaccrichses + I(year_contr_jaccrichses)
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-FD-SES-year-quad.rds")

## 4th set of models: Year-level (inter-annual) stability ~ Year-level (inter-annual) community change, BAT: site-level means
# Jaccard dissimilarity
form <- stabsc ~ year_jaccses
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-FD-SES-year-S-lin.rds")

form <- stabsc ~ year_jaccses + I(year_jaccses)
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jacc-FD-SES-year-S-quad.rds")

# Jaccard dissimilarity due to replacement
form <- stabsc ~ year_jaccreplses
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-FD-SES-year-S-lin.rds")

form <- stabsc ~ year_jaccreplses + I(year_jaccreplses)
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrepl-FD-SES-year-S-quad.rds")

# Jaccard dissimilarity due to richness
form <- stabsc ~ year_jaccrichses
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-FD-SES-year-S-lin.rds")

form <- stabsc ~ year_jaccrichses + I(year_jaccrichses)
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-jaccrich-FD-SES-year-S-quad.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccreplses
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-FD-SES-year-S-lin.rds")

form <- stabsc ~ year_contr_jaccreplses + I(year_contr_jaccreplses)
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrepl-FD-SES-year-S-quad.rds")

# Contribution of richness to Jaccard dissimilarity
form <- stabsc ~ year_contr_jaccrichses
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-FD-SES-year-S-lin.rds")

form <- stabsc ~ year_contr_jaccrichses + I(year_contr_jaccrichses)
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_stabsc-contrjaccrich-FD-SES-year-S-quad.rds")


### 2. FUNCTIONAL COMMUNITY CHANGE vs LATITUDE (total of 2- models)
## 1st set of models: Bout-level (intra-annual) community change vs latitude, BAT
# Jaccard dissimilarity
form <- bout_jaccses ~ latsc
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jacc-latsc-FD-SES-bout.rds")

# Jaccard dissimilarity due to replacement
form <- bout_jaccreplses ~ latsc
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccreplses-latsc-FD-SES-bout.rds")

# Jaccard dissimilarity due to richness
form <- bout_jaccrichses ~ latsc
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccrichses-latsc-FD-SES-bout.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- bout_contr_jaccreplses ~ latsc
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrepl-latsc-FD-SES-bout.rds")

# Contribution of richness to Jaccard dissimilarity
form <- bout_contr_jaccrichses ~ latsc
m1 <- inla(form, data=fish.bout.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrich-latsc-FD-SES-bout.rds")


## 2nd set of models: Bout-level (intra-annual) community change vs latitude, BAT: site-level means
# Jaccard dissimilarity
form <- bout_jaccses ~ latsc
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jacc-latsc-FD-SES-bout-S.rds")

# Jaccard dissimilarity due to replacement
form <- bout_jaccreplses ~ latsc
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccreplses-latsc-FD-SES-bout-S.rds")

# Jaccard dissimilarity due to richness
form <- bout_jaccrichses ~ latsc
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccrichses-latsc-FD-SES-bout-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- bout_contr_jaccreplses ~ latsc
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrepl-latsc-FD-SES-bout-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- bout_contr_jaccrichses ~ latsc
m1 <- inla(form, data=fish.bout.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrich-latsc-FD-SES-bout-S.rds")


## 3rd set of models: Year-level (inter-annual) community change vs latitude, BAT
# Jaccard dissimilarity
form <- year_jaccses ~ latsc
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jacc-latsc-FD-SES-year.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccreplses ~ latsc
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccreplses-latsc-FD-SES-year.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrichses ~ latsc
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccrichses-latsc-FD-SES-year.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- year_contr_jaccreplses ~ latsc
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrepl-latsc-FD-SES-year.rds")

# Contribution of richness to Jaccard dissimilarity
form <- year_contr_jaccrichses ~ latsc
m1 <- inla(form, data=fish.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrich-latsc-FD-SES-year.rds")


## 4th set of models: Year-level (inter-annual) community change vs latitude, BAT: site-level means
# Jaccard dissimilarity
form <- year_jaccses ~ latsc
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jacc-latsc-FD-SES-year-S.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccreplses ~ latsc
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccreplses-latsc-FD-SES-year-S.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrichses ~ latsc
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_jaccrichses-latsc-FD-SES-year-S.rds")

# Contribution of replacament to Jaccard dissimilarity
form <- year_contr_jaccreplses ~ latsc
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrepl-latsc-FD-SES-year-S.rds")

# Contribution of richness to Jaccard dissimilarity
form <- year_contr_jaccrichses ~ latsc
m1 <- inla(form, data=fish.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_contrjaccrich-latsc-FD-SES-year-S.rds")



### 3.YEAR-LEVEL FUNCTIONAL COMMUNITY CHANGE vs BOUT-LEVEL TAXONOMIC COMMUNITY CHANGE (total of 6 models)
fish.year.bat.sel <- fish.year.bat.m %>%
  select(siteID, year_jaccses, year_jaccreplses, year_jaccrichses)
fish.bout.year.bat <- fish.bout.bat %>%
  left_join(fish.year.bat.sel, by = c("siteID"))

fish.bout.year.bat.m <- fish.bout.bat.m %>%
  left_join(fish.year.bat.sel, by = c("siteID"))

## 1st set of models: BAT
# Jaccard dissimilarity
form <- year_jaccses ~ bout_jaccses
m1 <- inla(form, data=fish.bout.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearjacc-boutjacc-FD-SES.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccreplses ~ bout_jaccreplses
m1 <- inla(form, data=fish.bout.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearjaccrepl-boutjaccrepl-FD-SES.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrichses ~ bout_jaccrichses
m1 <- inla(form, data=fish.bout.year.bat, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearjaccrich-boutjaccrich-FD-SES.rds")


## 2nd set of models: BAT, Site-level means
# Jaccard dissimilarity
form <- year_jaccses ~ bout_jaccses
m1 <- inla(form, data=fish.bout.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearjacc-boutjacc-FD-SES-S.rds")

# Jaccard dissimilarity due to replacement
form <- year_jaccreplses ~ bout_jaccreplses
m1 <- inla(form, data=fish.bout.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearjaccrepl-boutjaccrepl-FD-SES-S.rds")

# Jaccard dissimilarity due to richness
form <- year_jaccrichses ~ bout_jaccrichses
m1 <- inla(form, data=fish.bout.year.bat.m, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE), 
           control.predictor=list(compute=TRUE, link=1))
saveRDS(m1, file="Models_out/fish_yearjaccrich-boutjaccrich-FD-SES-S.rds")










#Load packages
require(tidyverse)
require(dplyr)
'%>%' <- dplyr::`%>%`


###################REORGANIZE DATA: GET INTRA- AND INTER- BIODIVERSITY

################### small mammals
data_small_mammal <- readRDS(file = "data/mammal_processed.rds")

#pivot wide per day
data_small_mammal_group <- data_small_mammal %>%
  select(domainID, siteID, lat, lon, year, month, day, bout, scientificName) %>%
  group_by(domainID, siteID, lat, lon, year, month, day, bout, scientificName) %>%
  tally(name = "count") %>% 
  ungroup()

data_small_mammal_day <- data_small_mammal_group %>%
  pivot_wider(id_cols = c(domainID, siteID, lat, lon, year, month, day, bout), 
              names_from = scientificName, values_from = count,
              values_fill = list(count = 0)) %>%
  arrange(domainID, siteID, lat, lon, year, month, day, bout) 

saveRDS(data_small_mammal_day, file="data/mammal_abund-day.rds")

#pivot wide per bout
data_small_mammal_group <- data_small_mammal %>%
  select(domainID, siteID, lat, lon, year, month, bout, scientificName) %>%
  group_by(domainID, siteID, lat, lon, year, month, bout, scientificName) %>%
  tally(name = "count") %>% 
  ungroup()

data_small_mammal_bout <- data_small_mammal_group %>%
  pivot_wider(id_cols = c(domainID, siteID, lat, lon, year, month, bout), 
              names_from = scientificName, values_from = count,
              values_fill = list(count = 0)) %>%
  arrange(domainID, siteID, lat, lon, year, month, bout) 

saveRDS(data_small_mammal_bout, file="data/mammal_abund-bout.rds")

#pivot wide per year
data_small_mammal_group <- data_small_mammal %>%
  select(domainID, siteID, lat, lon, year, scientificName) %>%
  group_by(domainID, siteID, lat, lon, year, scientificName) %>%
  tally(name = "count") %>% 
  ungroup()

data_small_mammal_year <- data_small_mammal_group %>%
  pivot_wider(id_cols = c(domainID, siteID, lat, lon, year), 
              names_from = scientificName, values_from = count,
              values_fill = list(count = 0)) %>%
  arrange(domainID, siteID, lat, lon, year) 

saveRDS(data_small_mammal_year, file="data/mammal_abund-year.rds")



################### ground beetles
data_beetle <- readRDS(file = "data/beetle_processed.rds")

# Bout-level
# Calculate the total trapping days for a bout (abundance must be standardized by the duration of trapping)
bout_trap_days <- data_beetle %>%
  select(siteID, plotID, trapID, boutID, trappingDays) %>%
  distinct() %>%
  group_by(siteID, boutID) %>%
  mutate(totalTrapDays = sum(trappingDays)) %>%
  select(-c(trapID, plotID, trappingDays)) %>%
  distinct()

data_beetle_group <- data_beetle %>% 
  filter(!is.na(count)) %>%
  select(domainID, siteID, boutID, year, lat, lon, scientificName, count) %>%
  group_by(domainID, siteID, year, boutID, scientificName) %>%
  mutate(abundance = sum(count)) %>%
  select(-count) %>%
  distinct() %>%
  ungroup() %>%
  left_join(bout_trap_days, by = c("siteID", "boutID"))

data_beetle_group_trap <- data_beetle_group %>%
  mutate(cpue = abundance/totalTrapDays) %>%
  select(-c(abundance, totalTrapDays)) %>%
  mutate(bout = boutID)

#pivot wide per bout
data_beetle_bout <- data_beetle_group_trap %>%
  select(domainID, siteID, lat, lon, year, bout, scientificName, cpue) %>%
  pivot_wider(id_cols = c(domainID, siteID, lat, lon, year, bout), 
              names_from = scientificName, values_from = cpue,
              values_fill = list(cpue = 0)) %>%
  arrange(domainID, siteID, lat, lon, year, bout) 

saveRDS(data_beetle_bout, file="data/beetle_abund-bout.rds")


# Year-level
yearly_trap_days <- bout_trap_days %>%
  mutate(year = substr(boutID, 6,9)) %>%
  group_by(siteID, year) %>%
  mutate(totalTrapDays = sum(totalTrapDays)) %>%
  select(-boutID) %>%
  distinct()

data_beetle_group_trapyear <- data_beetle_group %>%
  mutate(year = substr(boutID, 6,9)) %>%
  select(-c(totalTrapDays, boutID)) %>%
  group_by(siteID, year, scientificName) %>%
  mutate(abundance = sum(abundance)) %>%
  distinct() %>%
  left_join(yearly_trap_days, by = c("siteID", "year")) %>%
  mutate(cpue = abundance/totalTrapDays) %>%
  select(-c(abundance, totalTrapDays)) %>%
  ungroup()

#pivot wide per year
data_beetle_year <- data_beetle_group_trapyear %>%
  select(domainID, siteID, lat, lon, year, scientificName, cpue) %>%
  pivot_wider(id_cols = c(domainID, siteID, lat, lon, year), 
              names_from = scientificName, values_from = cpue,
              values_fill = list(cpue = 0)) %>%
  arrange(domainID, siteID, lat, lon, year) 

saveRDS(data_beetle_year, file="data/beetle_abund-year.rds")



################### fish
data_fish <- readRDS(file = "data/fish_processed.rds")

#for beet sum abundance over the bout/year, and then divide by total trap days in that bout or year
#for fish...?

#pivot wide per bout
data_fish_group <- data_fish %>% 
  select(domainID, siteID, month, year, bout, lat, lon, scientificName, number_of_fish, mean_efishtime) %>%
  group_by(domainID, siteID, month, year, bout, lat, lon, scientificName) %>%
  mutate(mean_efishtime = sum(mean_efishtime, na.rm=TRUE)) %>%
  mutate(number_of_fish = sum(number_of_fish)) %>%
  mutate(cpue = number_of_fish/mean_efishtime * 3600) %>%
  mutate(cpue = ifelse(cpue == "Inf", 0, cpue)) %>%
  distinct() %>%
  ungroup() 

data_fish_bout <- data_fish_group %>%
  pivot_wider(id_cols = c(domainID, siteID, lat, lon, month, year, bout), 
              names_from = scientificName, values_from = cpue,
              values_fill = list(cpue = 0)) %>%
  arrange(domainID, siteID, lat, lon, year, month, bout) 

saveRDS(data_fish_bout, file="data/fish_abund-bout.rds")

#pivot wide per year
data_fish_group <- data_fish %>% 
  select(domainID, siteID, year, lat, lon, scientificName, number_of_fish, mean_efishtime) %>%
  group_by(domainID, siteID, year, lat, lon, scientificName) %>%
  mutate(mean_efishtime = sum(mean_efishtime, na.rm=TRUE)) %>%
  mutate(number_of_fish = sum(number_of_fish)) %>%
  mutate(cpue = number_of_fish/mean_efishtime * 3600) %>%
  mutate(cpue = ifelse(cpue == "Inf", 0, cpue)) %>%
  distinct() %>%
  ungroup() 

data_fish_year <- data_fish_group %>%
  pivot_wider(id_cols = c(domainID, siteID, lat, lon, year), 
              names_from = scientificName, values_from = cpue,
              values_fill = list(cpue = 0)) %>%
  arrange(domainID, siteID, lat, lon, year) 

saveRDS(data_fish_year, file="data/fish_abund-year.rds")



################### aquatic macroinvertebrates
data_macroinv <- readRDS(file = "data/macroinv_processed.rds")

#pivot wide per bout
data_macroinv_group <- data_macroinv %>% 
  select(domainID, siteID, month, year, bout, lat, lon, scientificName, estimatedTotalCount, benthicArea) %>%
  group_by(domainID, siteID, month, year, bout, lat, lon, scientificName) %>%
  mutate(estimatedTotalCount = sum(estimatedTotalCount, na.rm=TRUE)) %>%
  mutate(benthicArea = sum(benthicArea)) %>%
  mutate(den = estimatedTotalCount/benthicArea) %>%
  distinct() %>%
  ungroup() 

data_macroinv_bout <- data_macroinv_group %>%
  pivot_wider(id_cols = c(domainID, siteID, lat, lon, month, year, bout), 
              names_from = scientificName, values_from = den,
              values_fill = list(den = 0)) %>%
  arrange(domainID, siteID, lat, lon, year, month, bout) 

saveRDS(data_macroinv_bout, file="data/macroinv_abund-bout.rds")

#pivot wide per year
data_macroinv_group <- data_macroinv %>% 
  select(domainID, siteID, year, lat, lon, scientificName, estimatedTotalCount, benthicArea) %>%
  group_by(domainID, siteID, year, lat, lon, scientificName) %>%
  mutate(estimatedTotalCount = sum(estimatedTotalCount, na.rm=TRUE)) %>%
  mutate(benthicArea = sum(benthicArea)) %>%
  mutate(den = estimatedTotalCount/benthicArea) %>%
  distinct() %>%
  ungroup() 

data_macroinv_year <- data_macroinv_group %>%
  pivot_wider(id_cols = c(domainID, siteID, lat, lon, year), 
              names_from = scientificName, values_from = den,
              values_fill = list(den = 0)) %>%
  arrange(domainID, siteID, lat, lon, year) 
data_macroinv_year
saveRDS(data_macroinv_year, file="data/macroinv_abund-year.rds")



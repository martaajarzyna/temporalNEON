
# set options
options(stringsAsFactors = FALSE)

## getting all required packages
library(installr)
library(neonUtilities)
library(tidyverse)
library(googledrive)
library(devtools)
library(ecocommDP)
library(lubridate)
library(iNEXT)
library(xlsx)
library(geoNEON)
library(pkgbuild)
library(rhdf5)
library(httr)
library(jsonlite)
library(dplyr, quietly=T)
library(downloader)
library(magrittr)
library(sjmisc)
library(dbplyr)
library(DescTools)
library(data.table)
library(readr)
require(raster)

options(stringsAsFactors=F)

# fish dpid
my_dpid_fish <- 'DP1.20107.001'

# get taxon table from API, may take a few minutes to load
# Fish electrofishing, gill netting, and fyke netting counts 
# http://data.neonscience.org/data-product-view?dpCode=DP1.20107.001
full_taxon_fish <- neonUtilities::getTaxonTable(taxonType = 'FISH', recordReturnLimit = NA, stream = "true") 

# -- make ordered taxon_rank list for a reference (subspecies is smallest rank, kingdom is largest)
# a much simple table with useful levels of taxonomic resolution; 
# this might not be needed if taxon rank is extracted  from scientific names using stringr functions 
taxon_rank_fish <- c('superclass', 'class', 'subclass', 'infraclass', 'superorder',
                     'order', 'suborder', 'infraorder', 'section', 'subsection',
                     'superfamily', 'family', 'subfamily', 'tribe', 'subtribe', 'genus',
                     'subgenus','speciesGroup','species','subspecies') %>% rev() # get the reversed version


# get data FISH via api -- will take a while --  package = "basic" is also possible 
all_fish <- neonUtilities::loadByProduct(dpID = my_dpid_fish,  site = "all", startdate = NA, enddate = NA,
                                         package = "expanded", avg = "all", check.size = FALSE)

save(all_fish, file = "all_fish.RData" )
fsh_perFish = all_fish$fsh_perFish
save(fsh_perFish, file = "perFsh.RData")
load(file = "all_fish.RData")

# this joins reach-level data, you can just join by any two of these, and the results are the same, dropping the unique ID
fsh_dat1 <- dplyr::left_join(x = all_fish$fsh_perPass[-1], y = all_fish$fsh_fieldData[-1], by = c("reachID")) %>% 
  dplyr::filter(is.na(samplingImpractical) | samplingImpractical == "") #remove records where fish couldn't be collected, both na's and blank data is kept

# get rid of duplicate col names and .x suffix
fsh_dat1 <- fsh_dat1[,!grepl('\\.y',names(fsh_dat1))]
names(fsh_dat1) <- gsub('\\.x','',names(fsh_dat1))

# add individual fish counts; perFish data = individual level per each specimen captured
# this joins reach-level field data with individual fish measures
fsh_dat_indiv <- dplyr::left_join(all_fish$fsh_perFish, fsh_dat1, by = c("eventID")) 

# get rid of dupe col names and .x suffix
fsh_dat_indiv <- fsh_dat_indiv[,!grepl('\\.y',names(fsh_dat_indiv))]
names(fsh_dat_indiv) <- gsub('\\.x','',names(fsh_dat_indiv))

# fill in missing reachID with event ID
fsh_dat_indiv$reachID <- ifelse(is.na(fsh_dat_indiv$reachID), 
                                substr(fsh_dat_indiv$eventID, 1, 16), fsh_dat_indiv$reachID) 

# add bulk fish counts; bulk count data = The number of fish counted during bulk processing in each pass
# this will join all the reach-level field data with the species data with bulk counts
# fsh_dat_bulk <- dplyr::left_join(all_fish$fsh_bulkCount, fsh_dat1, by = c("eventID", "domainID", "siteID", "namedLocation"))
fsh_dat_bulk <- dplyr::left_join(all_fish$fsh_bulkCount, fsh_dat1, by = c("eventID"))

# get rid of dupe col names and .x suffix
fsh_dat_bulk <- fsh_dat_bulk[,!grepl('\\.y',names(fsh_dat_bulk))]
names(fsh_dat_bulk) <- gsub('\\.x','',names(fsh_dat_bulk))

#fill in missing reachID
fsh_dat_bulk$reachID <- ifelse(is.na(fsh_dat_bulk$reachID), 
                               substr(fsh_dat_bulk$eventID, 1, 16), fsh_dat_bulk$reachID) 

# combine indiv and bulk counts
fsh_dat <-dplyr::bind_rows(fsh_dat_indiv, fsh_dat_bulk)

# add count = 1 for indiv data
# before row_bind, the indiv dataset did not have a col for count
# after bind, the bulk count col has "NAs" need to add "1", since indiv col has individual fish per row
fsh_dat$count <- ifelse(is.na(fsh_dat$bulkFishCount), 1, fsh_dat$bulkFishCount)

# need to convert POSIXct format into as.character and then back to date-time format
# then, fill in missing site ID info and missing startDate into
fsh_dat$startDate <- dplyr::if_else(is.na(fsh_dat$startDate),
                                    lubridate::as_datetime(substr(as.character(fsh_dat$passStartTime), 1, 10)), fsh_dat$startDate) # 1-10: number of characters on date
fsh_dat$siteID <- dplyr::if_else(is.na(fsh_dat$siteID), 
                                 substr(fsh_dat$eventID, 1, 4), fsh_dat$siteID) # four characters on site

# get species and finer res, may not need this if wildcards are used with scientific names col
fsh_dat$taxonRank_ordered <- factor(fsh_dat$taxonRank, levels = taxon_rank_fish, ordered = TRUE) 

#############################################################################################

# in case the bulk count data set does not have a column on taxonomic rank, as such, need to add it here.
# in scientificName column-- species identified below species level (genus, family, order, phylum)-- appear as sp. or spp.
# both sp. and spp. identifications should be marked low-res identifications, aka above species level in ranking  
# and exclude from this analyses
fsh_dat2 <- fsh_dat %>% 
  dplyr::mutate(taxonRank = dplyr::case_when(stringr::str_detect(string = scientificName, pattern = " spp.$") ~ "not_sp_level1", 
                                             stringr::str_detect(string = scientificName, pattern = " sp.$") ~ "not_sp_level2", TRUE ~ "species")) 

save(fsh_dat2, file = "fsh_dat2.RData")

# get all records that have rank <= species; if bulkCount data set has taxonRank, skip this
fsh_dat_fine <- fsh_dat2 %>% dplyr::filter(taxonRank == "species")

# if bulkCount data set has taxonRank, use this this might be added now. 
fsh_dat_fine <- fsh_dat %>%
  dplyr::filter(taxonRank_ordered <= 'species')

# grouping vars for aggregating density measurements 
my_grouping_vars <- c('domainID','siteID','aquaticSiteType','namedLocation',
                      'startDate', 'endDate', 'reachID','eventID','samplerType', 'aquaticSiteType', 'netSetTime', 'netEndTime',
                      'netDeploymentTime', 'netLength', 'netDepth', 'fixedRandomReach','measuredReachLength','efTime', 
                      'efTime2', "passStartTime", "passEndTime", 'netDeploymentTime', 'scientificName', 'passNumber') 
# added a few metrics to quantify catch per unit effort such as 'passNumber', efish time, net deployment time 

# aggregate densities for each species group, pull out year and month from StartDate

fsh_dat_aggregate <- fsh_dat_fine %>%
  dplyr::select(!!c(all_of(my_grouping_vars), 'count')) %>%
  dplyr::group_by_at(dplyr::vars(all_of(my_grouping_vars))) %>%  
  dplyr::summarize(
    number_of_fish = sum(count),
    n_obs = dplyr::n()) %>%
  dplyr::mutate(
    year = startDate %>% lubridate::year(),
    month = startDate %>% lubridate::month()
  ) %>% dplyr::ungroup()


# some aquaticSiteType are NA, replace NAs if-based wildcarding namedLOcation 
# grepl is for wildcarding 
fsh_dat_aggregate$aquaticSiteType <- dplyr::if_else(is.na(fsh_dat_aggregate$aquaticSiteType),
                                                    if_else(grepl("fish.point", fsh_dat_aggregate$namedLocation), 'stream','lake'), fsh_dat_aggregate$aquaticSiteType)

# sampler type is also missing in a few cases, reaplace from eventID, with a wildcard
fsh_dat_aggregate$samplerType <- dplyr::if_else(is.na(fsh_dat_aggregate$samplerType), 
                                                dplyr::if_else(grepl("e-fisher", fsh_dat_aggregate$eventID), 'electrofisher', 
                                                               dplyr::if_else(grepl("gill", fsh_dat_aggregate$eventID), 'gill net', 'mini-fyke net')), fsh_dat_aggregate$samplerType) 


# make sure that e-fish samples have "" or NA for netset and netend times, this will leave actual missing data as na
# change netset/netend times to a datetime format if needed: lubridate::as_datetime(fsh_xx$netSetTime, format="%Y-%m-%d T %H:%M", tz="GMT")
## then calculate the netdeploymenttime (in hours) from netset and netend time
# to calculate pass duration (in mins) and also calculate average efish time (secs); also if efishtime is zero, --> na's 

fsh_aggregate_mod <- fsh_dat_aggregate %>% dplyr::mutate(netSetTime = dplyr::if_else(condition = samplerType ==  "electrofisher" | samplerType == "two electrofishers",
                                                                                     true = lubridate::as_datetime(NA), false = netSetTime), 
                                                         netEndTime = dplyr::if_else(condition = samplerType ==  "electrofisher" | samplerType == "two electrofishers",
                                                                                     true = lubridate::as_datetime(NA), false = netEndTime),
                                                         netDeploymentTime = dplyr::case_when(samplerType == grepl("electrofish", samplerType) ~ netDeploymentTime, 
                                                                                              is.na(netDeploymentTime) ~ as.numeric(difftime(netEndTime, netSetTime, tz = "GMT", units = "hours")), 
                                                                                              TRUE ~ netDeploymentTime),
                                                         mean_efishtime = base::rowMeans(dplyr::select(., c("efTime", "efTime2")), na.rm = T),
                                                         mean_efishtime = dplyr::case_when(mean_efishtime == 0 ~ NA_real_,TRUE ~ mean_efishtime))

# with the above changes, we have efish time and and net duration to calculate catch per unit effort before moving to wide format
# CPUE with efish time, calculated by = (total number of fish/average e-fish time in secs * 3600) as fish captured per 1-hr of e-fishing
# CPUE with gill nets and fykenets calculated by = (total number of fish/netDeployment time in hours * 24) as fish captured per 1-day-mong net deployment of e-fishing

fsh_aggregate_mod2 <- fsh_aggregate_mod %>% dplyr::mutate(CPUE = dplyr::if_else(condition = samplerType ==  "electrofisher" | samplerType == "two electrofishers",
                                                                                true = number_of_fish/mean_efishtime * 3600, 
                                                                                false = dplyr::if_else(condition = samplerType ==  "mini-fyke net" | samplerType == "gill net",
                                                                                                       true = number_of_fish/netDeploymentTime * 24, false = as.numeric(NA))))

# calculate abundance for all fish
fsh_dat_wide_total.obs2 = fsh_dat_wide_total.obs %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(abundance = dplyr::select(., `Agonostomus monticola`:`Xiphophorus hellerii`) %>% rowSums(.))

## species richness per row
fsh_dat_wide_total.obs3 = fsh_dat_wide_total.obs2 %>% 
  dplyr::mutate(species_richness = base::rowSums(dplyr::select(., `Agonostomus monticola`:`Xiphophorus hellerii`)  > 0))

# make wide for total observations without catch per unit efforts
fsh_dat_wide_total.obs <- fsh_aggregate_mod %>% 
  dplyr::group_by(year, month, siteID, namedLocation, reachID, fixedRandomReach, aquaticSiteType, samplerType) %>%
  tidyr::pivot_wider(names_from = scientificName, values_from = number_of_fish, names_repair = "unique",  values_fill = 0) %>%
  dplyr::select(-n_obs) # this col is misleading for wide format

# make wide for catch per unit efforts
fsh_dat_wide_CPUE <- fsh_aggregate_mod2 %>% 
  dplyr::group_by(year, month, siteID, namedLocation, reachID, fixedRandomReach, aquaticSiteType, samplerType) %>%
  tidyr::pivot_wider(names_from = scientificName, values_from = CPUE, names_repair = "unique",  values_fill = 0) %>%
  dplyr::select(-n_obs, -number_of_fish)

##############################################################################################################
# read the fieldData from NEON
fsh_fielddata <- readxl::read_excel(path = "fsh2020.xlsx", sheet = "fsh_fieldData")

# get fsh_fielddata from NEON and make a copy of it
fsh_f.data <- fsh_fielddata

####extract month, yr, and date into sep. columns, we might already have them
t_fsh_f.data <- t(as.data.frame(stringr::str_split( string = fsh_f.data$endDate, pattern = "-"))) # split date into yr, month, and date
head(t_fsh_f.data) # with date sep into three cols
fsh_f.data$sample_yr <- t_fsh_f.data[,1]
fsh_f.data$sample_mon <- t_fsh_f.data[,2]
fsh_f.data$sample_day <- t_fsh_f.data[,3]

# make a tibble
fsh_tbl <- tibble::as_tibble(fsh_f.data)

#how many temporal replicates for each site within each year?
fsh.data_site_yr_mon <- fsh_tbl %>%
  dplyr::group_by(siteID, sample_yr, sample_mon) %>%
  dplyr::select(siteID, sample_yr, sample_mon) %>%
  dplyr::summarise(count=n()) # shows number of sampling events for a give site across different yrs, and within a yr, across months
fsh.site_yr_mon <- dplyr::arrange(fsh.data_site_yr_mon, siteID, sample_yr, sample_mon) # count per month per year per site

# in each site, how many months were sampled in each yr?
fsh.data_site_yr <-  fsh.data_site_yr_mon %>%
  dplyr::group_by(siteID, sample_yr) %>%
  dplyr::select(siteID, sample_yr) %>%
  dplyr::summarise(count=n()) # no of months sampled per site per year

fsh.data_site_nyrs <- # number of yrs sampled for each site
  fsh.data_site_yr %>%
  dplyr::group_by(siteID) %>%
  dplyr::select(siteID) %>%
  dplyr::summarise(count=n())

##getting sites that have at least 3 years of data with at least 2 replicate-months within each year. two-month sampling bout per year
fsh.data_site_yr # this is the dataset

fsh.data.2mon.yr <- 
  fsh.data_site_yr %>%
  dplyr::filter(count >= 2) # what are the sites that at least have been sampled for 2 replicate months in a given year?.
# the count col gives you the number of months sampled in each year for each site
fsh.data.2mon.yr

fsh.data.2mon.3yr <- 
  fsh.data.2mon.yr %>%
  dplyr::group_by(siteID) %>%
  dplyr::select(siteID) %>%
  dplyr::summarise(count=n()) # per site, how many years have at least >=2 replicate months?
# count col gives you the number of years where each year was sampled at least 2 months

fsh.sites.2mon.3yr <- fsh.data.2mon.3yr %>%
  filter(count >= 3) #these are the sites with at least three yrs of data, each yr with at least 2-month replicates

# how many sites have at least three yrs of data collection, two bouts in each 
dplyr::n_distinct(fsh.sites.2mon.3yr$siteID)

save(fsh.sites.2mon.3yr, file = "fsh.sites.2mon.3yr.RData")

#############################

# concatenate year with month to make a separate column for bouts
fsh_dat_wide_CPUE2 <- fsh_dat_wide_CPUE %>% dplyr::mutate(bout = paste(year, month, sep = "_"))
save(fsh_dat_wide_CPUE2, file = "fsh_dat_wide_CPUE2.RData")

# filter out the sites that does not have temporal replicability
fsh_dat_wide_CPUE_3y_2mon <- fsh_dat_wide_CPUE2 %>% dplyr::filter(siteID %in% fsh.sites.2mon.3yr$siteID)
fsh_dat_wide_CPUE_3y_2mon
dplyr::n_distinct(fsh_dat_wide_CPUE_3y_2mon$siteID)
save(fsh_dat_wide_CPUE_3y_2mon, file = "fsh_dat_wide_CPUE_3y_2mon.RData")

# summarize by bout for each site for intra annual data
# first, omit all the unwanted columns
utils::View(colnames(fsh_dat_wide_CPUE_3y_2mon)) # get the column index
fsh_dat_wide_CPUE_3y_2mon2 <- fsh_dat_wide_CPUE_3y_2mon[,c(2, 22:23, 134, 25: 133)]

# now to summarize
fsh_dat_wide_CPUE_3y_2mon3 <- fsh_dat_wide_CPUE_3y_2mon2 %>% 
  dplyr::group_by(siteID, year, month, bout) %>% 
  dplyr::summarise_all(sum)

write.csv(x= fsh_dat_wide_CPUE_3y_2mon3, file = "fsh_dat_wide_CPUE_3y_2mon3.csv")

## change data structure
fsh_BAT_intra <- fsh_dat_wide_CPUE_3y_2mon3
fsh_BAT_intra$bout <- as.character(fsh_BAT_intra$bout)
fsh_BAT_intra$year <- as.double(fsh_BAT_intra$year)
fsh_BAT_intra$siteID <- as.factor(fsh_BAT_intra$siteID)

## save it to the desired format
saveRDS(fsh_BAT_intra, "fsh_BAT_intra2.rds") 
write.csv(x= fsh_BAT_intra, file = "fsh_BAT_intra2.csv")


########################################

# summarize by bout for each site for inter annual data
# first, omit all the unwanted columns
fsh_dat_wide_CPUE_3y_2mon4 <- fsh_dat_wide_CPUE_3y_2mon[,c(2, 22, 25: 133)]

# now to summarize
fsh_dat_wide_CPUE_3y_2mon5 <- fsh_dat_wide_CPUE_3y_2mon4 %>% 
  dplyr::group_by(siteID, year) %>% 
  dplyr::summarise_all(sum)

## change data str
fsh_BAT_inter <- fsh_dat_wide_CPUE_3y_2mon5
fsh_BAT_inter$year <- as.double(fsh_BAT_inter$year)
fsh_BAT_inter$siteID <- as.factor(fsh_BAT_inter$siteID)

## save it to the desired format
saveRDS(fsh_BAT_inter, "fsh_BAT_inter2.rds") 
write.csv(x= fsh_BAT_inter, file = "fsh_BAT_inter.csv")


###################################

# for codyn analyses
fsh_codyn <- fsh_aggregate_mod2[, c(2, 6, 21, 28)]

fsh_codyn$scientificName <- as.character(fsh_codyn$scientificName)
fsh_codyn$siteID <- as.character(fsh_codyn$siteID)

# change col names
fsh_codyn <- fsh_codyn %>%  dplyr::rename(abundance = CPUE) %>% dplyr::rename(collectDate = endDate)
  
# first, convert from datetime to date only
fsh_codyn$collectDate <- as.Date.character(fsh_codyn$collectDate, format= "%Y-%m-%d")

# now from date to chr
fsh_codyn$collectDate <- as.character(fsh_codyn$collectDate, format= "%Y-%m-%d")

# Save as an .rds file
# use readRDS to read it
saveRDS(object = fsh_codyn, file = "fsh_codyn2.rds")
write.csv(fsh_codyn, "fsh_codyn2.csv")



##############################################################################################
# Stephanie Parker (sparker@battelleecology.org)
# last updated 2019-11-04
# uses esokol's 'get_and_format_NEON_inv_data.R'
##############################################################################################

rm(list = ls())
gc()

# set options
options(stringsAsFactors = FALSE)

# load packages
library(neonUtilities)
library(tidyverse)


#################################################################################
# fish dpid
my_dpid <- 'DP1.20107.001'

# some aquatic sites
#my_site_list <- c('POSE', 'ARIK')

# get taxon table from API, may take a few minutes to load
full_taxon_table <- neonUtilities::getTaxonTable('FISH') 

#############################
# -- make ordered taxon_rank_list for a reference (subspecies is smallest rank, kingdom is largest)
taxon_rank_list_ordered <- c('kingdom', 'subkingdom',
  'infrakingdom', 'superphylum', 'phylum', 'subphylum', 'infraphylum',
  'superdivision', 'division', 'subdivision', 'infradivision', 'parvdivision',
  'superclass', 'class', 'subclass', 'infraclass', 'superorder',
  'order', 'suborder', 'infraorder', 'section', 'subsection',
  'superfamily', 'family', 'subfamily', 'tribe', 'subtribe',
  'genus','subgenus','speciesGroup','species','subspecies') %>% rev() #probably don't need all these levels, revise later
# #############################

# get data via api
all_tabs <- neonUtilities::loadByProduct(
  dpID = my_dpid,
  site = "all",
  check.size = FALSE)


# join field and pass tables
fsh_dat1 <- left_join(all_tabs$fsh_perPass, all_tabs$fsh_fieldData, 
  by = c('reachID')) %>% 
  filter(is.na(samplingImpractical)) #remove records where fish couldn't be collected

# get rid of dupe col names and .x suffix
fsh_dat1 <- fsh_dat1[,!grepl('\\.y',names(fsh_dat1))]
names(fsh_dat1) <- gsub('\\.x','',names(fsh_dat1))

# individual fish counts
fsh_dat_indiv <- left_join(all_tabs$fsh_perFish, fsh_dat1, by = "eventID") 

# get rid of dupe col names and .x suffix
fsh_dat_indiv <- fsh_dat_indiv[,!grepl('\\.y',names(fsh_dat_indiv))]
names(fsh_dat_indiv) <- gsub('\\.x','',names(fsh_dat_indiv))

# fill in missing reachID
fsh_dat_indiv$reachID <- ifelse(is.na(fsh_dat_indiv$reachID), 
  substr(fsh_dat_indiv$eventID, 1, 16), fsh_dat_indiv$reachID) 


# bulk fish counts
fsh_dat_bulk <- left_join(all_tabs$fsh_bulkCount, fsh_dat1, by = "eventID")

# get rid of dupe col names and .x suffix
fsh_dat_bulk <- fsh_dat_bulk[,!grepl('\\.y',names(fsh_dat_bulk))]
names(fsh_dat_bulk) <- gsub('\\.x','',names(fsh_dat_bulk))

#fill in missing reachID
fsh_dat_bulk$reachID <- ifelse(is.na(fsh_dat_bulk$reachID), 
  substr(fsh_dat_bulk$eventID, 1, 16), fsh_dat_bulk$reachID) 


# combine indiv and bulk counts
fsh_dat <- bind_rows(fsh_dat_indiv, fsh_dat_bulk)

# add count = 1 for indiv data
fsh_dat$count <- ifelse(is.na(fsh_dat$bulkFishCount), 1, fsh_dat$bulkFishCount)

#fill in missing startDate, siteID, NEON will QC portal data to find out why this is missing
fsh_dat$startDate <- ifelse(is.na(fsh_dat$startDate), 
  substr(fsh_dat$passStartTime, 1, 10), fsh_dat$startDate) 
fsh_dat$siteID <- ifelse(is.na(fsh_dat$siteID), 
  substr(fsh_dat$eventID, 1, 4), fsh_dat$siteID) 

# get species and finer res
fsh_dat$taxonRank_ordered <- factor(
  fsh_dat$taxonRank,
  levels = taxon_rank_list_ordered,
  ordered = TRUE) 

# get all records that have rank <= species
fsh_dat_fine <- fsh_dat %>%
  filter(taxonRank_ordered <= 'species')


# grouping vars for aggregating density measurements 
my_grouping_vars <- c('domainID','siteID','aquaticSiteType','namedLocation',
  'startDate','reachID','eventID','samplerType',
  'fixedRandomReach','measuredReachLength','efTime','netDeploymentTime',
  'scientificName') #could add in 'passNumber' if you want to

# aggregate densities for each species group, pull out year and month from collectDate
fsh_dat_aggregate <- fsh_dat_fine %>%
  select(!!c(my_grouping_vars, 'count')) %>%
  group_by_at(vars(my_grouping_vars)) %>%
  summarize(
    number_of_fish = sum(count),
    n_obs = n()) %>% 
  mutate(
    year = startDate %>% lubridate::year(),
    month = startDate %>% lubridate::month()
  ) %>% ungroup()

# may want to calculate as catch per unit effort before moving to wide format


# make wide
fsh_dat_wide <- fsh_dat_aggregate %>% 
  group_by(year, month, siteID, namedLocation, reachID, fixedRandomReach) %>%
  spread(scientificName, number_of_fish, fill = 0)

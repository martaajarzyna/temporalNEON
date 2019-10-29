rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

library(neonUtilities)
library(tidyverse)

#################################################################################
# -- user provided rank for desired taxonomic resolution
my_cutoff_taxon_rank <- 'genus'

#################################################################################
# macroinvert dpid
my_dpid <- 'DP1.20120.001'

# some aquatic sites
my_site_list <- c('COMO', 'ARIK')

# get taxon table from API, may take a few minutes to load
full_taxon_table <- neonUtilities::getTaxonTable('MACROINVERTEBRATE')

#############################
# -- make ordered taxon_rank_list for a reference (subspecies is smallest rank, kingdom is largest)
taxon_rank_list_ordered <- c('kingdom', 'subkingdom',
                     'infrakingdom', 'superphylum', 'phylum', 'subphylum', 'infraphylum',
                     'superdivision', 'division', 'subdivision', 'infradivision', 'parvdivision',
                     'superclass', 'class', 'subclass', 'infraclass', 'superorder',
                     'order', 'suborder', 'infraorder', 'section', 'subsection',
                     'superfamily', 'family', 'subfamily', 'tribe', 'subtribe',
                     'genus','subgenus','speciesGroup','species','subspecies') %>% rev()
# #############################

# get data via api
all_tabs <- neonUtilities::loadByProduct(
  dpID = my_dpid,
  site = my_site_list,
  check.size = TRUE)

# join taxon table and field table -- modified from Mariana's code
inv_dat <- left_join(all_tabs$inv_taxonomyProcessed, all_tabs$inv_fieldData, 
                     by = c('sampleID')) %>% 
  mutate(den = estimatedTotalCount/benthicArea) %>% 
  mutate(scientificName = forcats::fct_explicit_na(scientificName)) %>%
  filter(sampleCondition == "condition OK")

# get rid of dupe col names and .x suffix
inv_dat <- inv_dat[,!grepl('\\.y',names(inv_dat))]
names(inv_dat) <- gsub('\\.x','',names(inv_dat))


# get genus and finer res
inv_dat$taxonRank_ordered <- factor(
  inv_dat$taxonRank,
  levels = taxon_rank_list_ordered,
  ordered = TRUE) 

# get all records that have rank <= genus, where genus is not NA or blank
inv_dat_fine <- inv_dat %>%
  filter(taxonRank_ordered <= 'genus') %>%
  filter(!is.na(genus), genus != '')

# grouping vars for aggregating density measurements 
my_grouping_vars <- c('domainID','siteID','namedLocation','collectDate','sampleID', 
                      'aquaticSiteType','eventID','habitatType',
                      'substratumSizeClass', 'benthicArea', 'genus')

# aggregate densities for each genus group, pull out year and month from collectDate
inv_dat_aggregate <- inv_dat_fine %>%
  select(!!c(my_grouping_vars, 'den')) %>%
  group_by_at(vars(my_grouping_vars)) %>%
  summarize(
    den_sum = sum(den),
    n_obs = n()) %>% 
  mutate(
    year = collectDate %>% lubridate::year(),
    month = collectDate %>% lubridate::month()
  ) %>% ungroup()

# make wide
inv_dat_wide <- inv_dat_aggregate %>% 
  group_by(year, month, siteID, namedLocation, habitatType, substratumSizeClass) %>%
  spread(genus, den_sum, fill = 0)

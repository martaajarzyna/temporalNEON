#Load packages
require(tidyverse)
require(dplyr)
# remotes::install_github("EDIorg/ecocomDP")
require(ecocomDP)
#install.packages("neonUtilities")
require(neonUtilities)
require(lubridate)
library(devtools)
library(stringr)
library(here)
library(forcats)
'%>%' <- dplyr::`%>%`


###################DOWNLOAD AND PROCESS RAW DATA 

################### small mammals
# https://data.neonscience.org/data-products/DP1.10072.001
map_neon_data_to_ecocomDP.SMALL.MAMMAL <- function(
  neon.data.product.id = "DP1.10072.001",
  ...){
  # Author: Marta Jarzyna
  ### This code downloads, unzips and restructures small mammal data (raw abundances) from NEON
  ### for each NEON site, we provide mammal raw abundance aggregated per day (across all years),
  ### per month (bout) (across all years), and per year
  ### By Marta Jarzyna
  d = neonUtilities::loadByProduct(dpID = neon.data.product.id, site = "all", package = "basic")
  saveRDS(d, file = "data-raw/mammal_raw.rds")
  # d = readRDS(file = "data-raw/mammal_raw.rds")

  dat.mam = dplyr::mutate(d$mam_pertrapnight,
                          collectDate = lubridate::ymd(collectDate),
                          year = lubridate::year(collectDate),
                          month = lubridate::month(collectDate),
                          day = lubridate::day(collectDate),
                          # Since data collection usually takes place on several consecutive
                          # days within a given month, we consider each month to be a separate bout
                          bout = paste(year, month, sep = "_")) %>%
    tibble::as_tibble()
  table(dat.mam$plotType) # all distributed
  
  # remove 2020 because of incomplete sampling (covid-19) and Alaska (domains 18 & 19), Puerto Rico (D04), and Hawaii (D20)
  dat.mam <- dplyr::filter(dat.mam, year != 2020)
  dat.mam <- dplyr::filter(dat.mam, domainID != "D18")
  dat.mam <- dplyr::filter(dat.mam, domainID != "D19")
  dat.mam <- dplyr::filter(dat.mam, domainID != "D04") 
  dat.mam <- dplyr::filter(dat.mam, domainID != "D20")
  
  dat1 <- dat.mam %>%
    group_by(siteID, year, month) %>%
    select(siteID, year, month) %>%
    summarise(count=n()) %>%
    arrange(siteID, year, month) %>%
    ungroup
  
  dat2 <- dat1 %>%
    group_by(siteID, year) %>%
    select(siteID, year) %>%
    summarise(count=n()) %>%
    ungroup()
  
  # select sites with >=4 years of data and >=4 replicate-months (bouts) within each year
  dat3 <- dat2 %>%
    filter(count >= 4) %>%
  ungroup()
  
  dat4 <- 
    dat3 %>%
    group_by(siteID) %>%
    select(siteID) %>%
    summarise(count=n())
  
  dat5 <-
    dat4 %>%
    filter(count >= 4) #these are the final sites
  saveRDS(dat5, file="data/mammal_selsites.rds") #save sites that meet criteria
  
  ### Remove NA (no captures in the trap) altogether at this point
  dat.mam <- dplyr::filter(dat.mam, scientificName != "", !is.na(scientificName))

  ### Remove species that are bycatch (non-target), dead, or escapted while processing
  ### remove when fate <- 'dead' = dead, 'escaped' = escaped while handling, 'nontarget' = released, non-target species,
  dat.mam <- dplyr::filter(dat.mam, !fate %in% c("dead", "escaped", "nontarget"))

  ### Remove recaptures -- Y and U (unknown); only retain N
  # table(dat.mam$recapture)
  # sum(is.na(dat.mam$recapture))
  dat.mam <- dplyr::filter(dat.mam, recapture == "N" | is.na(recapture))

  ### Remove captures not id'ed to species 
  dat.mam <- dplyr::filter(dat.mam, !grepl(pattern = "sp.", scientificName))
  dat.mam <- dplyr::filter(dat.mam, !grepl(pattern = "/", scientificName))
  
  ### Retain sites that meet the replicatioon criteria
  dat.mam <- dplyr::filter(dat.mam, siteID %in% dat5$siteID)
  
  ### Get raw abundances per day
  data_small_mammal <- dat.mam %>%
    dplyr::select(domainID, siteID, plotID, namedLocation, nlcdClass, decimalLatitude,
           decimalLongitude, geodeticDatum, coordinateUncertainty,
           elevation, elevationUncertainty, trapCoordinate, trapStatus, year, month,
           day, bout, taxonID, scientificName, taxonRank, identificationReferences,
           nativeStatusCode, sex, lifeStage,
           pregnancyStatus, tailLength, totalLength, weight) %>%
    dplyr::distinct()
  
  ### Get mean latitude and longitude per neon site (certain measurements within a site have slightly different lat and lon, but our analysis is on the site level)
  mean.coord <- data_small_mammal %>%
    group_by(siteID) %>%
    summarize(lat = mean(decimalLatitude), lon = mean(decimalLongitude, na.rm=TRUE)) 
  data_small_mammal <- left_join(data_small_mammal, mean.coord, by="siteID") 
  
  return(data_small_mammal)
}

data_small_mammal = map_neon_data_to_ecocomDP.SMALL.MAMMAL(neon.data.product.id = "DP1.10072.001")
#save processed data
saveRDS(data_small_mammal, file = "data/mammal_processed.rds")
#data_small_mammal <- readRDS(file = "data/mammal_processed.rds")




################### ground beetles
map_neon_data_to_ecocomDP.BEETLE <- function(
  neon.data.product.id = "DP1.10022.001",
  ...){
  # author: Kari Norman
  
  # download data
  beetles_raw <- neonUtilities::loadByProduct(dpID = neon.data.product.id, site = "all", package = "basic")
  saveRDS(beetles_raw, file = "data-raw/beetle_raw.rds")
  # beetles_raw = readRDS(file = "data-raw/beetle_raw.rds")
  
  # helper function to calculate mode of a column/vector
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  ### Clean Up Sample Data ###
  
  data <- beetles_raw$bet_fielddata %>%
    dplyr::filter(sampleCollected == "Y") %>% #there's an entry for every trap, whether or not they got samples, only want ones with samples
    dplyr::select(sampleID, namedLocation, domainID, siteID, plotID, trapID, setDate, collectDate, eventID, trappingDays) %>%
    #eventID's are inconsistently separated by periods, so remove
    dplyr::mutate(eventID = str_remove_all(eventID, "[.]")) %>%
    dplyr::mutate(trappingDays = interval(ymd(setDate), ymd(collectDate)) %/% days(1))
  
  
  #Find the traps that have multiple collectDates/bouts for the same setDate
  #need to calculate their trap days from the previous collectDate, not the setDate
  adjTrappingDays <- data %>%
    select(plotID, trapID, setDate, collectDate, trappingDays, eventID) %>%
    group_by_at(vars(-collectDate, -trappingDays, -eventID)) %>%
    filter(n_distinct(collectDate) > 1) %>%
    group_by(plotID, trapID, setDate) %>%
    mutate(diffTrappingDays = trappingDays - min(trappingDays)) %>%
    mutate(adjTrappingDays = case_when(
      diffTrappingDays == 0 ~ trappingDays,
      TRUE ~ diffTrappingDays
    )) %>%
    select(-c(trappingDays, diffTrappingDays))
  
  data <- data %>%
    #update with adjusted trapping days where needed
    left_join(adjTrappingDays) %>%
    mutate(trappingDays = case_when(
      !is.na(adjTrappingDays) ~ adjTrappingDays,
      TRUE ~ trappingDays
    )) %>%
    select(-adjTrappingDays, -setDate) %>%
    #for some eventID's (bouts) collection happened over two days,
    #change collectDate to the date that majority of traps were collected on
    group_by(eventID) %>%
    mutate(collectDate = Mode(collectDate)) %>%
    ungroup() %>%
    #there are also some sites for which all traps were set and collect on the same date, but have multiple eventID's
    #we want to consider that as all one bout so we'll just create a new ID based on the site and collectDate
    unite(boutID, siteID, collectDate, remove = FALSE) %>%
    select(-eventID) %>%
    #and join to sample data
    left_join(beetles_raw$bet_sorting %>%
                filter(sampleType %in% c("carabid", "other carabid")) %>% #only want carabid samples, not bycatch
                select(sampleID, subsampleID, sampleType, taxonID, scientificName, taxonRank, individualCount,identificationQualifier),
              by = "sampleID") %>%
    filter(!is.na(subsampleID)) #even though they were marked a sampled, some collection times don't acutally have any samples
  
  ### Join taxonomic data from pinning with the sorting data ###
  
  # Replace sorting taxon info with pinning taxon info (people that pin specimens are more experienced with taxonomy), where available
  data_pin <- data %>%
    left_join(beetles_raw$bet_parataxonomistID %>% select(subsampleID, individualID, taxonID, scientificName, taxonRank,identificationQualifier), by = "subsampleID") %>%
    mutate_if(is.factor, as.character) %>%
    mutate(taxonID = ifelse(is.na(taxonID.y), taxonID.x, taxonID.y)) %>%
    mutate(taxonRank = ifelse(is.na(taxonRank.y), taxonRank.x, taxonRank.y)) %>%
    mutate(scientificName = ifelse(is.na(scientificName.y), scientificName.x, scientificName.y)) %>%
    mutate(identificationSource = ifelse(is.na(scientificName.y), "sort", "pin")) %>%
    mutate (identificationQualifier = ifelse(is.na(taxonID.y), identificationQualifier.x, identificationQualifier.y)) %>%
    select(-ends_with(".x"), -ends_with(".y"))
  
  #some subsamples weren't fully ID'd by the pinners, so we have to recover the unpinned-individuals
  lost_indv <- data_pin %>% 
    filter(!is.na(individualID)) %>%
    group_by(subsampleID, individualCount) %>%
    summarise(n_ided = n_distinct(individualID)) %>% 
    filter(n_ided < individualCount) %>%
    mutate(unidentifiedCount = individualCount - n_ided) %>%
    select(subsampleID, individualCount = unidentifiedCount) %>%
    left_join(data %>% select(-individualCount), by = "subsampleID") %>%
    mutate(identificationSource = "sort")
  
  #add unpinned-individuals back to the pinned id's, adjust the individual counts so pinned individuals have a count of 1
  data_pin <- data_pin %>%
    dplyr::mutate(individualCount = ifelse(identificationSource == "sort", individualCount, 1)) %>%
    dplyr::bind_rows(lost_indv)
  
  
  ### Join expert data to existing pinning and sorting data ###
  
  #There are ~10 individualID's for which experts ID'd more than one species (not all experts agreed), we want to exclude those expert ID's as per Katie Levan's suggestion
  ex_expert_id <- beetles_raw$bet_expertTaxonomistIDProcessed %>% 
    group_by(individualID) %>% 
    filter(n_distinct(taxonID) > 1) %>% 
    pull(individualID)
  
  # Add expert taxonomy info, where available
  data_expert <- left_join(data_pin, 
                           select(beetles_raw$bet_expertTaxonomistIDProcessed,
                                  individualID,taxonID,scientificName,taxonRank,identificationQualifier) %>%
                             filter(!individualID %in% ex_expert_id), #exclude ID's that have unresolved expert taxonomy
                           by = 'individualID', na_matches = "never") %>% distinct()
  
  # Replacement old taxon info with expert info, where available
  # NOTE - This is repetitive with the code snippet above, and if you want to do it this way you can just combine the calls into one chunk. BUT, you may
  #     want to do more than this, as it is *just* a replacement of IDs for individual beetles that an expert identified. If the expert identified
  #           a sample as COLSP6 instead of CARSP14, though, then all CARSP14 from that trap on that date should probably be updated to COLSP6…
  data_expert <- data_expert %>%
    mutate_if(is.factor, as.character) %>%
    mutate(taxonID = ifelse(is.na(taxonID.y), taxonID.x, taxonID.y)) %>%
    mutate(taxonRank = ifelse(is.na(taxonRank.y), taxonRank.x, taxonRank.y)) %>%
    mutate(scientificName = ifelse(is.na(scientificName.y), scientificName.x, scientificName.y)) %>%
    mutate(identificationSource = ifelse(is.na(scientificName.y), identificationSource, "expert")) %>%
    mutate (identificationQualifier = ifelse(is.na(taxonID.y), identificationQualifier.x, identificationQualifier.y)) %>%
    select(-ends_with(".x"), -ends_with(".y"))
  
  beetles_data <- data_expert
  #usethis::use_data(beetles_data)
  
  
  ### Retain captures id'ed to species 
  beetles_data <- dplyr::filter(beetles_data, taxonRank == "species")
  
  #Get raw counts table
  beetles_counts <- beetles_data %>%
    select(-c(subsampleID, sampleType, individualID, identificationSource, identificationQualifier)) %>%
    group_by_at(vars(-individualCount)) %>%
    summarise(count = sum(individualCount)) %>%
    mutate(year = substr(boutID, 6,9)) %>%
    ungroup()
  
  #usethis::use_data(beetles_counts)
  
  # remove 2020 because of incomplete sampling (covid-19) and Alaska (domains 18 & 19), Puerto Rico (D04), and Hawaii (D20)
  beetles_counts <- beetles_counts %>%
    dplyr::filter(year != 2020) %>%
    dplyr::filter(domainID != "D18") %>%
    dplyr::filter(domainID != "D19") %>%
    dplyr::filter(domainID != "D04") %>%
    dplyr::filter(domainID != "D20")
  
  ### Obtain sites that met our criteria (i.e., sampled for at least 4 years with 3 replicates) ###
  included_sites <- beetles_counts %>% 
    mutate(year = year(collectDate)) %>% 
    select(siteID, boutID, year) %>% 
    distinct() %>% count(siteID, year) %>% 
    filter(n > 2) %>% 
    group_by(siteID) %>% 
    filter(n_distinct(year) > 3) %>% 
    pull(siteID) %>% unique()
  
  beetles_counts <- beetles_counts %>%
    dplyr::filter(siteID %in% included_sites)
  
  # get relevant location info from the data
  table_location_raw <- beetles_raw$bet_fielddata %>%
    dplyr::select(domainID, siteID, namedLocation, plotType, nlcdClass, decimalLatitude,
                  decimalLongitude, geodeticDatum, coordinateUncertainty,
                  elevation, elevationUncertainty) %>%
    dplyr::distinct()
  
  
  data_beetle <-  dplyr::left_join(beetles_counts, table_location_raw,
                                   by = c("domainID", "siteID", "namedLocation"))
  all(paste0(data_beetle$plotID, ".basePlot.bet") == data_beetle$namedLocation)
  # data_beetle = dplyr::select(data_beetle, -namedLocation)
  
  ### Get mean latitude and longitude per neon site (certain measurements within a site have slightly different lat and lon, but our analysis is on the site level)
  mean.coord <- data_beetle %>%
    group_by(siteID) %>%
    summarize(lat = mean(decimalLatitude), lon = mean(decimalLongitude, na.rm=TRUE)) 
  data_beetle <- left_join(data_beetle, mean.coord, by="siteID") 
  
  return(data_beetle)
}

data_beetle = map_neon_data_to_ecocomDP.BEETLE(neon.data.product.id = "DP1.10022.001")
#save processed data
saveRDS(data_beetle, file = "data/beetle_processed.rds")
#data_beetle <- readRDS(file = "data/beetle_processed.rds")




################### fish
map_neon_data_to_ecocomDP.FISH <- function(
  neon.data.product.id = "DP1.20107.001",
  ...){
  # Authors: Stephanie Parker, Thilina Surasinghe (sparker@battelleecology.org, tsurasinghe@bridgew.edu)
  
  # get taxon table from API, may take a few minutes to load
  # Fish electrofishing, gill netting, and fyke netting counts
  # http://data.neonscience.org/data-product-view?dpCode=DP1.20107.001
  full_taxon_fish <- neonUtilities::getTaxonTable(taxonType = 'FISH', recordReturnLimit = NA, stream = "true")
  
  # -- make ordered taxon_rank list for a reference (subspecies is smallest rank, kingdom is largest)
  
  # a much simple table with useful levels of taxonomic resolution
  taxon_rank_fish <- c('superclass', 'class', 'subclass', 'infraclass', 'superorder',
                       'order', 'suborder', 'infraorder', 'section', 'subsection',
                       'superfamily', 'family', 'subfamily', 'tribe', 'subtribe', 'genus',
                       'subgenus','speciesGroup','species','subspecies') %>% rev()
  
  
  # get data FISH via api -- will take a while --  package = "basic" is also possible
  all_fish <- neonUtilities::loadByProduct(dpID = neon.data.product.id,
                                           site = "all", startdate = NA, enddate = NA,
                                           package = "expanded", avg = "all", check.size = FALSE)
  saveRDS(all_fish, file = "data-raw/fish_raw.rds")
  # all_fish = readRDS(file = "data-raw/fish_raw.rds")
  
  # this joins reach-level data, you can just join by any two of these, and the results are the same, dropping the unique ID
  # join field data table with the pass tables from all_fish
  fsh_dat1 <- dplyr::left_join(x = all_fish$fsh_perPass[-1], y = all_fish$fsh_fieldData[-1], by = c("reachID")) %>% 
    dplyr::filter(is.na(samplingImpractical) | samplingImpractical == "") %>%  #remove records where fish couldn't be collected, both na's and blank data is kept
    tibble::as_tibble()
  
  # get rid of duplicate col names and .x suffix
  fsh_dat1 <- fsh_dat1[,!grepl('\\.y',names(fsh_dat1))]
  names(fsh_dat1) <- gsub('\\.x','',names(fsh_dat1))
  
  # add individual fish counts; perFish data = individual level per each specimen captured
  # this joins reach-level field data with individual fish measures
  fsh_dat_indiv <- dplyr::left_join(all_fish$fsh_perFish, fsh_dat1, by = c("eventID")) %>%
    tibble::as_tibble()
  
  # get rid of dupe col names and .x suffix
  fsh_dat_indiv <- fsh_dat_indiv[,!grepl('\\.y',names(fsh_dat_indiv))]
  names(fsh_dat_indiv) <- gsub('\\.x','',names(fsh_dat_indiv))
  
  # fill in missing reachID with event ID
  fsh_dat_indiv$reachID <- ifelse(is.na(fsh_dat_indiv$reachID), 
                                  substr(fsh_dat_indiv$eventID, 1, 16), fsh_dat_indiv$reachID) 
  
  # add bulk fish counts; bulk count data = The number of fish counted during bulk processing in each pass
  # this will join all the reach-level field data with the species data with bulk counts
  fsh_dat_bulk <- dplyr::left_join(all_fish$fsh_bulkCount, fsh_dat1, by = c("eventID")) %>%
    tibble::as_tibble()
  
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
  
  # some aquaticSiteType are NA, replace NAs if-based wildcarding namedLOcation 
  # grepl is for wildcarding 
  fsh_dat$aquaticSiteType <- dplyr::if_else(is.na(fsh_dat$aquaticSiteType),
                                            dplyr::if_else(grepl("fish.point", fsh_dat$namedLocation), 'stream','lake'), 
                                            fsh_dat$aquaticSiteType)
  
  # sampler type is also missing in a few cases, reaplace from eventID, with a wildcard
  fsh_dat$samplerType <- dplyr::if_else(is.na(fsh_dat$samplerType), 
                                        dplyr::if_else(grepl("e-fisher", fsh_dat$eventID), 'electrofisher', 
                                                       dplyr::if_else(grepl("gill", fsh_dat$eventID), 'gill net', 
                                                                      'mini-fyke net')), fsh_dat$samplerType) 
  
  # in case the bulk count data set does not have a column on taxonomic rank, as such, need to add it here.
  # in scientificName column-- species identified below species level (genus, family, order, phylum)-- appear as sp. or spp.
  # both sp. and spp. identifications should be marked low-res identifications, aka above species level in ranking  
  # and exclude from this analyses
  fsh_dat2 <- fsh_dat %>% 
    dplyr::mutate(taxonRank = dplyr::case_when(stringr::str_detect(string = scientificName, pattern = " spp\\.$") ~ "not_sp_level1", 
                                               stringr::str_detect(string = scientificName, pattern = " sp\\.$") ~ "not_sp_level2", 
                                               TRUE ~ "species")) 
  
  # get all records that have rank <= species; if bulkCount data set has taxonRank, skip this
  fsh_dat_fine <- fsh_dat2 %>% dplyr::filter(taxonRank == "species")
  
  # drop the lake fish data and keep stream data only
  fsh_dat_stm <- fsh_dat_fine %>% dplyr::filter(aquaticSiteType == "stream")
  
  # some species are listed as two or more sub-species, go down to species level
  fsh_dat_stm <- fsh_dat_stm %>% dplyr::mutate(scientificName = stringr::word(string = scientificName, start = 1, end = 2))
  
  # there is a spelling error in Oncorhynchus clarkii, correct it
  fsh_dat_stm <-  fsh_dat_stm %>% dplyr::mutate(scientificName = stringr::str_replace_all(scientificName, pattern = "Oncorhynchus clarki.*$", 
                                                                                          replacement = "Oncorhynchus clarkii"))
  
  # aggregate densities for each species group, pull out year and month from StartDate
  
  fsh_dat_aggregate <-fsh_dat_stm %>%
    dplyr::select(!!c('domainID','siteID','aquaticSiteType','namedLocation', 'decimalLatitude', 'decimalLongitude',
                      'startDate', 'endDate', 'reachID','eventID','samplerType', 'fixedRandomReach','measuredReachLength','efTime', 
                      'efTime2', "passStartTime", "passEndTime", 'scientificName', 'passNumber', 'count')) %>%
    dplyr::group_by_at(dplyr::vars('domainID','siteID','aquaticSiteType','namedLocation','decimalLatitude', 'decimalLongitude',
                                   'startDate', 'endDate', 'reachID','eventID','samplerType', 'fixedRandomReach','measuredReachLength','efTime', 
                                   'efTime2', "passStartTime", "passEndTime", 'scientificName', 'passNumber')) %>%  
    dplyr::summarize(
      number_of_fish = sum(count),
      n_obs = dplyr::n()) %>%
    dplyr::mutate(
      year = startDate %>% lubridate::year(),
      month = startDate %>% lubridate::month()
    ) %>% dplyr::ungroup()
  
  
  # to calculate pass duration (in mins) and also calculate average efish time (secs); also if efishtime is zero, --> na's 
  fsh_aggregate_mod <- fsh_dat_aggregate %>% 
    dplyr::mutate(mean_efishtime = base::rowMeans(dplyr::select(., c("efTime", "efTime2")), na.rm = T),
                  mean_efishtime = dplyr::case_when(mean_efishtime == 0 ~ NA_real_, TRUE ~ mean_efishtime))
  
  # with the above changes, we have efish time to calculate catch per unit effort before moving to wide format
  # CPUE with efish time, calculated by = (total number of fish/average e-fish time in secs * 3600) as fish captured per 1-hr of e-fishing
  data_fish <- fsh_aggregate_mod %>% 
    dplyr::mutate(cpue = number_of_fish/mean_efishtime * 3600)
  
  
  ### Remove 2020 because of incomplete sampling (covid-19) and Alaska (domains 18 & 19), Puerto Rico (D04), and Hawaii (D20)
  data_fish <- data_fish %>%
    dplyr::filter(year != 2020) %>%
    dplyr::filter(domainID != "D18") %>%
    dplyr::filter(domainID != "D19") %>%
    dplyr::filter(domainID != "D04") %>%
    dplyr::filter(domainID != "D20")
  
  # concatenate year with month to make a separate column for bouts
  data_fish <- data_fish %>% dplyr::mutate(bout = paste(year, month, sep = "_"))
  
  ### Select sites that meet our criteria (at least 3 years of data with at least 2 replicate-months within each year)
   dat1 <- data_fish %>%
    group_by(siteID, year, month) %>%
    select(siteID, year, month) %>%
    summarise(count=n()) %>%
    arrange(siteID, year, month) %>%
    ungroup
  
  dat2 <- dat1 %>%
    group_by(siteID, year) %>%
    select(siteID, year) %>%
    summarise(count=n()) %>%
    ungroup()
  
  # select sites with >=3 years of data and >=2 replicate-months (bouts) within each year
  dat3 <- dat2 %>%
    filter(count >= 2) %>%
    ungroup()
  
  dat4 <- 
    dat3 %>%
    group_by(siteID) %>%
    select(siteID) %>%
    summarise(count=n())
  
  dat5 <-
    dat4 %>%
    filter(count >= 3) #these are the final sites
  
  data_fish <- data_fish %>%
    dplyr::filter(siteID %in% dat5$siteID)
  
  # concatenate siteID with year to make a separate column for combinations of site+year that might not meet criteria
  # because within final sites there still might be years that do not meet criteria of >1 bouts
  data_fish <- data_fish %>% dplyr::mutate(siteyear = paste(siteID, year, sep = "_"))
  
  dat6 <- data_fish %>%
    group_by(siteID, year, siteyear, bout) %>%
    select(siteID, year, siteyear, bout) %>%
    summarise(count=n()) %>%
    arrange(siteID, year, bout) %>%
    ungroup
  
  dat7 <- dat6 %>%
    group_by(siteID, year, siteyear) %>%
    select(siteID, year, siteyear) %>%
    summarise(count=n()) %>%
    ungroup()
  
  dat7 <- dat7 %>%
    filter(count >= 2) %>%
    ungroup() #this removes one site-year combination that has only 1 bout
  
  data_fish <- data_fish %>%
    dplyr::filter(siteyear %in% dat7$siteyear)
  
  ### Get mean latitude and longitude per neon site (certain measurements within a site have slightly different lat and lon, but our analysis is on the site level)
  mean.coord <- data_fish %>%
    group_by(siteID) %>%
    summarize(lat = mean(decimalLatitude, na.rm=TRUE), lon = mean(decimalLongitude, na.rm=TRUE)) 
  data_fish <- left_join(data_fish, mean.coord, by="siteID") 
  
  return(data_fish)
}

data_fish = map_neon_data_to_ecocomDP.FISH(neon.data.product.id = "DP1.10022.001")
#save processed data
saveRDS(data_fish, file = "data/fish_processed.rds")
#data_fish <- readRDS(file = "data/fish_processed.rds")




################### aquatic macroinvertebrates
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

map_neon_data_to_ecocomDP.MACROINVERTEBRATE <- function(
  neon.data.product.id ="DP1.20120.001",
  ...){
  # Authors: Mariana Perez Rocha, Matt Helmus
  # The below code used to download the raw macroinverts files is not run during the package build. Don't set `eval = TRUE` for a quick render!
  
    inv_allTabs <- neonUtilities::loadByProduct(dpID = neon.data.product.id, 
                                 package = "expanded", nCores = 4, check.size = FALSE)
  
  # saveRDS(inv_allTabs, file = "data-raw/macroinvert_raw.rds")
  # inv_allTabs = readRDS("data-raw/macroinvert_raw.rds")
  
  ### Add taxonomy
  # make ordered taxon_rank_list for a reference (subspecies is smallest rank, kingdom is largest)
  taxon_rank_list_ordered <- c('kingdom', 'subkingdom','infrakingdom', 'superphylum', 'phylum', 'subphylum', 'infraphylum',
                               'superdivision', 'division', 'subdivision', 'infradivision', 'parvdivision','superclass', 'subclass',
                               'superorder', 'order','suborder','infraorder','section','subsection', 'superfamily','family',
                               'subfamily','tribe','subtribe','genus','subgenus','speciesGroup','species','subspecies') %>% rev()
  
  #### Note
  # This code is nor really used but could be useful. It gets a taxon table from the API and do save it in the main list.
  
    full_taxon_table <- neonUtilities::getTaxonTable('MACROINVERTEBRATE')
  
  #### Save New Data
    inv_allTabs <- do.call(c, 
                           list(inv_allTabs, # NEON download
                                list(full_taxon_table = full_taxon_table), # NEON dl
                                list(date_of_download = Sys.time()))) # unique id
    macroinverts_raw <- inv_allTabs
    #use_data(macroinverts_raw, overwrite = TRUE)
  
  ### Load Data
  # If you did not download a new version of the data, then here, load the data. Note the the raw csv tables are in `data-raw`
  
    load_all()
    inv_allTabs <- macroinverts_raw
    data.frame(tables = labels(inv_allTabs))
  
    ### Check Data
    # Check to see if all match
    ifelse(
      length(
        setdiff(
          inv_allTabs$inv_taxonomyProcessed$sampleID,
          inv_allTabs$inv_fieldData$sampleID)
      )>0,
      "STOP There is an error",
      "All Good!")
  
    ### Join Data
    # Join the cleaned taxonomy (processed taxonomy) with the field data (sampling data). 
    # Make a data set with a density variable `den`. This `den` variable is what is being tested for variance and stability in the NEON temporal analyses. 
    
    #merge/join tables: processing macroinverts data to get to density/abundance 
    # join cleaned taxonomy to sample data
    inv_dat <- left_join(inv_allTabs$inv_taxonomyProcessed, 
                         inv_allTabs$inv_fieldData, 
                         by = c('sampleID')) %>% 
      mutate(den = estimatedTotalCount/benthicArea) %>% # make density
      mutate(scientificName = fct_explicit_na(scientificName)) %>% # explicit missing
      dplyr::filter(sampleCondition == "condition OK") # toss samples low quality
    
    ### Tidy Column Names
    # remove duplicate col names and .x suffix
    inv_dat <- inv_dat[,!grepl('\\.y',names(inv_dat))]
    names(inv_dat) <- gsub('\\.x','',names(inv_dat))

    ### Choose Taxonomic Scale
    # Toss all individuals not identified to the genus or lower taxonomic resolution (e.g., all individuals id-ed as Chironomidae are tossed).

    # get genus and finer resolution using ordered taxon_rank_list
    inv_dat$taxonRank_ordered <- factor(
      inv_dat$taxonRank,
      levels = taxon_rank_list_ordered,
      ordered = TRUE) 
    # get all records that have rank <= genus, where genus is not NA or blank
    inv_dat_fine <- inv_dat %>%
      filter(taxonRank_ordered <= 'genus') %>% # <= due to ordered factor
      filter(!is.na(genus), genus != '') # there are missing genera so toss them
    
    
    # remove 2020 because of incomplete sampling (covid-19) and Alaska (domains 18 & 19), Puerto Rico (D04), and Hawaii (D20)
    inv_dat_fine <- inv_dat_fine %>%
      mutate(year = substr(collectDate, 1, 4))
    inv_dat_fine <- dplyr::filter(inv_dat_fine, year != 2020)
    inv_dat_fine <- dplyr::filter(inv_dat_fine, domainID != "D18")
    inv_dat_fine <- dplyr::filter(inv_dat_fine, domainID != "D19")
    inv_dat_fine <- dplyr::filter(inv_dat_fine, domainID != "D04")
    inv_dat_fine <- dplyr::filter(inv_dat_fine, domainID != "D20")
    inv_dat_fine <- inv_dat_fine %>%
      mutate(month = substr(collectDate, 6,7))
    # drop the lake fish data and keep stream data only
    inv_dat_fine <- inv_dat_fine %>% dplyr::filter(aquaticSiteType == "stream")
    inv_dat_fine <- dplyr::filter(inv_dat_fine, !grepl(pattern = "/", genus))
    
    ### Select sites that meet our criteria (at least 3 years of data with at least 2 replicate-months within each year)
    dat1 <- inv_dat_fine %>%
      group_by(siteID, year, month) %>%
      select(siteID, year, month) %>%
      summarise(count=n()) %>%
      arrange(siteID, year, month) %>%
      ungroup
    
    dat2 <- dat1 %>%
      group_by(siteID, year) %>%
      select(siteID, year) %>%
      summarise(count=n()) %>%
      ungroup()
    
    # select sites with >=3 years of data and >=2 replicate-months (bouts) within each year
    dat3 <- dat2 %>%
      filter(count >= 2) %>%
      ungroup()
    
    dat4 <- 
      dat3 %>%
      group_by(siteID) %>%
      select(siteID) %>%
      summarise(count=n())
    
    dat5 <-
      dat4 %>%
      filter(count >= 3) #these are the final sites
    
    inv_dat_fine <- inv_dat_fine %>%
      dplyr::filter(siteID %in% dat5$siteID)
    
    # concatenate siteID with year to make a separate column for combinations of site+year that might not meet criteria
    # because within final sites there still might be years that do not meet criteria of >1 bouts
    inv_dat_fine <- inv_dat_fine %>% dplyr::mutate(siteyear = paste(siteID, year, sep = "_"))
    inv_dat_fine <- inv_dat_fine %>% dplyr::mutate(bout = paste(year, month, sep = "_"))
    
    dat6 <- inv_dat_fine %>%
      group_by(siteID, year, siteyear, month) %>%
      select(siteID, year, siteyear, month) %>%
      summarise(count=n()) %>%
      arrange(siteID, year, month) %>%
      ungroup
    
    dat7 <- dat6 %>%
      group_by(siteID, year, siteyear) %>%
      select(siteID, year, siteyear) %>%
      summarise(count=n()) %>%
      ungroup()
    
    dat7 <- dat7 %>%
      filter(count >= 2) %>%
      ungroup() #this removes one site-year combination that has only 1 bout
    
    inv_dat_fine <- inv_dat_fine %>%
      dplyr::filter(siteyear %in% dat7$siteyear)
    
    ### Choose Temporal Scale
    # Analyses are run at a specific unit of rime and that scale is the year and month (bout).
    
    ### Get mean latitude and longitude per neon site (certain measurements within a site have slightly different lat and lon, but our analysis is on the site level)
    mean.coord <- inv_dat_fine %>%
      group_by(siteID) %>%
      summarize(lat = mean(decimalLatitude, na.rm=TRUE), lon = mean(decimalLongitude, na.rm=TRUE)) 
    data_macroinvertebrate <- left_join(inv_dat_fine, mean.coord, by="siteID") 
    data_macroinvertebrate <- data_macroinvertebrate %>%
      mutate(scientificName = genus) #to get the same nomenclature as in other groups
    
    return(data_macroinvertebrate)
}

data_macroinvertebrate = map_neon_data_to_ecocomDP.MACROINVERTEBRATE(neon.data.product.id = "DP1.20120.001")
#save processed data
saveRDS(data_macroinvertebrate, file = "data/macroinv_processed.rds")
#data_macroinvertebrate <- readRDS(file = "data/macroinv_processed.rds")



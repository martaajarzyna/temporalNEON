##############################################################################################
#' @title NEON fish data cleaning for biodiversity 

#' @author
#' Stephanie Parker \email{sparker@battelleecology.org} \cr

#' @description Pull, join tables, and clean fish data for NEON Summit 2019

#' @return Cleaned dataframe with joined tables for fish data

# changelog and author contributions / copyrights
#   Stephanie Parker (2019-10-16)
#   Using base code form Matt in NEON Summit 2019 Working Group 14
##############################################################################################

rm(list=ls())
options(stringsAsFactors=F, strip.white=T)

library(dplyr)
library(neonUtilities)
library(lubridate)

#download NEON data using API
fsh_allTabs <- loadByProduct(dpID = "DP1.20107.001", 
  site = "all", package = "expanded", 
  check.size = FALSE)
list2env(fsh_allTabs,envir=.GlobalEnv)


########
names(fsh_allTabs)

field_dat <- dplyr::select(fsh_allTabs$fsh_fieldData, reachID, 
    measuredReachLength, fixedRandomReach, samplingImpractical)

#combine field and pass tables
fsh_dat <- left_join(fsh_allTabs$fsh_perPass, field_dat, by = "reachID") %>% 
    filter(is.na(samplingImpractical)) #remove records where fish couldn't be collected

#individual fish
fsh_dat_indiv <- left_join(fsh_allTabs$fsh_perFish, fsh_dat, by = "eventID") 

fsh_dat_indivSum <- fsh_dat_indiv %>%
    group_by(siteID.x, passStartTime.x, eventID, reachID, scientificName) %>% 
    count(scientificName) %>%
    mutate(year = year(passStartTime.x))
names(fsh_dat_indivSum) <- gsub("\\.x", "", names(fsh_dat_indivSum))

#fill in missing reachID
fsh_dat_indivSum$reachID <- ifelse(is.na(fsh_dat_indivSum$reachID), substr(fsh_dat_indivSum$eventID, 1, 16), fsh_dat_indivSum$reachID) 

#bulk counts
fsh_dat_bulk <- left_join(fsh_allTabs$fsh_bulkCount, fsh_perPass, by = "eventID")

fsh_dat_bulkSum <- fsh_dat_bulk %>%
  group_by(siteID.x, passStartTime.x, eventID, reachID, scientificName) %>% 
  count(scientificName) %>%
  mutate(year = year(passStartTime.x))
names(fsh_dat_bulkSum) <- gsub("\\.x", "", names(fsh_dat_bulkSum)) 


#fill in missing reachID
fsh_dat_bulkSum$reachID <- ifelse(is.na(fsh_dat_bulkSum$reachID), substr(fsh_dat_bulkSum$eventID, 1, 16), fsh_dat_bulkSum$reachID) 

#combine indiv and bulk data
fsh_dat_indivSum$siteID.x <- NA
fsh_dat_bulkSum$siteID.x <- NA
fsh_dat_indivSum$passStartTime.x <- NA
fsh_dat_bulkSum$passStartTime.x <- NA

fsh_dat_all <- bind_rows(fsh_dat_indivSum, fsh_dat_bulkSum) #'.x' colNames were causing errors if not included

fsh_dat_allSum <- fsh_dat_all %>%
  group_by(siteID, reachID, scientificName) %>% 
  count(scientificName)

fsh_dat_allSum$collectDate <- substr(fsh_dat_allSum$reachID, 6, 13)

fsh_dat_allSum <- fsh_dat_allSum[order(fsh_dat_allSum$reachID, 
    fsh_dat_allSum$scientificName),]




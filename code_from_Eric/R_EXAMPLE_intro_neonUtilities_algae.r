# contact: Eric Sokol (esokol@battelleecology.org)

# the NEON package from CRAN
library(neonUtilities)

# the algae data product ID:
# algae DP1.20166.001

# # all data -- could be big
# alg_allTabs <- loadByProduct(dpID = "DP1.20166.001", 
#                              site = "all", package = "expanded", 
#                              check.size = TRUE)

# Try this:
# only 2 sites, restricted in time
alg_allTabs <- loadByProduct(dpID = "DP1.20166.001", 
                             site = c("MAYF", "PRIN"),
                             startdate = "2016-1", enddate = "2018-11",  
                             package = "expanded", check.size = TRUE)

# pull out the field data into a data.frame
alg_field_data <- alg_allTabs$alg_fieldData

# want to use acceptedTaxonID
# consider taxon resolution
# consider algalParameterUnit
# algalParameterValue is the "count" or "density
# upll out the taxonomy data
alg_tax_long <- alg_allTabs$alg_taxonomyProcessed

# or pull by data table...
alg_fieldData_v2 <- neonUtilities::getDatatable(
  dpid = 'DP1.20166.001',
  data_table_name = 'alg_fieldData',
  sample_location_list = c("MAYF", "PRIN"))

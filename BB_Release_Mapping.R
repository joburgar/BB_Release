# Copyright 2019 Province of British Columbia
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at 
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#####################################################################
R_version <- paste0("R-",version$major,".",version$minor)
.libPaths(paste0("C:/Program Files/R/",R_version,"/library")) # to ensure reading/writing libraries from C drive
#####################################################################
## FLEX – Standalone runner (no SpaDES)
#####################################################################

## ---- Libraries ----------------------------------------------------
library(tidyverse)

# Load Packages
list.of.packages <- c("tidyverse","bcdata", "bcmaps","sp","sf", "Cairo", "terra")
# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)


###--- function to retrieve geodata from BCGW

retrieve_geodata_aoi <- function(ID, aoi) {
  
  aoi <- st_make_valid(aoi)
  
  aoi.geodata <- bcdc_get_data(ID) %>% 
    st_as_sf()
  
  aoi.geodata <- st_transform(aoi.geodata, st_crs(aoi))
  aoi.geodata <- st_intersection(aoi.geodata, aoi)
  
  aoi.geodata$Area_km2 <- units::drop_units(st_area(aoi.geodata)) * 1e-6
  
  return(aoi.geodata)
}

###--- function to retrieve raster data

clip_raster_to_aoi <- function(r, aoi) {
  
  r <- rast(r)
  
  # Convert AOI to SpatVector if it's sf
  if (inherits(aoi, "sf") | inherits(aoi, "sfc")) {
    aoi <- vect(aoi)
  }
  
  aoi <- project(aoi, crs(r))
  r_clip <- crop(r, aoi) |> mask(aoi)
  
  return(r_clip)
}

#####################################################################################

# if you have location data - first find the WMU based on the spatial coordinates
# aoi <- SP_cam[[1]] %>% st_transform(crs = 3005)
# aoi_utm <- SP_cam[[2]] %>% st_intersection(aoi %>% st_transform(26910))
# sa_points <- SP_cam[[3]]

# if you only have WMU data
bcdc_search("WHSE_WILDLIFE_MANAGEMENT", res_format = "wms")
bb_WMUs <- bcdc_query_geodata("028d4791-1241-437a-9f7b-fdf08b0d6dfb") %>% 
  filter(WILDLIFE_MGMT_UNIT_ID %in% c("2-3","2-4","2-5","2-8")) %>%
  collect()

ggplot()+
  geom_sf(data=bb_WMUs)

# load covariates from bcdata
# using the bc data warehouse option to clip to aoi
aoi <- bb_WMUs %>%
  st_make_valid() %>%
  st_transform(3005) %>%     # or match dataset CRS
  st_union()                 # optional but helpful


# watercourses layer (to ensure water source for bears)
# bcdc_search("NTS BC River", res_format = "wms")
aoi.RLW <- retrieve_geodata_aoi(aoi=aoi, ID = "414be2d6-f4d9-4f32-b960-caa074c6d36b") 

# transportation layer (Digital Road Atlas) (to consider access)
# bcdc_search("road", res_format = "wms")
aoi.DRA <- retrieve_geodata_aoi(aoi=aoi, ID = "bb060417-b6e6-4548-b837-f9060d94743e")

# parcelmap (to exclude bears from private land and be far from people)
# bcdc_search("parcelmap", res_format = "wms")
aoi.pcl <- retrieve_geodata_aoi(aoi=aoi, ID = "4cf233c2-f020-4f7a-9b87-1923252fbc24")

# protected areas (to exclude bears)
# bcdc_search("WHSE_TANTALIS.TA_PARK_ECORES_PA", res_format = "wms")
aoi.park <- retrieve_geodata_aoi(aoi=aoi, ID = "1130248f-f1a3-4956-8b2e-38d29d3e4af7")



# Load raster (food sources)
r.huck <- clip_raster_to_aoi(aoi=aoi, r="HuckProbability_2025_R2.tif")
plot(r.huck)
r.gbspring <- clip_raster_to_aoi(aoi=aoi, r="hsf_16quant_burgar_AOI_spring.tif")
plot(r.gbspring)

save.image("BB_Release_MapFiles.RData")
# load("BB_Release_MapFiles.RData")


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

## ---- Libraries ----------------------------------------------------
# Load Packages
list.of.packages <- c("tidyverse","bcdata", "sf", "terra")
# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)


###--- function to retrieve bcdata to aoi

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

###--- function to identify area per WMU

run_wmu_analysis <- function(wmu_ids,
                             res = 100,
                             water_dist = 100,
                             fsr_min = 2000,
                             fsr_max = 8000,
                             density_radius = 1000) {
  
  ## --- WMU ---
  # wmu_ids <- c("2-6")
  
  bb_WMUs <- bcdc_query_geodata("028d4791-1241-437a-9f7b-fdf08b0d6dfb") %>%
    filter(WILDLIFE_MGMT_UNIT_ID %in% wmu_ids) %>%
    collect() %>%
    st_make_valid() %>%
    st_transform(3005)
  
  aoi <- st_union(bb_WMUs)

  
  ## --- Raster template ---
  r_template <- rast(
    ext(vect(aoi)),
    resolution = as.numeric(res)[1],crs = "EPSG:3005")
  
  
  ## --- WATER
  water <- retrieve_geodata_aoi("414be2d6-f4d9-4f32-b960-caa074c6d36b", aoi) %>%
    filter(DESCRIPTION %in% c("Lake","River","Wetland"))
  
  ## --- Step 2: reduce AOI for heavy layers ---
  aoi_small <- water %>%
    st_buffer(5000) %>%   # adjust (3–10 km works well)
    st_union()
  
  ## --- ROADS
  roads <- st_read("DRA_MPAR_line.shp") %>% st_as_sf() %>%
    st_transform(st_crs(aoi_small)) %>%
    st_intersection(aoi_small)
  
  # roads %>% group_by(FTYPE) %>% count(ROAD_CLASS) %>% print(n=25)
  roads <- roads %>% filter(FTYPE=="Road") %>% 
    filter(ROAD_CLASS %in% c("recreation","resource","restricted","highway")) %>%
    select(ROAD_CLASS, FTYPE, RDNAME, RDSURFACE) 
  
  ## --- FSR extraction ---
  fsr <- roads %>% filter(ROAD_CLASS=="resource")
  
  ## --- LAND OWNERSHIP
  parcels <- st_read("PMBC_PF_O_polygon.shp") %>% st_as_sf() %>%
    st_transform(st_crs(aoi_small)) %>%
    st_intersection(aoi_small) %>%
    st_simplify(dTolerance = 20)
  parcels %>% group_by(OWNER_TYPE) %>% count(CLASS) %>% print(n=34)
  notcrown <- parcels %>% filter(OWNER_TYPE!="Crown Provincial")
  crown <- parcels %>% filter(OWNER_TYPE=="Crown Provincial")
  
  ## --- PARKS
  parks <- retrieve_geodata_aoi("1130248f-f1a3-4956-8b2e-38d29d3e4af7", aoi)
  
  ## --- CUTBLOCKS
  cutblocks <- st_read("CNS_CUT_BL_polygon.shp") %>% st_as_sf() %>%
    st_transform(st_crs(aoi)) %>%
    st_intersection(aoi)
  
  # Keep early seral (e.g., <= 20 yrs since harvest)
  cutblocks <- cutblocks %>%
    filter(!is.na(HRVST_MD_R)) %>%
    mutate(age = 2026 -HRVST_MD_R) %>%
    filter(age <= 20)

  
  ## --- Rasterize ---
  r_water  <- rasterize(vect(water),   r_template, field = 1, background = 0)
  r_roads  <- rasterize(vect(roads),   r_template, field = 1, background = 0)
  r_fsr    <- rasterize(vect(fsr),     r_template, field = 1, background = 0)
  r_crown <- rasterize(vect(crown), r_template, field = 1, background = 0)
  r_notcrown <- rasterize(vect(notcrown), r_template, field = 1, background = 0)
  r_parks  <- rasterize(vect(parks),   r_template, field = 1, background = 0)
  r_cut    <- rasterize(vect(cutblocks), r_template, field = 1, background = 0)
  
  ## --- Distance ---
  d_water  <- distance(r_water)
  d_fsr    <- distance(r_fsr)
  d_parcel <- distance(r_notcrown)
  
  ## --- ROAD DENSITY ---
  w <- focalMat(r_roads, d = density_radius, type = "circle")
  road_density <- focal(r_roads, w = w, fun = "sum", na.rm = TRUE)
  
  max_val <- as.numeric(global(road_density, "max", na.rm = TRUE)[1,1])
  road_density <- road_density / max_val
  road_density_suit <- 1 - road_density
  
  ## --- CRITERIA ---
  
  # 1 Water
  water_bin <- d_water <= water_dist
  
  # 2 Parks (mask)
  park_mask <- classify(r_parks, cbind(1, NA))
  
  # 3 Human
  human_pressure <- exp(-d_parcel / 1000)
  human_suit <- (1 - human_pressure) * road_density_suit
  
  # 4 FSR depth
  fsr_suit <- ifel(d_fsr >= fsr_min & d_fsr <= fsr_max, 1, 0)
  
  ## --- HABITAT ---
  
  # Huck
  r.huck <- resample(
    clip_raster_to_aoi("HuckProbability_2025_R2.tif", aoi),
    r_template
  )
  
  norm <- function(x) (x - minmax(x)[1]) /
    (minmax(x)[2] - minmax(x)[1])
  
  huck_suit <- norm(r.huck)
  
  # Early seral = cutblocks (binary → light smoothing helps)
  cut_suit <- focal(r_cut, w = focalMat(r_cut, 300, "circle"), fun = mean, na.rm = TRUE)
  
  # Combine habitat
  habitat <- (0.6 * huck_suit + 0.4 * cut_suit)
  
  ## --- FINAL MODEL ---
  suitability <- (0.25 * water_bin +
                    0.35 * habitat +
                    0.20 * human_suit +
                    0.20 * fsr_suit)
  
  suitability <- mask(suitability, park_mask)
  names(suitability) <- "BB_suitability"
  
  ## --- TOP AREAS ---
  vals <- values(suitability)
  thresh <- quantile(vals, 0.8, na.rm = TRUE)
  
  best <- suitability >= thresh
  best_poly <- as.polygons(best, dissolve = TRUE) |> st_as_sf()
  ggplot() + geom_sf(data=best_poly)
  ggplot() + geom_sf(data=aoi)
  
  return(list(
    suitability = suitability,
    best_areas = best_poly,
    cutblocks = r_cut,
    road_density = road_density
  ))
}
#####################################################################################


result <- run_wmu_analysis(c("2-3","2-4","2-5","2-6","2-8"))
st_write(result$best_areas, "BB_BestAreas_2026.shp")
writeRaster(result$road_density, "BB_RoadDensity_2026.tif", overwrite = TRUE)
writeRaster(result$cutblocks, "BB_Cutblocks_2026.tif", overwrite = TRUE)
writeRaster(result$suitability, "BB_suitability_2026.tif", overwrite = TRUE)


library(terra)
library(tidyverse)

prev <- rast("Prevalence-maps/Prevalence-median-GeoStat.tif")
suit <- rast("Bd-Suitability/DNC-175.tif")

crs(prev) <- crs("EPSG:24892")

spp <- lapply(list.files("Anuran-species/By-status/", ".shp", full.names = T), vect)

spp <- lapply(spp, function(x){project(x, crs(prev))})

r0 <- prev-prev

spp.r <- lapply(spp, function(x){
    x1 <- rasterize(x, prev, fun = "sum")
    x2 <- merge(x1, r0)
    return(x2)})

spp.r <- c(spp.r[[1]], spp.r[[2]], spp.r[[3]],
           spp.r[[4]], spp.r[[5]], spp.r[[6]],
           spp.r[[7]])

names(spp.r) <- c("CE", "DD", "En", "LC", "NI", "NT", "Vu")

writeRaster(spp.r, "Anuran-species/By-status/Richness-status.tif")

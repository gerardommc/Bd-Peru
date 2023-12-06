library(terra)
library(ntbox)
library(foreach)
library(tidyr)

# Occurrence data
bd.per <- read.csv("BD-occurrence/Peru/Bd-prevalence-data.csv")

bd.glob <- read.csv("BD-occurrence/Global/Bd-global-data.csv")

chels.pca <- rast(list.files("../Chelsa/PCA", "tif", full.names = T))

peru <- vect("Admin-layers/gadm41_PER.gpkg")

ref.r <- rast("Prevalence-maps/Prevalence-median-GeoStat.tif")
crs(ref.r) <- crs("EPSG:24892")

chels.peru <- crop(chels.pca, peru)

chels.peru <- project(chels.peru, crs(ref.r))
chels.peru <- resample(chels.peru, ref.r)
chels.peru <- mask(chels.peru, ref.r)

# Fitting ellipsoid
bd.glob.pres <- subset(bd.glob, SppDet > 0, Country != "Peru")

env.bd.glob <- terra::extract(chels.pca, bd.glob.pres[, c("longitude", "latitude")]) |> as.data.frame() |> na.omit()
env.bd.glob <- env.bd.glob[, -1]

combs <- combn(1:10, m = 4)

centres <- foreach(i = 1:ncol(combs)) %do% {
    cov_center(env.bd.glob, mve = T, vars = names(env.bd.glob)[combs[, i]], level = 0.95)
    }

dir.create("Bd-Suitability")
saveRDS(centres, "Bd-Suitability/Centres-covariances.rds")

library(raster)

chels.peru1 <- stack(chels.peru)

ellip <- foreach(j = seq_along(centres), .combine = c) %do% {
    lay <- dropLayer(chels.peru1, 
                             i = which(!names(chels.peru1) %in% names(env.bd.glob)[combs[, j]]))
    fit <- ellipsoidfit(envlayers = lay, 
                      centroid = centres[[j]]$centroid,
                      covar = centres[[j]]$covariance,
                      level = 0.95,
                      size = 1, plot = F)
    return(fit$suitRaster)
}

bd.pres <- subset(bd.per, N.positive > 0)
coordinates(bd.pres) <- ~ Longitude + Latitude
proj4string(bd.pres) <- CRS("+init=epsg:4326")
bd.pres <- spTransform(bd.pres, CRSobj = CRS(proj4string(chels.peru1)))

library(doParallel)
registerDoParallel(cores = 24)

roc.tests <- foreach(i = seq_along(ellip)) %dopar% {
    pROC(continuous_mod = ellip[[i]],
         test_data = coordinates(bd.pres),
         n_iter = 10000,
         E_percent = 5,
         boost_percent = 50,
         parallel = F,
         rseed = 63193)
}

roc.results <- foreach(i = seq_along(roc.tests), .combine = rbind) %do%{
    roc.tests[[i]]$pROC_summary
}
roc.results <-  data.frame(roc.results)
roc.results$Model <- paste0("DNC-", seq_along(roc.tests))

best.10 <- sort(roc.results$Mean_pAUC_ratio_at_5., decreasing = T)[1:10]

best.suits <- ellip[which(roc.results$Mean_pAUC_ratio_at_5. %in% best.10)]

write.csv(roc.results, "Bd-Suitability/Partial-roc-results.csv", row.names = F)
write.csv(roc.results[which(roc.results$Mean_pAUC_ratio_at_5. %in% best.10),], 
          "Bd-Suitability/Partial-roc-Best.csv", row.names = F)

for(i in 1:10){writeRaster(best.suits[[i]], 
                          paste0("Bd-Suitability/DNC-",
                                 which(roc.results$Mean_pAUC_ratio_at_5. %in% best.10)[i],
                                 ".csv"), "GTiff", overwrite = T)}

# Weighted mean model

best.suits.s <- stack(best.suits)

## Harmonic weighted mean
w.model <- 1/(weighted.mean(1/best.suits.s, exp(best.10)))

writeRaster(w.model, "Bd-Suitability/DNC-weighted-mean", "GTiff", overwrite = T)

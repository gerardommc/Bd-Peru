# Chelsa PCA
library(terra)
library(tidyr)

chelsa <- rast(list.files("../Chelsa/", "tif", full.names = T))
dir.create("../Chelsa/PCA")

rp <- data.frame(x = runif(50000, -180.0, 180.0), y = runif(50000, -90.0, 90.0))

library(doParallel)
registerDoParallel(cores = 4)

chelsa.ag <- aggregate(chelsa, 5)

samples <- foreach(i = 1:19, .combine = cbind) %dopar% {
    terra::extract(chelsa.ag[[i]], rp)[, -1] |> data.frame()
}

samples <- samples |> as.data.frame() |> na.omit()

names(samples) <- names(chelsa)

chels.m <- colMeans(samples)
chels.sd <- apply(samples, 2, sd)

samples.cent <- sweep(samples, 2, chels.m, FUN = "-")
samples.scale <- sweep(samples.cent, 2, chels.sd, FUN = "/")

chels.scale <- (chelsa.ag - chels.m)/chels.sd

pca <- princomp(samples.scale)

saveRDS(pca, "Environment-raster/Chelsa-PCA-results.rds")

rm(chelsa); gc(reset = T)

chels.pca <- predict(chels.scale, pca, index = 1:10)

for(i in 1:10){writeRaster(chels.pca[[i]], paste0("../Chelsa/PCA/Chelsa-PCA-", i, ".tif"))}

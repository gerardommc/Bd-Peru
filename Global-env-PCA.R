# Chelsa PCA
library(raster)

chelsa <- stack(list.files("../../Chelsa/", "tif", full.names = T))
dir.create("../../Chelsa/PCA")

rp <- dismo::randomPoints(chelsa, 10000)
samples <- na.omit(data.frame(extract(chelsa, rp)))
pca <- princomp(samples)

saveRDS(pca, "Environment-raster/Chela-PCA-resulds.rds")

chels.pca <- predict(chelsa, pca, index = 1:10)

for(i in 1:10)writeRaster(round(chels.pca[[i]], 2), paste0("../../Chelsa/PCA/Chelsa-PCA-", i), 
                         "GTiff", overwrite = T)
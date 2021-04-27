library(raster)
library(ntbox)

# Occurrence data
bd.per <- read.csv("BD-occurrence/Peru/Bd-prevalence-data.csv")

bd.glob <- read.csv("BD-occurrence/Global/Bd-global-data.csv")

# Environmental data
chelsa <- stack(list.files("../../Niche centroids/Chelsa/", "tif", full.names = T))

rp <- dismo::randomPoints(chelsa, 10000)
samples <- data.frame(extract(chelsa, rp))
pca <- princomp(samples)

chels.pca <- predict(chelsa, pca, index = 1:3)

# Fitting ellipsoid
bd.glob.pres <- subset(bd.glob, SiteBdDetected == "D")

env.bd.glob <- na.omit(data.frame(extract(chels.pca, bd.glob.pres[, c("longitude", "latitude")])))

centre <- cov_center(env.bd.glob, mve = T, vars = paste0("layer." ,c(1, 2, 3)), level = 0.95)

ellip <- ellipsoidfit(envlayers = chels.pca, 
                      centroid = centre$centroid,
                      covar = centre$covariance,
                      level = 0.99,
                      size = 1, plot = T)

ellip

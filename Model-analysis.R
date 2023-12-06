library(terra)

prev <- rast("Prevalence-maps/Prevalence-median-GeoStat.tif")
suit <- rast("Bd-Suitability/DNC-175.tif")

div <- rast("Anuran-species/Peru-anuran-species-PSAD56.tif")

crs(prev) <- crs("EPSG:24892")

suit <- project(suit, crs(prev))
div <- project(div, crs(prev))

suit <- resample(suit, prev)
div <- resample(div, prev)

model.data <- c(prev, suit, div)

model.data.df <- model.data |> as.data.frame() |> na.omit()
names(model.data.df) <- c("Prevalence", "Suitability", "Richness")

library(ggplot2)

library(viridis)

ggplot(model.data.df) + geom_hex(aes(x = Prevalence, y = Suitability, fill = after_stat(log10(count))), stat = "binhex") +
    geom_smooth(aes(x = Prevalence, y = Suitability), colour = "red")
ggplot(model.data.df) + geom_hex(aes(x = Prevalence, y = Richness, fill = after_stat(log10(count))), stat = "binhex")+
    geom_smooth(aes(x = Prevalence, y = Richness), colour = "red")
ggplot(model.data.df) + geom_hex(aes(x = Richness, y = Suitability, fill = after_stat(log10(count))), stat = "binhex")+
    geom_smooth(aes(x = Richness, y = Suitability), colour = "red")

# NPAs

npas <- vect("Protected-areas/Peru-NPAs.gpkg")

npas <- project(npas, crs(prev))

model.npas <- extract(model.data, npas, fun = "mean", na.rm = T)

npas <- cbind(npas, model.npas)

writeVector(npas,"Protected-areas/Peru-NPAs-BD-data.gpkg")

# Anuran species

spp <- vect("Anuran-species/Peru-Anurans.shp")
spp <- project(spp, crs("EPSG:24892"))

model.spp.mean <- extract(model.data, spp, fun = "mean", na.rm = T)[, -1]
model.spp.max <- extract(model.data, spp, fun = "max", na.rm = T)[, -1]
model.spp.sd <- extract(model.data, spp, fun = "sd", na.rm = T)[, -1]

names(model.spp.mean) <- paste0(names(model.spp.mean), ".mean")
names(model.spp.max) <- paste0(names(model.spp.max), ".max")
names(model.spp.sd) <- paste0(names(model.spp.sd), ".sd")

spp.model <- cbind(spp, model.spp.mean,
                   model.spp.max, model.spp.sd)

names(spp.model)

spp.mod.df <- spp.model |> as.data.frame()

cons.stat <-read.csv("Anuran-species/Conservation-status.csv")

spp.mod.df$Status <- NA

for(i in 1:nrow(spp.mod.df)){
    id <- which(cons.stat$scientificName == spp.mod.df$BINOMIAL[i])
    if(length(id) == 0){
        next
    } else {
        spp.mod.df$Status[i] <- cons.stat$redlistCategory[id]
    }
}

library(ggplot2)

ggplot(spp.mod.df) + geom_boxplot(aes(x = Status, y = prevalence.mean))
ggplot(spp.mod.df) + geom_boxplot(aes(x = Status, y = `DNC-175.mean`))



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

ggplot(model.data.df) + geom_hex(aes(x = Suitability, y = Prevalence, fill = after_stat(log10(count))), stat = "binhex") +
    geom_smooth(aes(x = Suitability, y = Prevalence), colour = "red") +
    scale_fill_gradientn(colours = viridis(100)) +
    labs(fill = expression(log[10] (count)))
ggplot(model.data.df) + geom_hex(aes(x = Prevalence, y = Richness, fill = after_stat(log10(count))), stat = "binhex")+
    geom_smooth(aes(x = Prevalence, y = Richness), colour = "red")+
    scale_fill_gradientn(colours = viridis(100))+
    labs(fill = expression(log[10] (count)))
ggplot(model.data.df) + geom_hex(aes(x = Suitability, y = Richness, fill = after_stat(log10(count))), stat = "binhex")+
    geom_smooth(aes(x = Suitability, y = Richness), colour = "red")+
    scale_fill_gradientn(colours = viridis(100))+
    labs(fill = expression(log[10] (count)))

pdf("Correlations.pdf", width = 6, height = 5)
ggplot(model.data.df) + geom_hex(aes(x = Suitability, y = Prevalence, fill = after_stat(log10(count))), stat = "binhex") +
    geom_smooth(aes(x = Suitability, y = Prevalence), colour = "red") +
    scale_fill_gradientn(colours = viridis(100)) +
    labs(fill = expression(log[10] (count)))
ggplot(model.data.df) + geom_hex(aes(x = Prevalence, y = Richness, fill = after_stat(log10(count))), stat = "binhex")+
    geom_smooth(aes(x = Prevalence, y = Richness), colour = "red")+
    scale_fill_gradientn(colours = viridis(100))+
    labs(fill = expression(log[10] (count)))
ggplot(model.data.df) + geom_hex(aes(x = Suitability, y = Richness, fill = after_stat(log10(count))), stat = "binhex")+
    geom_smooth(aes(x = Suitability, y = Richness), colour = "red")+
    scale_fill_gradientn(colours = viridis(100))+
    labs(fill = expression(log[10] (count)))
dev.off()

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

library(sf); library(tidyverse)

spp2 <- st_read("Anuran-species/Peru-Anurans.shp")

areas <- st_area(spp2)

spp.model <- cbind(spp, model.spp.mean,
                   model.spp.max, model.spp.sd)
spp.model$Area <- as.numeric(areas/1e+6)

cents <- centroids(spp, inside = T)
cents$Area <- as.numeric(areas/1e+6)

writeVector(cents, "Anuran-species/Centroid-distributions.shp", overwrite = T)

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

NAs <- is.na(spp.mod.df$Status)
spp.mod.df$Status[NAs] <- "No information"

st <- factor(spp.mod.df$Status)

st <- ordered(st, c("Least Concern", "No information", "Data Deficient", "Near Threatened",
                    "Vulnerable", "Endangered", "Critically Endangered"))
spp.mod.df$Status <- st

spp.mod.df <- na.omit(spp.mod.df)

library(ggplot2)

pdf("Suit-prev-conservation.pdf", width = 6, height = 5)
ggplot(spp.mod.df) + geom_boxplot(aes(x = Status, y = prevalence.mean)) +
    labs(y = "Prevalence", x = "IUCN Category") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
ggplot(spp.mod.df) + geom_boxplot(aes(x = Status, y = log10(Area))) +
    labs(y = expression(log[10](Area)), x = "IUCN Category") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
ggplot(spp.mod.df) + geom_boxplot(aes(x = Status, y = `DNC-175.mean`)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
    labs(y = "Suitability", x = "IUCN Category") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
dev.off()

spp.mod.df$OR <- with(spp.mod.df, (prevalence.mean+0.01)/(1 - (prevalence.mean + 0.01)))

pdf("Prev-Area-Status.pdf", width = 7, height = 5)
ggplot(spp.mod.df) + geom_point(aes(x = log((prevalence.mean+0.01)/(1-(prevalence.mean+0.01))),
                                    y = log10(Area), colour = Status, size = prevalence.sd), 
                                alpha = 0.5) +
    stat_ellipse(aes(x = log(prevalence.mean/(1-prevalence.mean)), y = log10(Area), colour = Status)) +
    labs(x = "log(Odds ratio)", y = expression(log[10](Area)), 
         size = "Prevalence S. D.", colour = "IUCN status") +
    theme_bw()

ggplot(spp.mod.df) + geom_point(aes(x = `DNC-175.mean`,
                                    y = log10(Area), colour = Status, size = prevalence.sd), 
                                alpha = 0.5) +
    stat_ellipse(aes(x = `DNC-175.mean`, y = log10(Area), colour = Status)) +
    labs(x = "Suitability", y = expression(log[10](Area)), 
         size = "Prevalence S. D.", colour = "IUCN status") +
    theme_bw()

ggplot(spp.mod.df) + geom_point(aes(x = `DNC-175.mean`,
                                    y = log(OR), colour = Status, size = prevalence.sd), 
                                alpha = 0.5) +
    stat_ellipse(aes(x = `DNC-175.mean`, y = log(OR), colour = Status)) +
    labs(x = "Suitability", y = "log(Odds ratio)", 
         size = "Prevalence S. D.", colour = "IUCN status") +
    theme_bw()
dev.off()

#Ternary plot

library(ggtern); library(tidyverse)

spp.mod.t <- data.frame(Status = spp.mod.df$Status)
spp.mod.t$Status <- ordered(spp.mod.t$Status, c("Critically Endangered", "Endangered", "Vulnerable", 
                            "Near Threatened", "Data Deficient",       
                            "No information", "Least Concern"))

reg <- function(x){(x - min(x))/(max(x) - min(x))}

spp.mod.t$OR <- log(spp.mod.df$OR) |> reg()
spp.mod.t$Suitability <- log(spp.mod.df$`DNC-175.mean`) |> reg()
spp.mod.t$Area <- log10(spp.mod.df$Area) |> reg()

spp.mod.t$Prev.sd <- spp.mod.df$prevalence.sd

breaks.L <- with(spp.mod.df, seq(min(log(OR)), max(log(OR)), len = 6)) 
breaks.R <- with(spp.mod.df, seq(min(log(Area)), max(log(Area)), len = 6)) 
breaks.T <- with(spp.mod.df, seq(min(log(`DNC-175.mean`)), max(log(`DNC-175.mean`)), len = 6)) 

breaks.L <- round(breaks.L, 1)
breaks.R <- round(breaks.R, 1)
breaks.T <- round(breaks.T, 1)

pdf("Ternary-plot.pdf", width = 7, height = 6)
ggtern(spp.mod.t, aes(x = OR, y = Suitability, z = Area)) +
    geom_point(aes(colour = Status, size = Prev.sd), alpha = 0.4) +
    theme_showarrows() +
    theme_showgrid_minor() +
    labs(x = "log(Odds ratio)",
         y = "log(Suitability)",
         z = "log(Area)") +
    scale_L_continuous(labels = breaks.L) +
    scale_T_continuous(labels = breaks.T) +
    scale_R_continuous(labels = breaks.R)
dev.off()

## Models for the analysis

library(tidymodels)

spp.mod.analysis <- data.frame(Status = spp.mod.df$Status,
                               OR = log(spp.mod.df$OR),
                               Area = log(spp.mod.df$Area),
                               Suitability = log(spp.mod.df$`DNC-175.mean`))
stand <- function(x) {(x - mean(x))/sd(x)}

spp.mod.analysis$OR <- spp.mod.analysis$OR |> stand()
spp.mod.analysis$Suitability <- spp.mod.analysis$Suitability |> stand()
spp.mod.analysis$Area <- spp.mod.analysis$Area |> stand()

mod.1 <- multinom_reg() |> fit(Status ~ 0 + OR + Suitability + Area, spp.mod.analysis)

sink("Results/Multinomial-model.txt", type = "output")
tidy(mod.1, exponentiate = TRUE, conf.int = TRUE) |> 
    mutate_if(is.numeric, round, 4) |> 
    select(-std.error, -statistic)
sink()

sink("Results/Multinomial-model-performance.txt", type = "output")
Predictions <- mod.1 |> 
    augment(new_data = spp.mod.analysis)

conf_mat(Predictions, truth = Status, estimate = .pred_class)

accuracy(Predictions, truth = Status, estimate = .pred_class)

roc_auc(Predictions, truth = Status, `.pred_Least Concern`, `.pred_No information`,
        `.pred_Data Deficient`, `.pred_Near Threatened`,
        .pred_Vulnerable, .pred_Endangered, `.pred_Critically Endangered`)
sink()

pdf("Results/Multinom-performance.pdf", width = 6, height = 5)
roc_curve(Predictions, truth = Status, `.pred_Least Concern`, `.pred_No information`,
          `.pred_Data Deficient`, `.pred_Near Threatened`,
          .pred_Vulnerable, .pred_Endangered, `.pred_Critically Endangered`) |> 
    ggplot(aes(x = 1 - specificity, y = sensitivity, color = .level)) +
    geom_line(size = 1, alpha = 0.7) +
    geom_abline(slope = 1, linetype = "dotted") +
    coord_fixed() +
    labs(color = NULL) +
    theme_light()
dev.off()

##
pdf("Status-cor-prev-suit-area.pdf", width = 7, height = 5)
ggplot(spp.mod.df) + geom_boxplot(aes(x = log(OR), y = Status), alpha = 0.5, outlier.size = 0) +
    labs(x = "log(Odds ratio)", y = "IUCN status", size = "EOO \nsize") +
    geom_point(aes(x = log(OR), y = Status, size = Area), alpha = 0.2) +
    geom_smooth(aes(x = log(OR), y = Status.num), method = "gam", colour = "darkgrey") +
    theme_bw()
ggplot(spp.mod.df) + geom_boxplot(aes(x = log(`DNC-175.mean`), y = Status), alpha = 0.5, outlier.size = 0) +
    labs(x = "log(Suitability)", y = "IUCN status", size = "EOO \nsize") +
    geom_point(aes(x = log(`DNC-175.mean`), y = Status, size = Area), alpha = 0.2) +
    geom_smooth(aes(x = log(`DNC-175.mean`), y = Status.num), method = "gam", colour = "darkgrey") +
    theme_bw()
ggplot(spp.mod.df) + geom_bar(aes(x = Status), alpha = 0.5) +
    labs(x = "IUCN status", y = "No. of species") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
dev.off()

prev.mod <- readRDS("Model-results/Bd-GeoStat-Prevalence.rds")

library(PrevMap)

sink("Model-results/Prev-model-summary.txt")
summary(prev.mod)
sink()

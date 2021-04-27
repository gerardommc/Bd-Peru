# Rscript to select and fit model based geostatistical model for Bd

# Reading Bd data
bd.dat <- read.csv("BD-occurrence/Peru/Bd-prevalence-data.csv")
bd.dat$logit.prev <- with(bd.dat, log(Prevalence/(1 - Prevalence)))
bd.dat$logit.prev[bd.dat$logit.prev == -Inf] <- -5
bd.dat$logit.prev[bd.dat$logit.prev == Inf] <- 5

# Extracting environmental raster data
library(raster)
lst <- stack(list.files("Environment-raster/LST/", full.names = T))
ndvi <- stack(list.files("Environment-raster/NDVI/", full.names = T))
topo <- stack(list.files("Environment-raster/Topo/", full.names = T))

bd.dat.psad <- bd.dat
coordinates(bd.dat.psad) <- ~ Longitude + Latitude
proj4string(bd.dat.psad) <- CRS("+init=epsg:4326")

bd.dat.psad <- spTransform(bd.dat.psad, CRS(proj4string(lst)))

# Fitting glm model
lst.dat <- data.frame(extract(lst, bd.dat.psad))
ndvi.dat <- data.frame(extract(ndvi, bd.dat.psad))
topo.dat <- data.frame(extract(topo, bd.dat.psad))

lst.reg <- cbind(bd.dat[, c("N.sampled", "N.positive")], lst.dat)
ndvi.reg <- cbind(bd.dat[, c("N.sampled", "N.positive")], ndvi.dat)
topo.reg <- cbind(bd.dat[, c("N.sampled", "N.positive")], topo.dat)

# lst model
m1.lst <- glm(cbind(N.positive, N.sampled - N.positive) ~ LST.PCA.1 + LST.PCA.2 +
              LST.PCA.3+LST.PCA.4+
              I(LST.PCA.1^2) + 
              I(LST.PCA.3^2) + 
              I(LST.PCA.5^2), 
          lst.reg,
          family = binomial)
m2.lst <- step(m1.lst)
summary(m2.lst)
plot(m2.lst)
anova(m2.lst)

# ndvi model

m1.ndvi <- glm(cbind(N.positive, N.sampled - N.positive) ~ NDVI.PCA.1 + NDVI.PCA.2 +
                  NDVI.PCA.3+
                  I(NDVI.PCA.1^2) + 
                  I(NDVI.PCA.2^2) + 
                  I(NDVI.PCA.3^2), 
              ndvi.reg,
              family = binomial)
m2.ndvi <- step(m1.ndvi)
summary(m2.ndvi)
plot(m2.ndvi)
anova(m2.ndvi)

# ndvi and lst

joint.reg <- cbind(bd.dat[, c("N.sampled", "N.positive")], lst.dat, ndvi.dat)

m1.joint <-  glm(cbind(N.positive, N.sampled - N.positive) ~ LST.PCA.1 + 
                     LST.PCA.3 + 
                     LST.PCA.4 +
                     I(LST.PCA.1^2) +
                     I(LST.PCA.5^2),
           joint.reg,
           family = binomial)
m2.joint <- step(m1.joint)
summary(m2.joint)

# ndvi, lst, topo
all.reg <- cbind(bd.dat[, c("N.sampled", "N.positive")], lst.dat, ndvi.dat, topo.dat)

m1.all <- glm(cbind(N.positive, N.sampled - N.positive) ~ LST.PCA.1 + 
                  LST.PCA.3 + 
                  LST.PCA.4 +
                  NDVI.PCA.1 +
                  Topo.PCA.2 + Topo.PCA.3 + Topo.PCA.4 +
                  I(LST.PCA.1^2) +
                  I(LST.PCA.5^2) +
                  I(NDVI.PCA.1^2) +
                  I(Topo.PCA.2^2) + I(Topo.PCA.3^2) + I(Topo.PCA.4^2),
              all.reg,
              family = binomial)
m2.all <- step(m1.all)

summary(m2.all)
anova(m2.all, m2.lst)
plot(m2.all)

anova(m2.all, m2.lst)
AIC(m2.all, m2.lst)

#All interactions with NDVI
m1.ints <- glm(cbind(N.positive, N.sampled - N.positive) ~ NDVI.PCA.1 * (LST.PCA.1 + 
                  LST.PCA.3 + 
                  LST.PCA.4 +
                  Topo.PCA.2 + Topo.PCA.3 + Topo.PCA.4 +
                  I(LST.PCA.1^2) +
                  I(LST.PCA.5^2) +
                  I(Topo.PCA.2^2) + I(Topo.PCA.3^2) + I(Topo.PCA.4^2)) +
                  I(NDVI.PCA.1^2),
              all.reg,
              family = binomial)

m2.ints <- step(m1.ints)
summary(m2.ints)

m3.ints <- update(m2.ints, .~. - LST.PCA.4)
summary(m3.ints)

m4.ints <- update(m3.ints, .~. - LST.PCA.3)
summary(m4.ints)


m5.ints <- update(m4.ints, .~. -I(Topo.PCA.2^2))
summary(m5.ints)

# Fitting variogram
library(geoR)

xy <- coordinates(bd.dat.psad)
vari <- variog(coords=xy,data=bd.dat.psad$logit.prev, angles = F)

plot(vari,xlab="distance (decimal degrees)")

vari.fit <- variofit(vari,ini.cov.pars=c(0.1,0.5),cov.model="matern",
                     fix.nugget=FALSE,nugget=5,
                     fix.kappa=FALSE,kappa=0.1,
                     minimisation.function = "nlm")
lines(vari.fit)

vari.fit

# Fitting geostatistical model
library(PrevMap)

m.fin <- glm(cbind(N.positive, N.sampled - N.positive) ~ NDVI.PCA.1 : LST.PCA.1 + LST.PCA.1 +  
                  I(LST.PCA.1^2) + I(LST.PCA.5^2) + 
                  I(Topo.PCA.3^2) + I(Topo.PCA.4^2), 
             all.reg,
             family = binomial)

params <- c(coef(m.fin),vari.fit$cov.pars,vari.fit$nugget)     #c(beta,sigma2,phi,tau2)

data <- cbind(xy, all.reg)
mcmc.1 <- control.mcmc.MCML(n.sim=50000,burnin=5000,thin=45,h=(1.65)/(nrow(data)^(1/6)))

start.1 <- c(params[length(params)],params[length(params) - 1]/params[length(params)-1])     #starting values of phi and the relative variance of the nugget effect nu2 respectively

# Fitting geo-statistical model
fit.MCML.1 <- binomial.logistic.MCML(formula=N.positive ~  NDVI.PCA.1 : LST.PCA.1 + LST.PCA.1 +  
                                         I(LST.PCA.1^2) + I(LST.PCA.5^2) + 
                                         I(Topo.PCA.3^2) + I(Topo.PCA.4^2),
                                     units.m=~ N.sampled,
                                     par0= params,
                                     coords=~Longitude+ Latitude,
                                     data=data,
                                     control.mcmc=mcmc.1,
                                     kappa=0.5,
                                     start.cov.pars=start.1)

summary(fit.MCML.1)

library(splancs)

predictors.r <- stack(lst, ndvi, topo)
predictors.r <- aggregate(predictors.r, 2)

predictors <- data.frame(rasterToPoints(predictors.r))
names(predictors) <- c("Longitude", "Latitude", names(predictors)[-c(1, 2)])
predictors <- na.omit(predictors)

grid.pred <- predictors[, c("Longitude", "Latitude")]

pred.MCML <- spatial.pred.binomial.MCML(fit.MCML.1,
                                        grid.pred=grid.pred,
                                        predictors=predictors,
                                        control.mcmc=mcmc.1,
                                        type="joint",
                                        scale.predictions="prevalence",
                                        standard.errors=TRUE,
                                        scale.thresholds="prevalence",
                                        thresholds = c(0.025, 0.5, 0.975))

plot(pred.MCML, type = "prevalence")

prevalence <- rasterFromXYZ(data.frame(grid.pred, prevalence = pred.MCML$prevalence$predictions))
exceed <- rasterFromXYZ(data.frame(grid.pred, pred.MCML$exceedance.prob))

proj4string(prevalence) <- CRS(proj4string(lst))
proj4string(exceed) <- CRS(proj4string(lst))

dir.create("Prevalence-maps")

#writing up results
writeRaster(prevalence, "Prevalence-maps/Prevalence-median-m2", "GTiff")
for(i in 1:3)writeRaster(exceed[[i]], paste0("Prevalence-maps/Exceedence-prob-", c("025", "50" , "975")[i], "-m2"), "GTiff")

library(dplyr)

bd.data <- read.csv("../Chytrid data/BD-Presence-absence-Peru.csv")

bd.sum <- bd.data %>% group_by(Latitude, Longitude, Date) %>% summarise(N, Bd.Positive, "sum")

full.bd.dat <- with(bd.sum,
                    data.frame(Longitude = Longitude,
                               Latitude = Latitude,
                               Date = Date,
                               N.sampled = N,
                               N.positive = as.integer(Bd.Positive),
                               Prevalence = Bd.Positive/N))

write.csv(full.bd.dat, "BD-occurrence/Peru/Bd-prevalence-data.csv", row.names = F)

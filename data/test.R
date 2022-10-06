#-*- coding: utf-8 -*-

### File: test.R
### Time-stamp: <2022-10-06 13:26:21 a23579>
###
### Created: 05/10/2022	15:42:50
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################



(test <- range2logNparams(min = 100, max = 300))

(test2 <- range2logNparams(mean = 100, max = 300))

## validation mean:
range2logNparams(min = test2["min"], max = test2["max"])
range2logNparams(mean = test["mean"], max = test["max"])


qlnorm(p = c(0.025, 0.975), meanlog = log(test["mean"]), sdlog = test["logsd"])

qlnorm(p = c(0.025, 0.975), meanlog = log(test2["mean"]), sdlog = test2["logsd"])


range2logNparams(1, 180)

range2logNparams(mean = 100, max = 300)
range2logNparams(min = 20, mean = 100, max = 300)
range2logNparams(min = 20, max = 300)

## ###########################################################################
## Parametrisation depth:

## Haddock:
range2logNparams(mean = 100, max = 300) # NS
range2logNparams(min = 25, max = 300)   # NS alternative

## Saithe:
range2logNparams(mean = 100, max = 300) # NS
range2logNparams(min = 30, max = 300)   # NEA, partially more coastal stock
                                        # (not realistic to use min = 0-10)

## Cod:
range2logNparams(mean = 250, max = 300) # NS... seems very restrictive. May use
range2logNparams(min = 30, max = 300) # NS alternative

## Golden redfish:
range2logNparams(min = 100, max = 400) # NEA

## Atlantic wolffish:
range2logNparams(min = 30, max = 900) # NEA (min 30 instead of zero,
                                      # abundant up to 350, but found up to 900, so
                                      # taking this value as indicative of large
                                        # tolerance)
range2logNparams(min = 20, max = 900)
range2logNparams(min = 30, max = 350, prob = 0.75) # What if we consider preferential rather than limits
range2logNparams(min = 50, max = 350, prob = 0.75) # very sensitive to lower bound.
range2logNparams(min = 50, max = 150, prob = 0.75) # Preferential range for spawning.

## ###########################################################################
## Spatial patchiness parametrisation from the data:

load("Data/[confidential]_wolffish_detailed_data.RData")

dfSBcatchDF0subCW %>%
    head(3) %>%
    as.data.frame()

table(dfSBcatchDF0subCW$gear)

dataSub <- dfSBcatchDF0subCW %>%
    mutate(gear_group = as.factor(trunc(as.numeric(gear)/100))) %>%
    filter(gear_group %in% c(41),
           !gear %in% c(4129,4149)) %>%
    mutate(X = ifelse(is.na(longitudestart), longitudeend, longitudestart),
           Y = ifelse(is.na(latitudestart), latitudeend, latitudestart),
           X = X - min(X, na.rm = TRUE),
           Y = Y - min(Y, na.rm = TRUE))

summary(dataSub$BottomDepth)
summary(dataSub$X)
summary(dataSub$Y)

sum(! is.na(dataSub$BottomDepth))
mean(! is.na(dataSub$BottomDepth))

dataSub2 <- dataSub %>%
    filter(! is.na(X) & ! is.na(Y) & ! is.na(BottomDepth)) %>%
    as.data.frame()

colnames(dataSub2)


map_grid <- expand.grid(X=1:40, Y=1:20)
data.bathym1 <- data.frame(ID = 1:nrow(map_grid), map_grid , BottomDepth=50)

Dat <- Detect_simulate_variogram(inputdata = dataSub2,
                                 outputdata=data.bathym1, response_var = "BottomDepth")

ggplot(Dat, aes(x = X, y = Y, fill = depth)) +
    geom_raster() +
    scale_fill_viridis_c()
## debugonce(Detect_simulate_variogram)


### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 100
### End:

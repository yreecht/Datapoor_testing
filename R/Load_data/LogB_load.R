#-*- coding: utf-8 -*-

### File: LogB_load.R
### Time-stamp: <2022-12-05 10:02:13 a23579>
###
### Created: 14/10/2022	14:42:00
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################

library(RstoxData)
library(dplyr)
library(tidyr)

salesNv <- RstoxData::readLssFile(file = "./FDIR_HI_LSS_VERDI_2021_PR_2022-10-03.psv", strict = FALSE)
salesNc <- RstoxData::readLssFile(file = "./FDIR_HI_LSS_FANGST_2021_PR_2022-10-03.psv", strict = TRUE)

logB <- RstoxData::readErsFile(file = "./FDIR_HI_ERS_2021_PR_2022-10-03.psv")

head(logB)

head(salesNc)

c("Hovedområde (kode)", "Redskap (kode)", "Redskap - hovedgruppe (kode)",
  "Redskap", "Rundvekt", "Art FAO", "Art FAO (kode)", "Art - FDIR")

table(salesNc$`Redskap`, salesNc$`Redskap (kode)`)
table(salesNc$`Redskap - hovedgruppe (kode)`)

table(salesNc$`Art FAO`)

areaSub <- c("00", "02", "03", "04", "05", "06", "07", "08", "09", "28")

spSub <- c("Gråsteinbit", "Torsk", "Sei", "Hyse", "Uer (vanlig)")

gearSub <- c("20", "22") ## , "21"

salesNc %>%
    filter(`Redskap (kode)` %in% gearSub,
           `Art FAO` %in% spSub,
           `Hovedområde (kode)` %in% areaSub) %>%
    group_by(`Art FAO`) %>%
    summarise(Rundvekt = sum(Rundvekt, na.rm = TRUE)) %>%
    mutate(landings_t = Rundvekt * 1e-3) %>%
    bind_cols(model_est = c(35.02, 1762, 33704, 84666, 1416),
    ## bind_cols(model_est = c(4.40, 25393, 55663, 171176, 4639),
              ## model_q = c(2e-4, rep(0.05, 4))) %>%
              model_q = c(0.00158, 0.0035, 0.0275, 0.0245,
                          0.0165)) %>%
    mutate(ratio = landings_t / model_est,
           new_q = ratio * model_q )
## head(2)

salesNc2 <- RstoxData::readLssFile(
                           ## file = "./FDIR_HI_LSS_FANGST_2016_PR_2017-10-31.psv",
                           ## file = "./FDIR_HI_LSS_FANGST_2011_PR_2016-12-11.psv",
                           ## file = "./FDIR_HI_LSS_FANGST_2021_PR_2022-10-03.psv",
                           file = "./FDIR_HI_LSS_FANGST_2020_PR_2021-11-02.psv",
                           strict = TRUE)

areaSub <- c("40", "41", "42", "08", "09", "28")


spSub <- c("Torsk")

table(salesNc2$Områdegruppering)



salesNc2 %>%
    filter(##`Redskap (kode)` %in% gearSub,
           `Art FAO` %in% spSub,
           `Hovedområde (kode)` %in% areaSub) %>%
    group_by(`Art FAO`, `Kyst/hav (kode)`) %>%
    summarise(Rundvekt = sum(Rundvekt, na.rm = TRUE)) %>%
    mutate(landings_t = Rundvekt * 1e-3)

### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 100
### End:

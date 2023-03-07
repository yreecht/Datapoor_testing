#-*- coding: utf-8 -*-

### File: S3_def.R
### Time-stamp: <2023-03-07 14:31:25 a23579>
###
### Created: 07/03/2023	10:10:02
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################

setGlobalParam <- function(simSettings,
                           parameter,
                           value,
                           constrained = TRUE,
                           constr.length = 1,
                           verbose = TRUE)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  7 Mar 2023, 11:15

    ## Possible parameters:
    parameter <- match.arg(parameter,
                           c("year.depth.restriction", "depth.restriction",
                             "fish.mvt", "Nregion", "Nvessels", "CV_effort",
                             "CV_vessel", "vessel_seeds", "vessel_concentration_factor",
                             "do.tweedie", "Preference", "Changing_preference",
                             "PSU", "SSU", "PSU_prob", "SSU_prob",
                             "seed", "sampling_months", "sampling_startyear",
                             "SSU_change_interval", "samp_mincutoff",
                             "plotting", "parallel", "Interactive"))

    ## Abord if inconsistent length
    if (constrained)
    {
        stopifnot(length(value) == constr.length)
    }

    simSettings[[parameter]] <- value

    return(simSettings)
}

setMonthlyParam <- function(simSettings,
                            parameter,
                            value,
                            verbose = TRUE)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  7 Mar 2023, 11:15

    parameter <- match.arg(parameter,
                           c("Effort_alloc_month"))

    stopifnot("simSettings" %in% class(simSettings))

    if (length(value) == 1)
    {
        if (verbose) message(parameter, " initialised with constant value across months")
        value <- rep(value, 12)
    }else{
        ## Abord if inconsistent length
        stopifnot(length(value) == 12)
    }

    simSettings[[parameter]] <- value

    return(simSettings)
}

setYearlyParam <- function(simSettings,
                            parameter,
                            value,
                            verbose = TRUE)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  7 Mar 2023, 11:15

    parameter <- match.arg(parameter,
                           c("Tot_effort_year"))

    stopifnot("simSettings" %in% class(simSettings))

    if (length(value) == 1)
    {
        if (verbose) message(parameter, " initialised with constant value across years")
        value <- rep(value, simSettings$n_years)
    }else{
        ## Abord if inconsistent length
        stopifnot(length(value) == simSettings$n_years)
    }

    simSettings[[parameter]] <- value

    return(simSettings)
}


setSpeciesParam <- function(simSettings,
                            parameter,
                            value,
                            verbose = TRUE)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  7 Mar 2023, 11:15

    parameter <- match.arg(parameter,
                           c("func_mvt_dist", "func_mvt_depth",
                             "func_mvt_lat", "B0", "Pop_ratio_start",
                             "r", "sigma_p", "sigma_p_timing",
                             "qq_original", "catch_trunc", "xi",
                             "phi", "Discard_rate_beta1", "Discard_rate_beta2",
                             "Discard_survival_rate1", "Discard_survival_rate2"))

    stopifnot("simSettings" %in% class(simSettings))

    if (length(value) == 1)
    {
        value <- rep(value, simSettings$n_species)
    }else{
        ## In case values are named after species, can subset and order:
        if (! is.null(names(value)) &&
            all(simSettings$species_names) %in% names(value))
        {
            value <- value[simSettings$species_names]
        }

        ## Abord if inconsistent length
        stopifnot(length(value) == simSettings$n_species)
    }

    simSettings[[parameter]] <- value

    return(simSettings)
}

setSpeciesMonthlyParam <- function(simSettings,
                                   parameter,
                                   value,
                                   verbose = TRUE)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  7 Mar 2023, 11:15

    parameter <- match.arg(parameter,
                           c("price_fish", "Fish_dist_par1",
                             "Fish_depth_par1", "Fish_depth_par2",
                             "Rangeshift_catchability_adjust"))

    stopifnot("simSettings" %in% class(simSettings))

    if (is.vector(value))
    {
        if (length(value) == 1)
        {
            if (verbose) message(parameter, " initialised with the same value for all species and months")
            value <- rep(value, simSettings$n_species)
        }else{
            ## In case values are named after species, can subset and order:
            if (! is.null(names(value)) &&
                all(simSettings$species_names %in% names(value)))
            {
                value <- value[simSettings$species_names]
            }

            ## Abord if inconsistent length
            stopifnot(length(value) == simSettings$n_species)

            if (verbose) message(parameter, " constant across months")
        }

        value <- matrix(value, ncol = length(value), nrow = 12, byrow = TRUE)
    }

    if (is.matrix(value))
    {
        ## In case value columns are named after species, can subset and order:
        if (! is.null(colnames(value)) &&
            all(simSettings$species_names %in% colnames(value)))
        {
            value <- value[ , simSettings$species_names]
        }

        ## Abord if inconsistent length
        stopifnot(dim(value) == c(12, simSettings$n_species))
    }

    simSettings[[parameter]] <- value

    return(simSettings)
}

new_simSettings <- function(n_species = 4,
                            species_names = paste0("Sp", seq(length.out = as.numeric(n_species))),
                            n_years = 15,
                            Range_X = 1:20,
                            Range_Y = 1:40,
                            SD_O = 100,
                            SpatialScale = 3,
                            input_habitat = NULL,
                            ...)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  7 Mar 2023, 10:14

    ## New object to be populated:
    simSet <- structure(list(), class = "simSettings")

    ## browser()

    ## Check consistency in case both n_species and species_names
    ##   are passed by the user:
    if(! missing(n_species) && ! missing(species_names))
    {
        stopifnot(n_species == length(species_names))
    }

    Nsp <- length(species_names)

    simSet$n_species <- Nsp
    simSet$species_names <- as.character(species_names)
    simSet$n_years <- as.numeric(n_years)

    ## Habitat characteristics
    if (missing(input_habitat) || is.null(input_habitat))
    {
        simSet$Range_X <- Range_X
        simSet$Range_Y <- Range_Y
        simSet$SD_O <- SD_O
        simSet$SpatialScale <- SpatialScale
    }else{
        ## In case a habitat matrix is provided:
        simSet$Range_X <- seq(length.out = ncol(input_habitat))
        simSet$Range_Y <- seq(length.out = nrow(input_habitat))
    }

    ## Getting optional supplementary parameters:
    dots <- list(...)

    ## Supplied or default value for further parameters
    if ("price_fish" %in% names(dots))
    {
        price_fish <- dots$price_fish
    }else{
        price_fish <- matrix(1, ncol = Nsp, nrow = 12)
    }

    if ("Fish_dist_par1" %in% names(dots))
    {
        Fish_dist_par1 <- dots$Fish_dist_par1
    }else{
        Fish_dist_par1 <- matrix(3, ncol = Nsp, nrow = 12)
    }

    if ("func_mvt_dist" %in% names(dots))
    {
        func_mvt_dist <- dots$func_mvt_dist
    }else{
        func_mvt_dist <- "Exp"
    }

    if ("func_mvt_depth" %in% names(dots))
    {
        func_mvt_depth <- dots$func_mvt_depth
    }else{
        func_mvt_depth <- "Lognorm"
    }

    if ("func_mvt_lat" %in% names(dots))
    {
        func_mvt_lat <- dots$func_mvt_lat
    }else{
        func_mvt_lat <- "Unif"
    }

    if ("Fish_dist_par1" %in% names(dots))
    {
        Fish_dist_par1 <- dots$Fish_dist_par1
    }else{
        Fish_dist_par1 <-  matrix(1, ncol = Nsp, nrow = 12)
    }

    if ("Fish_dist_par2" %in% names(dots))
    {
        Fish_dist_par2 <- dots$Fish_dist_par2
    }else{
        Fish_dist_par2 <- rep(0, Nsp)
    }

    if ("parm" %in% names(dots))
    {
        parm <- dots$parm
    }else{
        parm <- <def>
    }

    if ("Fish_depth_par1" %in% names(dots))
    {
        Fish_depth_par1 <- dots$Fish_depth_par1
    }else{
        Fish_depth_par1 <- matrix(100, ncol = Nsp, nrow = 12)
    }

    if ("Fish_depth_par2" %in% names(dots))
    {
        Fish_depth_par2 <- dots$Fish_depth_par2
    }else{
        Fish_depth_par2 <- matrix(0.5, ncol = Nsp, nrow = 12)
    }

    if ("Rangeshift_catchability_adjust" %in% names(dots))
    {
        Rangeshift_catchability_adjust <- dots$Rangeshift_catchability_adjust
    }else{
        Rangeshift_catchability_adjust <- matrix(1, ncol = Nsp, nrow = 12)
    }

    if ("B0" %in% names(dots))
    {
        B0 <- dots$B0
    }else{
        B0 <- rep(1e6, Nsp)
    }

    if ("Pop_ratio_start" %in% names(dots))
    {
        Pop_ratio_start <- dots$Pop_ratio_start
    }else{
        Pop_ratio_start <- rep(0.5, Nsp)
    }

    if ("r" %in% names(dots))
    {
        r <- dots$r
    }else{
        r <- rep(0.01, Nsp)
    }

    if ("sigma_p" %in% names(dots))
    {
        sigma_p <- dots$sigma_p
    }else{
        sigma_p <- rep(0.5, Nsp)
    }

    if ("sigma_p_timing" %in% names(dots))
    {
        sigma_p_timing <- dots$sigma_p_timing
    }else{
        sigma_p_timing <- rep(2, Nsp)
    }

    if ("fish.mvt" %in% names(dots))
    {
        fish.mvt <- dots$fish.mvt
    }else{
        fish.mvt <- TRUE
    }

    if ("Nregion" %in% names(dots))
    {
        Nregion <- dots$Nregion
    }else{
        Nregion <- c(2,2)
    }

    if ("Nvessels" %in% names(dots))
    {
        Nvessels <- dots$Nvessels
    }else{
        Nvessels <- 50
    }

    if ("CV_effort" %in% names(dots))
    {
        CV_effort <- dots$CV_effort
    }else{
        CV_effort <- 0.2
    }

    ## if ("parm" %in% names(dots))
    ## {
    ##     parm <- dots$parm
    ## }else{
    ##     parm <- <def>
    ## }

    ## Set parameters:
    simSet <- setSpeciesMonthlyParam(simSet, "price_fish", price_fish, verbose = FALSE)

    simSet <- setSpeciesParam(simSet, "func_mvt_dist", func_mvt_dist, verbose = FALSE)
    simSet <- setSpeciesParam(simSet, "func_mvt_depth", func_mvt_depth, verbose = FALSE)
    simSet <- setSpeciesParam(simSet, "func_mvt_lat", func_mvt_lat, verbose = FALSE)
    simSet <- setSpeciesMonthlyParam(simSet, "Fish_dist_par1", Fish_dist_par1, verbose = FALSE)


    return(simSet)
}



test <- new_simSettings(species_names = c("haddock", "saithe",
                                  "cod", "redfish",
                                  "wolfish_f", "wolfish_m"),
                func_mvt_dist = "Exp",
                func_mvt_depth = c(rep("Lognorm", 6)),
                price_fish = c(13.33, # 2021 average (round, fresh), https://www.rafisklaget.no/statistikk-detaljer
                               8.56,
                               23.89,
                               12.61,
                               rep(12.16, 2)),
                Fish_dist_par1 = matrix(c(5, 10, 6, 5, 3, 5, # mvt distance mean - their mobility within month1 only values >0
                                          5, 12, 6, 5, 3, 3, # month2   only values >0 Values for target species somewhat arbitrary.
                                          5, 12, 6, 5, 3, 3, # month3   only values >0
                                          5, 12, 6, 5, 3, 3, # month4   only values >0
                                          5, 10, 6, 5, 3, 3, # month5   only values >0
                                          5, 10, 6, 5, 5, 5, # month6   only values >0 Increased mobility before (and after) reproduction
                                          5, 10, 6, 5, 2, 2, # month7   only values >0
                                          5, 10, 6, 5, 2, 0.5, # month8   only values >0 males A. wolffish not moving while caring for offspring.
                                          5, 10, 6, 5, 2, 0.5, # month9   only values >0
                                          5, 10, 6, 5, 2, 0.5, # month10  only values >0
                                          5, 10, 6, 5, 5, 0.5, # month11  only values >0
                                          5, 10, 6, 5, 3, 2),# month12  only values >0
                                        nrow=12, ncol = 6, byrow = TRUE))

test

test2 <- new_simSettings(species_names = c("haddock", "saithe",
                                           "cod", "redfish",
                                           "wolfish_f", "wolfish_m"))

test2

test3 <- setSpeciesMonthlyParam(test2, "price_fish",
                                c(13.33, 8.56, 23.89, 12.61, rep(12.16, 2)))

test4 <- setSpeciesMonthlyParam(simSettings = test2,
                                param = "price_fish",
                                value = c("haddock" = 13.33, "saithe" = 8.56,
                                          "redfish" = 12.61, "cod" = 23.89,
                                          "wolfish_f" = 12.16, "wolfish_m" = 12.16,
                                          "Other" = 2))

test4

identical(test3, test4)

### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 100
### End:

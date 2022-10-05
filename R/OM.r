############################################################################################################
###
### Spatial explicit operating model for modeling fishery and species "monitoring"
###
### The most important part is to define the fish abundance, price, and habitat overlap --> need to make it more explicit
### --> because this will create the basis for trade-offs and targeting change
### think about changing correlation by year (or random effect / random walk for the cor par)
### streamline with furrr the run of scenario & models
### Make more documentation as in sdmTMB

############################################################################################################

rm(list=ls())
gc()

#### Loading libraries
pacman::p_load(parallel, MASS, RandomFields, fields, geoR, gtools, tweedie, ggplot2, tidyverse, ggnewscale, TMB, TMBhelper, sdmTMB)


#### sourcing codes
	source("R/Generate_scenario_data.R")  ## File that runs the simulation
  source("R/Functions.R")               ## File that include the main functions that are called in various part of the simulation
  source("R/Scenario_setup.R")            ## The main file to set-up the operating model i.e. configuring the population parameters and the fleet dynamics


	# some ideas of setting
	# Sim2$Fish_depth_par1 = c(120, 250, 150, 500) # bycatch species has similar distribution as another species
	# Sim_Settings3$Fish_depth_par1 = c(120, 250, 450, 200) # bycatch species is distributed deeper than the main target
	Sim2 <- Sim1
	Sim2$SD_O = 150
	Sim2$SpatialScale = 15
	Sim2$sigma_p= c(1, 1, 1, 0.5)
	Sim2$CV_vessel= 0.1

	## Running the simulation model with the user-specified configurations
	system.time(
	  Data <- Generate_scenario_data(Sim_Settings = Sim2, seed_input=12)
  )


	Data$Data <- as.data.frame(Data$Data)

	# The range of values of catch for each species
	apply(Data$Data, 2, range)

	# Making a histogram of the catch
	Data_long <- pivot_longer(Data$Data, cols = starts_with("Sp"), names_to = "Species", values_to = "Catch")
	ggplot(Data_long, aes(x=Catch)) + facet_grid(.~Species) + geom_histogram() + theme_bw()

	# Looking at the catch composition
	Comps_data <- Data$Data %>% dplyr::select(starts_with("Sp")) %>% mutate( . / rowSums(.))
	apply(Comps_data,2,mean)


	# Simulated depth distribution
	ggplot(Data$bathym) + geom_raster(aes(x=X, y=Y, fill=depth)) + scale_fill_viridis_c()

	# Simulated depth distribution
	# ggplot(Data$Biomass) + geom_raster(aes(x=X, y=Y, fill=Sp1)) + facet_wrap(~year) + scale_fill_viridis_c()
  # Simulated effort distribution in space
	ggplot(data = Data$Data %>% filter(year == 10, month==11)) + geom_raster(data=Data$bathym, aes(x=X, y=Y, fill=depth)) + scale_fill_viridis_c() +
	  geom_point(aes(x=X, y=Y))


#### Calculate the true index
	Average_monthly_biomass <- apply(Data$Biomass, c(1,2,4), sum)
	plot(1:Sim2$n_years, Average_monthly_biomass[,1,1]/Average_monthly_biomass[1,1,1], type="l", ylim=c(0,1.5))
	lines(1:Sim2$n_years, Average_monthly_biomass[,1,2]/Average_monthly_biomass[1,1,2], col = "red")
	lines(1:Sim2$n_years, Average_monthly_biomass[,1,3]/Average_monthly_biomass[1,1,3], col = "blue")
	lines(1:Sim2$n_years, Average_monthly_biomass[,1,4]/Average_monthly_biomass[1,1,4], col = "green")

  true_index <- apply(Average_monthly_biomass[,1,], 2, function(x) x/x[1]) %>% as.data.frame()
  colnames(true_index) <- paste0("Sp", 1:Sim2$n_species)


#### Perform the sample selection
  Final_data <- sampling_select(data = Data$Data %>% as.data.frame(), percent = Sim2$samp_prob, unit=Sim2$samp_unit, seed=Sim2$samp_seed, months = c(11:12))
  Final_data <- sampling_select(data = Data$Data %>% as.data.frame(), percent = 0.2, unit=Sim2$samp_unit, seed=2, months = c(11:12))

  Year_adj <- 1   # if the data is taken towards the end of the year, add the year adjustement factor

  nb_x = 2
  nb_y = 2

  Final_data <- equal_partition(Final_data, nb_x, nb_y, Range_X = Sim2$Range_X, Range_Y = Sim2$Range_Y)
  Final_data$depth_scl <- scale(Final_data$depth)

  truncate_fn <- function(x) ifelse(x < Sim2$samp_mincutoff, 0, x)
  Final_data <- Final_data %>% filter(year >= Sim2$start_year) %>% # to remove the initial year effect of the simulator ("burn-in" time)
    mutate_at(vars(starts_with("Sp")), truncate_fn) # truncate variables as in real world
  Final_data <- Final_data %>% mutate(year_fct = as.factor(year),
                          area_fct = as.factor(area),
                          month_fct = as.factor(month),
                          vessel_fct = as.factor(vessel))


  Final_data_df <- Final_data %>% pivot_longer(cols = starts_with("Sp"), names_to = "Species", values_to = "CPUE_trunc")
  table(Final_data_df$CPUE_trunc==0)/nrow(Final_data_df)
  ggplot(Final_data_df %>% filter(Species == "Sp1"), aes(x=CPUE_trunc)) + geom_histogram(bins=100) + theme_bw()
  ggplot(Final_data_df %>% filter(Species == "Sp2"), aes(x=CPUE_trunc)) + geom_histogram(bins=100) + theme_bw()
  ggplot(Final_data_df %>% filter(Species == "Sp3"), aes(x=CPUE_trunc)) + geom_histogram(bins=100) + theme_bw()
  ggplot(Final_data_df %>% filter(Species == "Sp4"), aes(x=CPUE_trunc)) + geom_histogram(bins=100) + theme_bw()
  plot_fishing(Final_data, years = seq(1, Sim2$n_years, by=5))


#### projection grid
  qcs_grid <- do.call(rbind, replicate(Sim2$n_years, Data$bathym, simplify=FALSE))
  qcs_grid$year <- as.numeric(rep(1:Sim2$n_years, each=nrow(Data$bathym)))
  qcs_grid$vessel <- unique(Final_data_df$vessel)[1]
  qcs_grid$vessel_fct <- unique(Final_data_df$vessel)[1]
  qcs_grid <-  qcs_grid %>% mutate(X = as.numeric(X), Y= as.numeric(Y), depth_scl = (depth-mean(qcs_grid$depth))/sd(qcs_grid$depth))
  qcs_grid$CPUE_trunc <- 1
  qcs_grid <- equal_partition(qcs_grid, nb_x, nb_y, Range_X = Sim2$Range_X, Range_Y = Sim2$Range_Y)
  qcs_grid$year_area_fct <- as.factor(apply(qcs_grid[, c('year','area')], 1, function(x) paste(x, collapse="_")))
  qcs_grid <- qcs_grid %>% filter(year >= Sim2$start_year)

  # now adding the average species catch in each year, area combination
  qcs_grid <- qcs_grid %>% left_join(
    Final_data %>% pivot_longer(cols = starts_with("Sp"), names_to = "Species", values_to = "CPUE") %>%
      group_by(year, area, Species) %>% summarize(Mean = mean(CPUE)) %>% ungroup() %>%
      pivot_wider(names_from = "Species", values_from = "Mean"))
  qcs_grid <- qcs_grid %>% mutate(area_fct = factor(area), year_fct = as.factor(year), yy=1)



##### generating single species data (for single species models)
  Final_data_df$vessel_fct <- as.factor(Final_data_df$vessel)

  Which_sp <- "Sp4"

  Final_data_bycatch = Final_data_df %>% filter(Species == Which_sp)
  Final_data_bycatch$year_fct <- as.factor(Final_data_bycatch$year)
  Final_data_bycatch <- Final_data_bycatch %>% mutate(depth_scl = (depth-mean(qcs_grid$depth))/sd(qcs_grid$depth))
  Final_data_bycatch$year_area_fct <- as.factor(apply(Final_data_bycatch[, c('year','area')], 1, function(x) paste(x, collapse="_")))




########### Now that we have generated the data, it is time to fit different models
## List of models to test:
# gam
# glmmTMB (tweedie)
# glmmTMB - zeroinflated models
# student-t distribution with long tails with possible zero-inflation (code in TMB)
# mixture model
# randomforest: ranger
# boosted regression tree: gbm
# sdmTMB
# multispecies model with fishing behavior (tweedie)

## to include
  # flexSDM package to test (not a priority but the types of model included)
  # neural network: neuralnet & keras
  # spatial model with fishing tactics

####### gam:
  library(mgcv)
  gam1 <- gam(CPUE_trunc ~ 0 + as.factor(year) + as.factor(area) +  s(depth_scl) +
                s(vessel_fct, bs="re"), data=Final_data_bycatch, family = tw)
  gam1 <- gam(CPUE_trunc ~ 0 + as.factor(year) + s(X,Y) +  s(depth_scl) +
                s(vessel_fct, bs="re"), data=Final_data_bycatch, family = tw)
  # gam2 <- gam(Sp3 ~ 0 + as.factor(year) + as.factor(area) +  s(depth_scl) +
                # s(Sp1,Sp2), data=Final_data, family = tw)
  gam.check(gam1)
  predicted <- predict(gam1, newdata = qcs_grid, type="response", se.fit=F)
  pred <- qcs_grid
  pred$pred = predicted
  pred <- pred %>% mutate(year_or = year, year = year_or + Year_adj)
  IA_gam <- pred %>% group_by(year) %>% summarize(IA = sum(pred, na.rm=T)) %>% mutate(Method="GAM", IA = IA/IA[1])
  ggplot(IA_gam, aes(x=year, y=IA)) + labs(x = "Year", y = "Index") +
    geom_line(lwd = 1, colour = "grey30")  + ggtitle("GAM") +
    geom_line(data = data.frame(IA=true_index[-c(1:(Sim2$start_year-1)),Which_sp]/true_index[Sim2$start_year+Year_adj,Which_sp], year = Sim2$start_year:Sim2$n_years), col="red", size=2)+
    theme_bw()+ coord_cartesian(ylim=c(0,1.5))


####### random forest:
  library(ranger)
  rf1 <- ranger(CPUE_trunc ~ ., data=Final_data_bycatch %>% dplyr::select(year, area, depth, vessel, CPUE_trunc),
                num.trees = 5000)
  predicted <- predict(rf1, data = qcs_grid, type="response")
  pred <- qcs_grid
  pred$pred = predicted$predictions
  pred <- pred %>% mutate(year_or = year, year = year_or + Year_adj)
  IA_rf <- pred %>% group_by(year) %>% summarize(IA = sum(pred)) %>% mutate(Method="RF", IA = IA/IA[1])
  ggplot(IA_rf, aes(x=year, y=IA)) + labs(x = "Year", y = "Index") +
    geom_line(lwd = 1, colour = "grey30")  + ggtitle("Random Forest") +
    geom_line(data = data.frame(IA=true_index[-c(1:(Sim2$start_year-1)),Which_sp]/true_index[Sim2$start_year+Year_adj,Which_sp], year =Sim2$start_year:Sim2$n_years), col="red", size=2)+
    theme_bw()+ coord_cartesian(ylim=c(0,1.5))


####### glmmTMB
  library(glmmTMB)
  library(splines)
  if (Which_sp != "Sp4") lm1 <- glmmTMB::glmmTMB(CPUE_trunc ~ 0 + as.factor(year) + as.factor(area) + ns(depth_scl) + (1|vessel_fct), #+ (1|year_area_fct)
                          data=Final_data_bycatch, family = tweedie(link = "log"))
  if (Which_sp == "Sp4") lm1 <- glmmTMB::glmmTMB(CPUE_trunc ~ 0 + as.factor(year) + as.factor(area) + ns(depth_scl) + (1|vessel_fct), #+ (1|year_area_fct)
                          data=Final_data_bycatch, family = compois(link = "log"))
  if (Which_sp == "Sp4") lm1 <- glmmTMB::glmmTMB(CPUE_trunc ~ as.factor(year) + ns(depth_scl) + (1|vessel_fct), #+ (1|year_area_fct)
                          data=Final_data_bycatch, family = nbinom2(link = "log"))

  sim <- DHARMa::simulateResiduals(lm1)
  plot(sim)
  predicted <- predict(lm1, newdata = qcs_grid, type="response", se.fit=F)
  pred <- qcs_grid
  pred$pred = predicted
  pred <- pred %>% mutate(year_or = year, year = year_or + Year_adj)
  IA_glmmTMB <- pred %>% group_by(year) %>% summarize(IA = sum(pred)) %>% mutate(Method="GlmmTMB", IA = IA/IA[1])
  ggplot(IA_glmmTMB, aes(x=year, y=IA)) + labs(x = "Year", y = "Index") +
    geom_line(lwd = 1, colour = "grey30")  + ggtitle("GLMMTMB") +
    geom_line(data = data.frame(IA=true_index[-c(1:(Sim2$start_year-1)),Which_sp]/true_index[Sim2$start_year+Year_adj,Which_sp], year =Sim2$start_year:Sim2$n_years), col="red", size=2)+
    theme_bw()+ coord_cartesian(ylim=c(0,1.5))


####### generalized boosted regression tree
  library(gbm)
  gbm1 <- gbm(CPUE_trunc ~ ., data=Final_data_bycatch %>% select(year, area, depth, vessel, CPUE_trunc),
              distribution = "gaussian",  n.trees = 5000)
  predicted <- predict(gbm1, newdata = qcs_grid, type="response")
  pred <- qcs_grid
  pred$pred = predicted
  IA_gbm <- pred %>% group_by(year) %>% summarize(IA = sum(pred)) %>% mutate(Method="GBM", IA = IA/IA[1])
  ggplot(IA_gbm, aes(x=year, y=IA)) + labs(x = "Year", y = "Index") +
    geom_line(lwd = 1, colour = "grey30") + ggtitle("GBM") +
    geom_line(data = data.frame(IA=true_index[-c(1:(Sim2$start_year-1)),Which_sp]/true_index[Sim2$start_year+Year_adj,Which_sp], year =Sim2$start_year:Sim2$n_years), col="red", size=2)+
    theme_bw()+ coord_cartesian(ylim=c(0,1.5))


  ####### Spatial model: sdmTMB to the data
  library(sdmTMB)
  mesh <- make_mesh(Final_data_bycatch, xy_cols = c("X", "Y"), cutoff = 5)
  plot(mesh)

  fit <- sdmTMB(CPUE_trunc ~ 0 + as.factor(year) + s(depth_scl, k=3) ,
                data = Final_data_bycatch, spatial = "on", spatiotemporal = "iid", time = "year", family = tweedie(link="log"), mesh = mesh)
  fit
  res_test <- residuals(fit)
  qqnorm(res_test);qqline(res_test)

  p_st <- predict(fit, newdata = qcs_grid,
                  return_tmb_object = TRUE)
  index <- get_index(p_st)
  index <- index %>% mutate(Method="sdmTMB", IA = est/est[1])
  ggplot(index, aes(x=year, y=IA)) +
    # geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey90") +
    geom_line(lwd = 1, colour = "grey30") +
    labs(x = "Year", y = "Biomass (kg)") +
    geom_line(data = data.frame(IA=true_index[-c(1:(Sim2$start_year-1)),Which_sp]/true_index[Sim2$start_year+Year_adj,Which_sp], year = Sim2$start_year:Sim2$n_years), col="red", size=2)+
    theme_bw()#+coord_cartesian(ylim=c(0, 20000))


####### spatial DFA (similar to VAST)
#### Using the multispecies model and can include targeting behavior

  source("R/utils.R")
  source("R/Prepare_data.R")
  source("R/fit_sdm.R")

  # prepare model configuration file
  conf <- list(
    model_formula      = formula(yy ~ year_fct + s(depth, k=3) + as.factor(area) + (1|vessel_fct)),  # this works
    spatial_model      = 0L,    # 0 = no, non-spatial model, 1 = yes, make it a spatial model
    include_st         = 0L,    # include the spatio-temporal component? 0 = no, make a simple spatial model with only the average spatial field, 1: yes, st follows IID, 2: yes, st follows AR1
    barrier            = 0L,    # using a barrier model?
    do_spatialDFA      = 1L,    # apply DFA method to model the spatial effects? 0 = no, 1 = yes
    incl_target        = 0L,    # including or not the small scale targeting behavior
    incl_sp_int        = 0L,    # including or not the multispecies catch interaction term (not estimable at the moment)
    Nfactor            = 2L,    # how many factor you want to use to reduce the dimensionality of the problem
    Do_predict         = 1L,    # do you want to predict the CPUE values on the prediction grid cells
    se_fit             = 0L,    # do you want to calculate se on the above prediction?
    sim                = 0L,    # 0: no simulation, 1: yes
    INLAmesh_cutoff    = 5      # the minimum distance between meshes (default to 5), lower it when more data points and vice versa
  )

  # prepare data
  dat_use = Prepare_data(data=Final_data, conf=conf, Pred_data=qcs_grid)

  # Fit model
    startTime = Sys.time()
    run = fit_sdm(data=dat_use, conf)
    endTime = Sys.time()
    endTime-startTime


  run$opt$diagnostics
  run$opt$max_gradient

  betas <- t(run$obj$report(par = run$obj$env$last.par)$beta)
  colnames(betas) <- colnames(dat_use$tmb_data$X)

  IAs <- run$obj$report(par = run$obj$env$last.par)$mu_proj
  colnames(IAs) <- paste0("Sp", 1:Sim2$n_species)

  # now plotting
  pred <- qcs_grid
  pred$pred = IAs[,Which_sp]
  IA_multsp <- pred %>% group_by(year) %>% summarize(IA = sum(pred)) %>% mutate(Method="Multsp", IA = IA/IA[1]) %>% mutate(year_or = year, year = year_or + Year_adj)
  ggplot(IA_multsp, aes(x=year, y=IA/IA[1])) + labs(x = "Year", y = "Index") +
    geom_line(lwd = 1, colour = "grey30")  + ggtitle("Multsp") +
    geom_line(data = data.frame(IA=true_index[-c(1:(Sim2$start_year-1)),Which_sp]/true_index[Sim2$start_year+Year_adj,Which_sp], year = Sim2$start_year:Sim2$n_years), col="red", size=2)+
    theme_bw() + coord_cartesian(ylim=c(0,1.5))  # +
  # geom_line(data = data.frame(year = Sim2$start_year:Sim2$n_years, IA = IA_multsp$IA[1]*exp(betas[as.numeric(gsub("Sp", "", Which_sp)),1]+c(0,as.numeric(betas[as.numeric(gsub("Sp", "", Which_sp)), grep("year", colnames(betas))])))), col="green")


#### Using the multispecies model via the species correlation

  source("R/utils.R")
  datdat <- Final_data
  datdat$yy <- as.numeric(datdat[,'Sp1'])
  datdat <- datdat %>% mutate(year_fct = as.factor(year),
                              area_fct = as.factor(area),
                              month_fct = as.factor(month),
                              vessel_fct = as.factor(vessel))

  model_formula <- formula(yy ~ area_fct + year_fct + s(depth, k=3)) # + (1|vessel_fct))  # this works
  fixed_effect_formula = nobars(model_formula)
  RE_effects = barnames(findbars(model_formula))
  mgcv_mod <- mgcv::gam(fixed_effect_formula, data = datdat)
  X <- model.matrix(mgcv_mod)


  RE_indexes = as.matrix(0L)
  nobs_RE = 0L
  ln_tau_G_index= 0L
  if (length(RE_effects)>0)
  {
    RE_indexes <- sapply(RE_effects, function(x) as.numeric(factor(datdat[[x]], labels=1:length(unique(datdat[[x]]))))) - 1L
    RE_indexes <- as.matrix(RE_indexes)
    nobs_RE <- unname(apply(RE_indexes, 2L, max)) + 1L
    if (length(nobs_RE) == 0L) nobs_RE <- 0L
    ln_tau_G_index <- unlist(lapply(seq_along(nobs_RE), function(i) rep(i, each = nobs_RE[i]))) - 1L
  }

  yobs= datdat %>% select(starts_with("sp"))
  Nspecies = ncol(yobs)


  ## Creating the projection data i.e. combination of different variables. Need to choose a reference category
  # Pred_data <- expand.grid(area_fct = levels(datdat$area_fct), year_fct = levels(datdat$year_fct),
  #                          month_fct = levels(datdat$month_fct)[1],
  #                          vessel_fct = levels(datdat$vessel_fct)[1],
  #                          depth = mean(datdat$depth), yy=1)
  qcs_grid <- qcs_grid %>% mutate(area_fct = factor(area), year_fct = as.factor(year), yy=1)
  Pred_data <- qcs_grid
  X_pred <- mgcv::predict.gam(mgcv_mod, type = "lpmatrix", newdata = Pred_data)
  RE_indexes_proj = as.matrix(0L)
  if (length(RE_effects)>0)
  {
    temp <- data.frame(val=unique(datdat[, RE_effects]), new_val = 1:length(unique(datdat[, RE_effects])) - 1)
    RE_indexes_proj <- as.matrix(temp$new_val[match(Pred_data[, RE_effects], temp$val)])
    RE_indexes_proj[is.na(RE_indexes_proj)] <- 999
    RE_indexes_proj = as.matrix(RE_indexes_proj)
  }

  tmb_data <- list(
    Nfactor = 3,
    X = as.matrix(X),
    yobs = as.matrix(yobs),
    RE_indexes = RE_indexes,
    nobs_RE = nobs_RE,
    ln_tau_G_index = ln_tau_G_index,
    # to_keep = matrix(1, nrow(yobs), Nspecies),
    family = 0,
    link = 1,
    sim = 0,
    incl_target = 0, # including or not the small scale targeting behavior
    incl_sp_int = 1,  # including or not the multispecies catch interaction term (not estimable at the moment)
    Do_predict = 1,
    se_fit = 0,
    Npred = nrow(X_pred),
    X_proj = as.matrix(X_pred),
    RE_indexes_proj = RE_indexes_proj,
    target_index = 69
  )


  vals = rnorm(Nspecies*(tmb_data$Nfactor))
  vals = rnorm(Nspecies*(tmb_data$Nfactor) - (tmb_data$Nfactor)*(tmb_data$Nfactor-1)/2)

  tmb_param <- list(
    beta        = matrix(0, nrow=ncol(X), ncol=Nspecies),
    L_val_target= vals,
    L_val_sp    = rnorm(Nspecies*(Nspecies-1)/2) ,
    RE_species_latent  = matrix(0.1, nrow=nrow(yobs), ncol=tmb_data$Nfactor) ,
    RE_species_int  = matrix(0.1, nrow=nrow(yobs), ncol=Nspecies) ,
    logsds      = rep(0.1, Nspecies),
    thetaf      = rep(0.1, Nspecies) ,
    ln_phi      = rep(0.1, Nspecies) ,
    ln_tau_G    = matrix(0.1, nrow = length(nobs_RE), ncol=Nspecies) ,
    RE          = matrix(0.01, nrow=sum(nobs_RE), ncol=Nspecies)
  )


  library(TMB)
  version <- paste0(getwd(), "/src/mult_species"); dllversion = "mult_species"

  #dyn.unload(dynlib(version))
  if (.Platform$OS.type == "unix") compile(paste0(version, ".cpp"),"-O0 -g")
  if (.Platform$OS.type == "windows") compile(paste0(version, ".cpp"), "-O1 -g", DLLFLAGS="")
  dyn.load(dynlib(version))

  Map_phase1 <- list()
  if (tmb_data$incl_target == 0) {
    Map_phase1$L_val_target <- factor(rep(NA, length(tmb_param$L_val_target)))
    Map_phase1$RE_species_latent <- factor(rep(NA, length(tmb_param$RE_species_latent)))
    tmb_param$RE_species_latent = matrix(0, nrow=nrow(yobs), ncol=tmb_data$Nfactor)
  }
  if (tmb_data$incl_sp_int == 0) {
    Map_phase1$logsds <- factor(rep(NA, length(tmb_param$logsds)))
    Map_phase1$L_val_sp <- factor(rep(NA, length(tmb_param$L_val_sp)))
    Map_phase1$RE_species_int <- factor(rep(NA, length(tmb_param$RE_species_int)))
    tmb_param$RE_species_int = matrix(0, nrow=nrow(yobs), ncol=Nspecies)
  }

  if(sum(nobs_RE)==0) Map_phase1$ln_tau_G <- factor(rep(NA, length(tmb_param$ln_tau_G)))
  random_list <- NULL
  if(sum(nobs_RE)>0)  random_list <- c(random_list, "RE")
  if(tmb_data$incl_target == 1)  random_list <- c(random_list, "RE_species_latent")
  if(tmb_data$incl_sp_int == 1)  random_list <- c(random_list, "RE_species_int")

  tmb_obj_phase1 <- MakeADFun(data = tmb_data, parameters = tmb_param, map = Map_phase1,
                              random = random_list, DLL = dllversion, silent = TRUE)


  opt_phase1 <- fit_tmb(tmb_obj_phase1, lower=-15, upper=15, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))

  # sdrep <- sdreport(tmb_obj_phase1)
  # fixed <- summary(sdrep, "fixed")
  # IAs <- summary(sdrep, "report")
  # IAs <- data.frame(IAs, dplyr::bind_rows(replicate(Nspecies, Pred_data, simplify = FALSE)), Species = rep(1:Nspecies, each = nrow(Pred_data)))
  tmb_obj_phase1$report(par = tmb_obj_phase1$env$last.par)$Lt
  head(tmb_obj_phase1$report(par = tmb_obj_phase1$env$last.par)$RE_species_latent)

  tmb_obj_phase1$report(par = tmb_obj_phase1$env$last.par)$Cov

  temp <- c()
  for (i in levels(datdat$area_fct)){
    asd <- sapply(which(datdat$area == i), function(x) tmb_obj_phase1$report(par = tmb_obj_phase1$env$last.par)$Lt %*% tmb_obj_phase1$report(par = tmb_obj_phase1$env$last.par)$RE_species_latent[x,])
    temp <- rbind(temp, apply(asd, 1, mean))
  }

  IAs <- tmb_obj_phase1$report(par = tmb_obj_phase1$env$last.par)$mu_proj
  colnames(IAs) <- paste0("Sp", 1:Sim2$n_species)

  # now plotting
  pred <- qcs_grid
  pred$pred = IAs[,Which_sp]
  IA_multsp <- pred %>% group_by(year) %>% summarize(IA = sum(pred))
  ggplot(IA_multsp, aes(x=year, y=IA/IA[1])) + labs(x = "Year", y = "Index") +
    geom_line(lwd = 1, colour = "grey30")  + ggtitle("Multsp") +
    geom_line(data = data.frame(IA=true_index[-c(1:(Sim2$start_year-1)),Which_sp]*IA_multsp$IA[1], year =Sim2$start_year:Sim2$n_years), col="red", size=2)+
    theme_bw() + coord_cartesian(ylim=c(0,1.5))   #+
   # geom_line(data = data.frame(year = 1:Sim2$n_years, IA = IA_multsp$IA[1]*exp(opt_phase1$par[1]+c(0,opt_phase1$par[grep("year", colnames(X))]))), col="green")









############################################################################################################
###
### Operating model for the rare species analysis
###
### To do:
### 1. vessel fishing should be not equal
###
############################################################################################################

rm(list=ls())
gc()

#### Loading libraries
pacman::p_load(MASS, RandomFields, fields, geoR, gtools, tweedie, ggplot2, tidyverse, ggnewscale)


#### sourcing codes
	source("R/Generate_scenario_data.r")

#### Setting OM parameters

### Scenario 1:
### equal price to all species (attractiveness based on expected catch)
	Sim_Settings1 <- list(
	"n_years" = 30,

	### Species parameters
		"n_species"  = 4,														# petrale, dover, darkblotched, sablefish
		"price_fish" = matrix(rep(c(1,1,-1,1),30), ncol=4, byrow=T),

		## Fish habitat preference/movement control paramters
			"func_mvt_dist" = c("Exp", "Exp", "Exp", "Exp"),					# shape of mvt p_function
			"func_mvt_depth" = c("Lognorm", "Lognorm", "Lognorm", "Lognorm"),	# shape of depth preference p_function (choice of Normal, Exp, Lognormal, Uniform)
			"func_mvt_lat" = c("Unif", "Unif", "Unif", "Unif"),					# shape of range preference p_function (choice of Normal, Exp, Lognormal, Uniform)

			"Fish_dist_par1" = c(rep(0.4, 3), 1),							# mvt distance
			"Fish_dist_par2" = rep(0, 4), 										# mvt distance
			"Fish_depth_par1" = c(220, 450, 250, 1200),       # depth preference
			"Fish_depth_par2" = c(1, 1, 0.3, 0.5),            # depth preference
			"Fish_range_par1" = rep(0, 4),                    # X-axis range
			"Fish_range_par2" = rep(50,4),                    # X axis range

		## Species Schaefer pop dyn parameters
			"B0" = c(50000,820000,36000,540000),								# petrale, dover, darkblotched, sablefish (mt)
			"r" = c(0.22, 0.17, 0.08, 0.08),
			"no.mvt" = FALSE,				# whether animals redistribute based on habitat preference

	### Parameters controlling the vessel dynamics
		"Nvessels" = 30,					# nb vessels
		"Tot_effort" = 200,					# end year nb of effort
		"CV_effort" = 0.2, 					# CV of effort around the yearly mean effort
		"qq_original" = 0.2,				# the average catchability coef for ALL vessels
		"CV_vessel" = 0.05,					# the between vessel variability in mean catchability
		"do.tweedie" = TRUE,				# include observation error in catchability
		"xi" = c(1.9,1.9,1.2,1.2),   # power of the tweedie distribution
		"phi" = c(0.2,0.2,0.2,0.2),  # the scaler of the tweedie distributin
		"Preference" = 1, 					# controls how effort concentrate into areas
		"Changing_preference" = FALSE, 		# whether effort concentration changes over time

	### Habitat (depth) parameters
		"Range_X" = 1:30,
		"Range_Y" = 1:30,
		"SD_O" = 100,						    # SD of spatial observation
		"SpatialScale" = 3,					# scale of spatial variation
		"year.depth.restriction" = 1000,	# If MPA, then years when MPA implemented (default > n_years)
		"depth.restriction" = c(50,150),	# location of MPA

	### Other OM control features
		"plotting" = FALSE,
		"Interactive" = FALSE
	)										# End of simulation setting


	Sim_Settings1$Fish_depth_par1 = c(120, 250, 150, 500)
	Data <- Generate_scenario_data(Sim_Settings1, seed=13)
	Sim_Settings2 <- Sim_Settings1
	Sim_Settings2$SD_O = 200
	Sim_Settings2$SpatialScale = 10
	Sim_Settings2$Fish_depth_par1 = c(120, 250, 150, 500) # bycatch species has similar distribution as another species
	Data <- Generate_scenario_data(Sim_Settings2, seed=32)
	Sim_Settings3 <- Sim_Settings2
	Sim_Settings3$Fish_depth_par1 = c(120, 250, 450, 200) # bycatch species is distributed deeper than the main target
	Data <- Generate_scenario_data(Sim_Settings3, seed=32)


	qwe <- apply(Data$Biomass, 2, function(x) x/x[1])
	matplot(qwe, type="l", ylim=c(0,1))


	library(tidyverse)
	sp_sel = 3

	asd <- Data$Data %>% as.data.frame() %>% subset(Species == sp_sel) %>% mutate(CPUE = ifelse(CPUE < 1, 0, CPUE))
	ggplot(asd, aes(x=CPUE)) + geom_histogram(bins=100) + theme_bw()
	table(asd$CPUE==0)/nrow(asd)


#### Sampling selection --> need to implement some impartial coverage issue
	## e.g. cost function with ports, divide the coast in areas and include fishing areas by vessels
	sampling_select <- function(data, percent, unit="vessel", seed = 123){
	  set.seed(123)
	  if (unit == "vessel") {
	    vessel_sel = sample(unique(data$vessel), floor(length(unique(data$vessel))*percent))
	    dat <- data %>% filter(vessel %in% vessel_sel)
	  }
	  if (unit == "fishing") {
	    obs_sel = sample(1:nrow(data), floor(nrow(data)*percent))
	    dat <- data[obs_sel,]
	  }
    return(dat)
	}

  Final_data <- sampling_select(data = asd, percent = 0.1, unit="vessel", seed=123)
  Final_data <- sampling_select(data = asd, percent = 0.1, unit="fishing", seed=123)

  plot_fishing <- function(data, years, ...){
    data %>% filter(year %in% years) %>% group_by(X,Y,year) %>% summarize(total_effort = n()) %>%
      ggplot(., aes(x=X, y=Y)) +
      facet_wrap(.~year)  + theme_bw() +
      geom_tile(data = Data$bathym, aes(x=X, y=Y, col=depth, fill=depth)) +
      scale_fill_viridis_c() + scale_color_viridis_c()+ geom_point()
  }
  plot_fishing(Final_data, years = c(1,5,10,15,20,25,30))

##### Now we need to perform the CPUE standardization

	#### Method 1: Applying sdmTMB to single species data
  library(sdmTMB)

  Final_data$vessel_fct <- as.factor(Final_data$vessel)
  mesh <- make_mesh(Final_data, xy_cols = c("X", "Y"), cutoff = 2)
  # plot(mesh)
  fit <- sdmTMB(CPUE ~ 0 + as.factor(year) + s(depth, k=3) , data = Final_data, time = "year", family = tweedie(link="log"), spde = mesh)
  fit
  # res_test <- residuals(fit)
  # qqnorm(res_test);qqline(res_test)

  qcs_grid <- do.call(rbind, replicate(30, data.bathym, simplify=FALSE))
  qcs_grid$year <- as.numeric(rep(1:30, each=nrow(data.bathym)))
  qcs_grid$vessel <- 8
  qcs_grid <-  qcs_grid %>% mutate(X = as.numeric(X), Y= as.numeric(Y))
  p_st <- predict(fit, newdata = qcs_grid,
                  return_tmb_object = TRUE)
  index <- get_index(p_st)
  ggplot(index, aes(year, est)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey90") +
    geom_line(lwd = 1, colour = "grey30") +
    labs(x = "Year", y = "Biomass (kg)") +
    geom_line(data = data.frame(est=qwe[,sp_sel]*index$est[1], year =1:30), col="red", size=2)+
    theme_bw()+coord_cartesian(ylim=c(0, 20000))


#### glmmTMB
  library(glmmTMB)
  lm1 <- glmmTMB(CPUE ~ 0 + as.factor(year) + depth + I(depth^2), data=Final_data, family=tweedie)
  predicted <- predict(lm1, newdata = qcs_grid, type="response", se.fit=F)
  pred <- qcs_grid
  pred$pred = predicted[-length(predicted)]
  IA <- pred %>% group_by(year) %>% summarize(IA = sum(pred))
  ggplot(IA, aes(x=year, y=IA)) + labs(x = "Year", y = "Index") +
    geom_line(lwd = 1, colour = "grey30") +
    geom_line(data = data.frame(IA=qwe[,sp_sel]*IA$IA[1], year =1:30), col="red", size=2)+
    theme_bw()


#### spatial DFA or VAST


#### VAST



#### Scenario 2: targetting changes
### year 1:10 	= base price
### year 11:20	= target petrale, avoid rockfish and sablefish, same for dover
### year 21:30	= back to "normal"


	Sim_Settings2 <- list(
	"n_years" = 30,

	### Species parameters
		"n_species"  = 4,														# petrale, dover, darkblotched, sablefish
		"price_fish" = matrix(c(rep(c(1,1,1,1),10),								# important because controls vessel dynamics
								rep(c(2,1,0.5,0.5),10),
								rep(c(1,1,1,1),10)), ncol=4, byrow=T),

		## Fish habitat preference/movement control paramters
			"func_mvt_dist" = c("Exp", "Exp", "Exp", "Exp"),					# shape of mvt p_function
			"func_mvt_depth" = c("Lognorm", "Lognorm", "Lognorm", "Lognorm"),	# shape of depth preference p_function (choice of Normal, Exp, Lognormal, Uniform)
			"func_mvt_lat" = c("Unif", "Unif", "Unif", "Unif"),					# shape of range preference p_function (choice of Normal, Exp, Lognormal, Uniform)

			"Fish_dist_par1" = c(rep(0.4, 3), 1),								# mvt distance
			"Fish_dist_par2" = rep(0, 4), 										# mvt distance
			"Fish_depth_par1" = c(220, 450, 250, 1200),         				# depth preference
			"Fish_depth_par2" = c(0.5, 0.5, 0.3, 0.5),        					# depth preference
			"Fish_range_par1" = rep(0, 4),                    					# X-axis range
			"Fish_range_par2" = rep(50,4),                    					# X axis range

		## Species Schaefer pop dyn parameters
			"B0" = c(50000,820000,36000,540000),								# petrale, dover, darkblotched, sablefish (mt)
			"r" = c(0.22, 0.17, 0.08, 0.08),
			"no.mvt" = FALSE,				# whether animals redistribute based on habitat preference

	### Parameters controlling the vessel dynamics
		"Nvessels" = 30,					# nb vessels
		"Tot_effort" = 200,					# end year nb of effort
		"CV_effort" = 0.2, 					# CV of effort around the yearly mean effort
		"qq_original" = 0.2,				# the average catchability coef for ALL vessels
		"CV_vessel" = 0.2,					# the between vessel variability in mean catchability
		"do.tweedie" = TRUE,				# include obsrevation error in catchability
		"xi" = 1.2,							    # controls observation error in vessel catchability
		"phi" = 0.5,						    # controls observation error in vessel catchability
		"Preference" = 1, 					# controls how effort concentrate into areas
		"Changing_preference" = FALSE, 		# whether effort concentration changes over time

	### Habitat (depth) parameters
		"Range_X" = 1:20,
		"Range_Y" = 1:20,
		"SD_O" = 150,						# SD of spatial observation
		"SpatialScale" = 2,					# scale of spatial variation
		"year.depth.restriction" = 1000,	# If MPA, then years when MPA implemented (default > n_years)
		"depth.restriction" = c(50,150),	# location of MPA

	### Other OM control features
		"plotting" = TRUE,
		"Interactive" = TRUE
	)										# End of simulation setting


#### Generate data
	Data <- Generate_scenario_data(Sim_Settings1, seed=123)
	Data <- Generate_scenario_data(Sim_Settings2)


	library(tidyverse)

	asd <- Data$Data %>% as.data.frame() %>% subset(Species == 3) %>% mutate(CPUE = ifelse(CPUE < 0.1, 0, CPUE))
	ggplot(asd, aes(x=CPUE)) + geom_histogram(bins=100) + theme_bw()
  table(asd$CPUE==0)/nrow(asd)

	Dat_df <- Data$Data %>% as.data.frame() %>% pivot_wider(., names_from=Species, values_from = CPUE)




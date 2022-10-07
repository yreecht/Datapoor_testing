
## Main function to generate data (the OM)
Generate_scenario_data <- function(Sim_Settings, seed_input = 123, parallel = FALSE){


	#### Set the bathymetry field
  ## SD_O sets the absolute scale of variability in depth e.g. from 0-100, from 0-500, etc
  ## scale sets the coarseness of spatial variation: the bigger, the more things are correlated in space.
  ## a small value e.g. 1 will make the map very patchy (each cell is independent of its neighboor)
	set.seed(seed_input)
	model_bathym <- RMgauss(var=Sim_Settings$SD_O^2, scale=Sim_Settings$SpatialScale)
	map_grid <- expand.grid(X=Sim_Settings$Range_X, Y=Sim_Settings$Range_Y)
	Bathym <- RFsimulate(model_bathym, x=map_grid)
	data.bathym <- data.frame(ID = 1:nrow(map_grid), map_grid , depth=Bathym@data[,1])
	data.bathym$depth = data.bathym$depth+abs(min(data.bathym$depth))	# to avoid getting negative depth here depth takes positive real number
	# image.plot(Sim_Settings$Range_X, Sim_Settings$Range_Y, matrix(data.bathym$depth,nrow=length(Sim_Settings$Range_X), ncol=length(Sim_Settings$Range_Y)))

	#### Set the pop mvt param from the simulation settings for each season (here month)
	Par_mvt_adult <- lapply(1:12, function(x)
	  rbind(Sim_Settings$Fish_dist_par1[x,], Sim_Settings$Fish_dist_par2, Sim_Settings$Fish_depth_par1[x,], Sim_Settings$Fish_depth_par2[x,]));

	#### Calculate the movement matrices
	if (Sim_Settings$parallel == FALSE) Mvt_mat_adult <- lapply(1:12, function(x)
	  lapply(1:Sim_Settings$n_species, function(xxx) Mvt_matrix(Pop_param=Par_mvt_adult[[x]][,xxx], Depth_eff="TRUE", Dist_eff="TRUE", Lat_eff="TRUE", Dist_func=Sim_Settings$func_mvt_dist[xxx], Depth_func=Sim_Settings$func_mvt_depth[xxx], Lat_func=Sim_Settings$func_mvt_lat[xxx], Eastings=Sim_Settings$Range_X, Northings=Sim_Settings$Range_Y, data.bathym)))
	if (Sim_Settings$parallel == TRUE) {
	  cl <- makeCluster(getOption("cl.cores", 4))
	  clusterExport(cl, varlist=list("Sim_Settings","Mvt_upd1","Mvt_matrix","Par_mvt_adult","data.bathym"))
	  Mvt_mat_adult <- parLapply(cl, 1:12, function(x)
	    lapply(1:Sim_Settings$n_species, function(xxx)
	      Mvt_matrix(Pop_param=Par_mvt_adult[[x]][,xxx], Depth_eff="TRUE", Dist_eff="TRUE", Lat_eff="TRUE", Dist_func=Sim_Settings$func_mvt_dist[xxx], Depth_func=Sim_Settings$func_mvt_depth[xxx], Lat_func=Sim_Settings$func_mvt_lat[xxx], Eastings=Sim_Settings$Range_X, Northings=Sim_Settings$Range_Y, data.bathym)))
	}

	#### Biomass dynamics: similarly to the work from Thorson et al but with more features
	### 1. Determining the initial population distribution
	Pop_adult <- sapply(1:Sim_Settings$n_species, function(x) Stable_pop_dist(Mvt_mat_adult[[1]][[x]], Sim_Settings$B0[x], Sim_Settings$Range_X, Sim_Settings$Range_Y))

	if(Sim_Settings$plotting==TRUE){
		X11()
		nf <- layout(matrix(1:6, nrow=3))
		par(mar=c(1,1,1,1))
		fields::image.plot(Sim_Settings$Range_X, Sim_Settings$Range_Y, matrix(data.bathym$depth,nrow=length(Sim_Settings$Range_X), ncol=length(Sim_Settings$Range_Y)))
		for (sp in 1:Sim_Settings$n_species) fields::image.plot(Sim_Settings$Range_X, Sim_Settings$Range_Y, matrix(Pop_adult[,sp],length(Sim_Settings$Range_X),length(Sim_Settings$Range_Y)))
	}

	Biomass <- array(NA, dim=c(Sim_Settings$n_years, 12, nrow(data.bathym), Sim_Settings$n_species))
	Biomass[1,1,,] <- t(apply(Pop_adult,1,function(x) x*Sim_Settings$Pop_ratio_start))
	Biomass[1,1,,] <- apply(Biomass[1,1,,], 2,function(x) replace(x, which(x<=0), 0))

	### 2. Structuring the vessel dynamics

	## generate total effort by year & month (total effort will be randomly sampled with log normal error around a logistic type function)
	Effort <- array(0, dim=c(Sim_Settings$n_years, 12));
	map.size <- length(Sim_Settings$Range_X)*length(Sim_Settings$Range_Y)
	Effort_area_year <- array(0, dim=c(Sim_Settings$n_years, 12, map.size));

	Sig_effort <- sqrt(log(Sim_Settings$CV_effort^2+1))
	Sig_vessel <- sqrt(log(Sim_Settings$CV_vessel^2+1))
	qq_vessel <- sapply(1:Sim_Settings$n_species, function(x) rlnorm(Sim_Settings$Nvessels, log(Sim_Settings$qq_original[x])-Sig_vessel^2/2, Sig_vessel))

	# use of tweedie distribution to generate catch data for individual boat
	CPUE_tot <- array(NA, dim=c(Sim_Settings$n_years, 12, nrow(data.bathym)))
	Catch_area_year <- array(0, dim=c(Sim_Settings$n_years, 12, map.size, Sim_Settings$n_species));
	Catch_area_year_ind <- c();

	# Define historic fishing regions for each vessel randomly
	which_regions <- t(sapply(rep(sum(Sim_Settings$Nregion),Sim_Settings$Nvessels), function(x) sample(1:sum(Sim_Settings$Nregion), x, replace=TRUE)))
	areas <- equal_partition( expand.grid(X = Sim_Settings$Range_X, Y = Sim_Settings$Range_Y), Sim_Settings$Nregion[1], Sim_Settings$Nregion[2], Sim_Settings$Range_X, Sim_Settings$Range_Y)$area

	# Keep track of the population abundance/number that left the fishing ground
# 	Gone_outside <- array(0, dim=c(Sim_Settings$n_years, 12, 1, Sim_Settings$n_species))
#   Gone_total <- rep(0, Sim_Settings$n_species)  ## a dummy variable to keep track of the total population that left

	#### Updating the population and generate catch data according to the above specifications
	for (iyear in 1:Sim_Settings$n_years){
	  for (month in 1:12){

	    set.seed(seed_input + iyear*12 + month)
  		### do vessel change their fishing preference over the year (this ONLY controls vessel concentration factor)
	    Preference=1
	    if(Sim_Settings$Changing_preference==TRUE) Preference <- 2-exp(-0.1*iyear)

  		### Assigning total effort in a year
	    Effort[iyear,] <- Sim_Settings$Tot_effort_year[iyear]*Sim_Settings$Effort_alloc_month/sum(Sim_Settings$Effort_alloc_month)		# mean effort by year

  		### Adjusting the catchability by month (if any deviation in catchability due to changes in the availability of a population)
  		qq <- Sim_Settings$qq_original*Sim_Settings$Rangeshift_catchability_adjust[month,]

  		### Continuous catch equation - calculating the expected revenue
  		CPUE_tot[iyear,month,] <- apply(Biomass[iyear,month,,]*matrix((1-exp(-qq)), nrow=map.size, ncol=Sim_Settings$n_species, byrow=T) * matrix(Sim_Settings$price_fish[iyear,], nrow=map.size, ncol=Sim_Settings$n_species, byrow=T),1,sum)

  		### If a site is not worth going (expected revenue is negative)
  		CPUE_tot[iyear,month,][which(CPUE_tot[iyear,month,] < 0 )]<- 0

  		### Distribute the total effort according to the above expected revenue
  		p <- CPUE_tot[iyear,month,]^Preference/sum(CPUE_tot[iyear,month,]^Preference);min(p);max(p)
  		Effort_area_year[iyear,month,] <- rmultinom(1, size=Effort[iyear,], p=p)

  		### Then calculate the total catch in each area
  		aa <- c(); AA <- c()
  		Catch_year_area_mat <- matrix(0, nrow=map.size, ncol=Sim_Settings$n_species)

  		vessel_seed <- sample(1:Sim_Settings$vessel_seeds, 1)
  		set.seed(vessel_seed)     # this piece of code is to create heterogeneity in the vessel fishing times (without it all vessels fish about the same number of times)
  		for (xyz in which(Effort_area_year[iyear,month,]>0)){
  			vessel_pool_area <- (1:Sim_Settings$Nvessels)[apply(which_regions, 1, function(x) areas[xyz] %in% x)]
  		  which.vessel.fished <- sample(vessel_pool_area, Effort_area_year[iyear,month,xyz], replace=TRUE)

  			# This is the code using the generating random obs error per vessel & species
  			if(Sim_Settings$do.tweedie==TRUE) rand.vessel <- matrix(sapply(1:Sim_Settings$n_species, function(x) rtweedie( Effort_area_year[iyear,month,xyz], mu= qq_vessel[which.vessel.fished,x], power=Sim_Settings$xi[x], phi=Sim_Settings$phi[x])),ncol=Sim_Settings$n_species)
  		  # if(Sim_Settings$do.tweedie==TRUE) rand.vessel <- matrix(sapply(1:Sim_Settings$n_species, function(x) rtweedie_new(n=Effort_area_year[iyear,month,xyz], mu= qq_vessel[which.vessel.fished,x], power=Sim_Settings$xi[x], phi=Sim_Settings$phi[x])),ncol=Sim_Settings$n_species)
  		  if(Sim_Settings$do.tweedie==FALSE) rand.vessel <- matrix(sapply(1:Sim_Settings$n_species, function(x) qq_vessel[which.vessel.fished,x]),ncol=Sim_Settings$n_species)
  			catch_area_vessel <- sapply(1:Sim_Settings$n_species, function(Y) { hum <- Biomass[iyear,month,xyz,Y]*(1-exp(-sum(rand.vessel[,Y]))); ifelse(hum==0, catch_area_vessel <- rep(0, Effort_area_year[iyear,month,xyz]), catch_area_vessel <- hum*rand.vessel[,Y]/sum(rand.vessel[,Y])); return(catch_area_vessel)})
  			# now truncate value when necessary
  			if (any(Sim_Settings$catch_trunc == 1)) {
  			  to_trunc <- which(Sim_Settings$catch_trunc == 1)
  			  if (length(catch_area_vessel) == Sim_Settings$n_species)  catch_area_vessel[to_trunc ] <- trunc(catch_area_vessel[to_trunc ])
  			  if (length(catch_area_vessel) > Sim_Settings$n_species)  catch_area_vessel[,to_trunc ] <- trunc(catch_area_vessel[,to_trunc ])
  			}
  			catch_area_vessel <- data.frame(iyear, month, areas[xyz], matrix(rep(c(data.bathym[xyz,-1]), each= Effort_area_year[iyear,month,xyz]),nrow= Effort_area_year[iyear,month,xyz]), which.vessel.fished, matrix(catch_area_vessel,ncol=Sim_Settings$n_species));
  			catch_area_vessel <- matrix(sapply( catch_area_vessel, FUN=as.numeric ), nrow=nrow(catch_area_vessel))
  			colnames(catch_area_vessel) <- c("year", "month", "fishing_district", "X", "Y", "depth", "vessel", paste0("Sp", 1:Sim_Settings$n_species))
  			Catch_year_area_mat[xyz,] <- colSums(catch_area_vessel[,grep("Sp", colnames(catch_area_vessel)),drop=FALSE])

  			AA <- rbind(AA, catch_area_vessel)
  		}

  		Catch_area_year_ind <- rbind(Catch_area_year_ind, AA)
  		Catch_area_year_ind <- Catch_area_year_ind %>% as.data.frame()
  		### Now dividing the above catch based on the discard rate to "landed" (kept) and "discarded"
  		### then the survival rate determines how much of the discard rejoins the population

  		  ## Including the discard rate
    		Catch_area_year_ind_discard <- Catch_area_year_ind %>% as.data.frame() %>%
    		    dplyr::select(starts_with("Sp"))
    		for (sp in 1:Sim_Settings$n_species){
    		  Mean = Sim_Settings$Discard_rate_beta1[sp]
    		  Var = Sim_Settings$Discard_rate_beta2[sp]
    		  shape1 = (Var - Mean*(1-Mean))/(Mean * Var)
    		  shape2 = (Var - Mean*(1-Mean))/(Mean * Var)*(1-Mean)
    		  Catch_area_year_ind_discard[,sp] <- Catch_area_year_ind_discard[,sp]*ifelse(Var==0, Mean, rbeta(nrow(Catch_area_year_ind_discard),shape1, shape2))
    		  if (Sim_Settings$catch_trunc[sp] == 1) Catch_area_year_ind_discard[,sp] <- trunc(Catch_area_year_ind_discard[,sp])
    		}
    		colnames(Catch_area_year_ind_discard) <- paste0("Sp", 1:Sim_Settings$n_species, "_disc")

    		Catch_area_year_ALL <- cbind(Catch_area_year_ind, Catch_area_year_ind_discard)

  		  ## Including the survival rate to calculate amount (either mass or number) that survived
    		Catch_area_year_ind_survived <- Catch_area_year_ind_discard
    		for (sp in 1:Sim_Settings$n_species){
    		  Mean_surv = Sim_Settings$Discard_survival_rate1[sp]
    		  Var_surv = Sim_Settings$Discard_survival_rate2[sp]
    		  shape1_surv = (Var_surv - Mean_surv*(1-Mean_surv))/(Mean_surv * Var_surv)
    		  shape2_surv = (Var_surv - Mean_surv*(1-Mean_surv))/(Mean_surv * Var_surv)*(1-Mean_surv)
    		  Catch_area_year_ind_survived[,sp] <- Catch_area_year_ind_discard[,sp]*ifelse(Var_surv==0, Mean_surv, rbeta(nrow(Catch_area_year_ind_discard),shape1_surv, shape2_surv))
    		  if (Sim_Settings$catch_trunc[sp] == 1) Catch_area_year_ind_survived[,sp] <- trunc(Catch_area_year_ind_survived[,sp])
    		}

    		## Putting back these surviving individual to the sea = biomass
    		for (sp in 1:Sim_Settings$n_species){
    		  where_survived <- which(Catch_area_year_ind_survived[,sp]>0)
    		  if (length(where_survived) > 0) {
    		    for (iii in where_survived) {
    		      grid_cell <- Convert_grid_coordinates(X=Catch_area_year_ind$X[iii], Y=Catch_area_year_ind$Y[iii], Settings=Sim_Settings)
    		      Biomass[iyear,month,grid_cell,sp] <- Biomass[iyear,month,grid_cell,sp] + Catch_area_year_ind_survived[iii,sp]
    		    }
    		  }
    		}


  		Catch_area_year[iyear,month,,] <- Catch_year_area_mat

     	### Then update the population = growth & catch, then at a specific month, we have the process error (e.g. rec dev)
  		if (month < 12){
    		temp <- (Biomass[iyear,month,,]+Sim_Settings$r*Biomass[iyear,month,,]*(1-Biomass[iyear,month,,]/Biomass[1,1,,])-Catch_area_year[iyear,month,,])
    		temp <- apply(temp,2,function(x) replace(x, which(is.na(x)==TRUE | x<=0),0))
    		temp <- sapply(1:Sim_Settings$n_species, function(x) {if (month %in% Sim_Settings$sigma_p_timing[x]) out <- sapply(seq_along(temp[,x]), function(y) rlnorm(1, log(temp[y,x])- Sim_Settings$sigma_p[x]^2/2, Sim_Settings$sigma_p[x])) else temp[,x]})
    		Biomass[iyear,month+1,,] <- apply(temp,2,function(x) replace(x, which(is.na(x)==TRUE | x<=0),0))
  		}
  		if (month == 12 & iyear < Sim_Settings$n_years){
    		temp <- Biomass[iyear,month,,]+Sim_Settings$r*Biomass[iyear,month,,]*(1-Biomass[iyear,month,,]/Biomass[1,1,,])-Catch_area_year[iyear,month,,]
    		temp <- apply(temp,2,function(x) replace(x, which(is.na(x)==TRUE | x<=0),0))
    		temp <- sapply(1:Sim_Settings$n_species, function(x) {if (month %in% Sim_Settings$sigma_p_timing[x]) out <- sapply(seq_along(temp[,x]), function(y) rlnorm(1, log(temp[y,x])- Sim_Settings$sigma_p[x]^2/2, Sim_Settings$sigma_p[x])) else temp[,x]})
    		Biomass[iyear+1,1,,] <- apply(temp,2,function(x) replace(x, which(is.na(x)==TRUE | x<=0),0))
  		}

  		### population mvt at the end of each month (to their favorite habitat)
  		if (month < 12){
  		  if (Sim_Settings$fish.mvt==TRUE)	Biomass[iyear,month+1,,] <- sapply(1:Sim_Settings$n_species, function(x) t(Mvt_mat_adult[[month]][[x]])%*%Biomass[iyear,month+1,,x])
  		  if (Sim_Settings$fish.mvt==FALSE)	Biomass[iyear,month+1,,] <- Biomass[iyear,month+1,,]
  		  Biomass[iyear,month+1,,] <- apply(Biomass[iyear,month+1,,],2,function(x) replace(x, which(is.na(x)==TRUE | x<=0),0))
  		}
  		if (month == 12 & iyear < Sim_Settings$n_years){
  		  if (Sim_Settings$fish.mvt==TRUE)	Biomass[iyear+1,1,,] <- sapply(1:Sim_Settings$n_species, function(x) t(Mvt_mat_adult[[month]][[x]])%*%Biomass[iyear+1,1,,x])
  		  if (Sim_Settings$fish.mvt==FALSE)	Biomass[iyear+1,1,,] <- Biomass[iyear+1,1,,]
  		  Biomass[iyear+1,1,,] <- apply(Biomass[iyear+1,1,,],2,function(x) replace(x, which(is.na(x)==TRUE | x<=0),0))
  		}

	  }
	  if (Sim_Settings$Interactive == TRUE) print(iyear)

	 }

	#### Return data

	Return = list(Data= Catch_area_year_ALL, Biomass=Biomass, bathym=data.bathym)
	return(Return)
}

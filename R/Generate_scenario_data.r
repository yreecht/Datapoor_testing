
## Main function to generate data (the OM)
Generate_scenario_data <- function(Sim_Settings, seed = 123, parallel = FALSE){

  ## generate mvt matrix of an animal from a cell to the other
	Mvt_matrix <- function(Pop_param, Depth_eff, Dist_eff, Lat_eff, Depth_func, Dist_func="Exp", Lat_func, x,y,data.bathym){
		Mvt_mat <- matrix(0,(length(x)*length(y)),(length(x)*length(y)))
		for (loc in 1:(length(x)*length(y))){
			Mvt_mat[loc,] <- Mvt_upd1(loc, Pop_param, Depth_eff=Depth_eff, Dist_eff=Dist_eff, Lat_eff=Lat_eff, Depth_func=Depth_func, Dist_func=Dist_func, Lat_func=Lat_func,x,y, data.bathym)
		}
		return(Mvt_mat)
	}

  ## for updating the population distribution according to different different covariates (depth and/or dist), using joint prob
	Mvt_upd1 <- function(loc, Pop_param, Depth_eff, Dist_eff, Lat_eff, Depth_func, Dist_func="Exp", Lat_func, x,y, data.bathym){
		Di1 <- Pop_param[1]  	# 1st value for dist function
		Di2 <- Pop_param[2]  	# 2nd value for dist function
		De1 <- Pop_param[3]  	# 1st value for depth function
		De2 <- Pop_param[4]	 	# 2nd value for depth function
		Lat1 <- Pop_param[5]  # 1st value for lat/long function
		Lat2 <- Pop_param[6]	# 2nd value for lat/long function

		Dist <- sqrt((data.bathym$X-data.bathym$X[loc])^2+(data.bathym$Y-data.bathym$Y[loc])^2)

		Prob1 <- 1
		if(Dist_eff=="T"){
			Prob1 <- dexp(Dist, rate=1/Di1)
		}

		Prob2 <- 1
		if(Depth_eff=="T"){
			if(Depth_func=="Norm") Prob2 <- dnorm(data.bathym$depth, De1, De2*De1)
			if(Depth_func=="Exp") Prob2 <- dexp(data.bathym$depth, rate=1/De1)
			if(Depth_func=="Lognorm") { SD <- sqrt(log(De2^2+1)); Prob2 <- dlnorm(data.bathym$depth, (log(De1)-SD^2/2), SD) }
			if(Depth_func=="Unif") Prob2 <- dunif(data.bathym$depth, De1, De2)
		}

		Prob3 <- 1
		if(Lat_eff=="T"){
			if(Lat_func=="Norm") Prob3 <- dnorm(data.bathym$X, Lat1, Lat2)
			if(Lat_func=="Exp") Prob3 <- dexp(data.bathym$X, rate=1/Lat1)
			if(Lat_func=="Lognorm") Prob3 <- dlnorm(data.bathym$X, (log(Lat1)-Lat2^2/2), Lat2)
			if(Lat_func=="Unif") Prob3 <- dunif(data.bathym$X, Lat1, Lat2)
		}

		Prob <- Prob1*Prob2*Prob3/sum(Prob1*Prob2*Prob3)
		return(Prob)
	}

  ## Calculate the initial population distribution (start from uniform distribution in space)
	Stable_pop_dist <- function(Mvt_mat_adult, B0){
		Pop1 <- rep(B0/(length(Sim_Settings$Range_X)*length(Sim_Settings$Range_Y)), (length(Sim_Settings$Range_X)*length(Sim_Settings$Range_Y)))
		Mvt_mat_adult1 <- t(Mvt_mat_adult)
		for (Time in 1:1000){
			Pop <- Mvt_mat_adult1%*%Pop1
			Pop1 <- Pop
		}
		Initpop_adult_mat <- matrix(Pop1, nrow=length(Sim_Settings$Range_X), ncol=length(Sim_Settings$Range_Y))
		return(Initpop_adult_mat)
	}


	#### Set the bathymetry field
	set.seed(seed)
	model_bathym <- RMgauss(var=Sim_Settings$SD_O^2, scale=Sim_Settings$SpatialScale)
	map_grid <- expand.grid(X=Sim_Settings$Range_X, Y=Sim_Settings$Range_Y)
	Bathym <- RFsimulate(model_bathym, x=map_grid)
	data.bathym <- data.frame(ID = 1:nrow(map_grid), map_grid , depth=Bathym@data[,1])
	data.bathym$depth = data.bathym$depth+abs(min(data.bathym$depth))	# to avoid getting negative depth here depth takes positive real number
	image.plot(Sim_Settings$Range_X, Sim_Settings$Range_Y, matrix(data.bathym$depth,nrow=length(Sim_Settings$Range_X), ncol=length(Sim_Settings$Range_Y)))

	#### Set the pop mvt param
	Par_mvt_adult <- lapply(1:12, function(x)
	  rbind(Sim_Settings$Fish_dist_par1, Sim_Settings$Fish_dist_par2, Sim_Settings$Fish_depth_par1[x,], Sim_Settings$Fish_depth_par2[x,], Sim_Settings$Fish_range_par1, Sim_Settings$Fish_range_par2));
	if (Sim_Settings$parallel == FALSE) Mvt_mat_adult <- lapply(1:12, function(x)
	  lapply(1:Sim_Settings$n_species, function(xxx) Mvt_matrix(Pop_param=Par_mvt_adult[[x]][,xxx], Depth_eff="T", Dist_eff="T", Lat_eff="T", Dist_func=Sim_Settings$func_mvt_dist[xxx], Depth_func=Sim_Settings$func_mvt_depth[xxx], Lat_func=Sim_Settings$func_mvt_lat[xxx], x=Sim_Settings$Range_X, y=Sim_Settings$Range_Y, data.bathym)))
	if (Sim_Settings$parallel == TRUE) {
	  cl <- makeCluster(getOption("cl.cores", 4))
	  clusterExport(cl, varlist=list("Sim_Settings","Mvt_upd1","Mvt_matrix","Par_mvt_adult","data.bathym"))
	  Mvt_mat_adult <- parLapply(cl, 1:12, function(x)
	    lapply(1:Sim_Settings$n_species, function(xxx)
	      Mvt_matrix(Pop_param=Par_mvt_adult[[x]][,xxx], Depth_eff="T", Dist_eff="T", Lat_eff="T", Dist_func=Sim_Settings$func_mvt_dist[xxx], Depth_func=Sim_Settings$func_mvt_depth[xxx], Lat_func=Sim_Settings$func_mvt_lat[xxx], x=Sim_Settings$Range_X, y=Sim_Settings$Range_Y, data.bathym)))
	}

	#### Biomass dynamics: similarly to the work from Caruthers, we assume that there is a regional carrying capacity
	# Initial population distribution
	Pop_adult <- sapply(1:Sim_Settings$n_species, function(x) Stable_pop_dist(Mvt_mat_adult[[1]][[x]], Sim_Settings$B0[x]))

	if(Sim_Settings$plotting==TRUE){
		windows()
		nf <- layout(matrix(1:6, nrow=3))
		par(mar=c(1,1,1,1))
		fields::image.plot(Sim_Settings$Range_X, Sim_Settings$Range_Y, matrix(data.bathym$depth,nrow=length(Sim_Settings$Range_X), ncol=length(Sim_Settings$Range_Y)))
		fields::image.plot(Sim_Settings$Range_X, Sim_Settings$Range_Y, matrix(Pop_adult[,1],length(Sim_Settings$Range_X),length(Sim_Settings$Range_Y)))
		fields::image.plot(Sim_Settings$Range_X, Sim_Settings$Range_Y, matrix(Pop_adult[,2],length(Sim_Settings$Range_X),length(Sim_Settings$Range_Y)))
		fields::image.plot(Sim_Settings$Range_X, Sim_Settings$Range_Y, matrix(Pop_adult[,3],length(Sim_Settings$Range_X),length(Sim_Settings$Range_Y)))
		fields::image.plot(Sim_Settings$Range_X, Sim_Settings$Range_Y, matrix(Pop_adult[,4],length(Sim_Settings$Range_X),length(Sim_Settings$Range_Y)))
		fields::image.plot(Sim_Settings$Range_X, Sim_Settings$Range_Y, matrix(Pop_adult[,5],length(Sim_Settings$Range_X),length(Sim_Settings$Range_Y)))
	}

	Biomass <- array(NA, dim=c(Sim_Settings$n_years, 12, nrow(data.bathym), Sim_Settings$n_species))
	Biomass[1,1,,] <- Pop_adult
	Biomass[1,1,,] <- apply(Biomass[1,1,,], 2,function(x) replace(x, which(x<=0), 0))

	#### Vessel dynamics

	### Determining the predicted area specific CPUE
	map.size <- length(Sim_Settings$Range_X)*length(Sim_Settings$Range_Y)

	# generate total effort by year & month (total effort is randomly sampled with log normal error around a logistic type function)
	Effort <- array(0, dim=c(Sim_Settings$n_years, 12));
	Mean_Effort <- array(0, dim=c(Sim_Settings$n_years, 12));
	Effort_area_year <- array(0, dim=c(Sim_Settings$n_years, 12, map.size));

	Sig_effort <- sqrt(log(Sim_Settings$CV_effort^2+1))
	Sig_vessel <- sqrt(log(Sim_Settings$CV_vessel^2+1))
	qq_vessel <- sapply(1:Sim_Settings$n_species, function(x) rlnorm(Sim_Settings$Nvessels, log(Sim_Settings$qq_original[x])-Sig_vessel^2/2, Sig_vessel))


	# use of tweedie distribution to generate catch data for individual boat
	CPUE_tot <- array(NA, dim=c(Sim_Settings$n_years, 12, nrow(data.bathym)))
	Catch_area_year <- array(0, dim=c(Sim_Settings$n_years, 12, map.size, Sim_Settings$n_species));
	Catch_area_year_ind <- c();

	# Define historic fishing regions for each vessel randomly
	# nb_regions <- sample(c(1:sum(Sim_Settings$Nregion)), Sim_Settings$Nvessels, replace=T)
	# which_regions <- sapply(nb_regions, function(x) sample(1:Sim_Settings$Nregion, x, replace=F))
	which_regions <- t(sapply(rep(sum(Sim_Settings$Nregion),Sim_Settings$Nvessels), function(x) sample(1:sum(Sim_Settings$Nregion), x, replace=TRUE)))
	# which_regions <- do.call(rbind, which_regions)
	areas <- equal_partition( expand.grid(X = Sim2$Range_X, Y = Sim2$Range_Y), Sim_Settings$Nregion[1], Sim_Settings$Nregion[2])$area

	#### Updating the population and generate catch data according to the above specifications
	for (iyear in 1:Sim_Settings$n_years){
	  for (month in 1:12){

	    set.seed(seed + iyear*12 + month)
  		### do vessel change their fishing preference over the year (this ONLY controls vessel concentration factor)
	    Preference=1
	    if(Sim_Settings$Changing_preference==TRUE) Preference <- 2-exp(-0.1*iyear)

  		### Assigning total effort in a year
  		Mean_Effort[iyear,] <- (Sim_Settings$Tot_effort/(1+exp(-0.1*iyear)))*Sim_Settings$Effort_alloc_month/sum(Sim_Settings$Effort_alloc_month)		# mean effort by year
  		Effort[iyear,] <- trunc(rlnorm(12, log(Mean_Effort[iyear,]-Sig_effort^2/2), Sig_effort))		# effort by year

  		### changing the catchability "qq" every year depending on whether a specific area is closed or not
  		if (iyear>=Sim_Settings$year.depth.restriction){
  			catchability <- rep(1,map.size)
  			if(length(depth.restriction)==2) catchability <- replace(catchability, which(bathym >= Sim_Settings$depth.restriction[1] & bathym <= Sim_Settings$depth.restriction[2]), 0)
  			if(depth.restriction== "random") catchability <- replace(catchability, sample(1:map.size, size=20), 0)
  			qq <- catchability
  		} else { qq <- Sim_Settings$qq_original }

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
  		set.seed(vessel_seed)
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

  		if (Sim_Settings$Interactive == TRUE) print(iyear)
	  }
	 }

	#### Plotting
  if(Sim_Settings$plotting==TRUE) {
  	# plot total catch 1st year and last year
		# windows()
		nf <- layout(1:2)
		par(mar=c(1,1,2,1), oma=c(4,4,2,2))
		hist(Catch_area_year_ind[,'Sp1'], breaks=100)
		hist(Catch_area_year_ind[,'Sp2'], breaks=100)

	  # plot biomass change over time
		# windows()
		for (iyear in 1:Sim_Settings$n_years){
			nf <- layout(1:5)
			par(mar=c(1,1,2,1), oma=c(4,4,2,2))
			fields::image.plot(matrix(Biomass[iyear,1,,1],length(Sim_Settings$Range_X),length(Sim_Settings$Range_X)), main=paste("Year", iyear))
			fields::image.plot(matrix(Biomass[iyear,1,,2],length(Sim_Settings$Range_X),length(Sim_Settings$Range_Y)))
			fields::image.plot(matrix(Biomass[iyear,1,,3],length(Sim_Settings$Range_X),length(Sim_Settings$Range_Y)))
			fields::image.plot(matrix(Biomass[iyear,1,,4],length(Sim_Settings$Range_X),length(Sim_Settings$Range_Y)))
			fields::image.plot(matrix(Effort_area_year[[iyear]] , length(Sim_Settings$Range_X),length(Sim_Settings$Range_X)))
			if(Sim_Settings$Interactive==TRUE) gtools::ask()
		}

	  # Time series of effort, raw CPUE, depletion
		# windows()
		# nf <- layout(1:4)
		# par(mar=c(4,4,1,1))
		# # Effort
		# matplot(Effort, ylim=c(min(Effort*0.8), max(Effort*1.2)), lwd=2, type="l", xlab="Year")
		# lines(Mean_Effort, lwd=2)
		# legend("topleft", legend=c("Sp1","Sp2","Sp3","Sp4"), col=c("black","red","green","blue"), lty=1, bty="n")
		# # Total catch
		# aa <- t(sapply(1:n_years, function(x) apply(Catch_area_year[[x]],2,sum)))   # total catch
		# matplot(aa, type="l", ylab="Total catch", lty=1, xlab="Year", lwd=2)
		# # Raw CPUE (dash line) and biomass (solid line)
		# CPUE_raw <- aa/Effort
		# bb <- t(sapply(1:n_years, function(x) apply(Biomass[[x]],2,sum)/B0))
		# matplot(bb, type="l",ylim=c(0,1), lty=1, lwd=2, ylab="Relative level") 	# the true biomass trajectory in solid line
		# RAW_CPUE <- cbind(CPUE_raw[,1]/CPUE_raw[1,1], CPUE_raw[,2]/CPUE_raw[1,2], CPUE_raw[,3]/CPUE_raw[1,3], CPUE_raw[,4]/CPUE_raw[1,4])
		# matlines(RAW_CPUE, type="l",ylim=c(0,1), lty=2, lwd=2) 	# the raw CPUE trajectory in dash line
		# legend("topright", c("Biomass", "Raw CPUE"), lty=c(1,2), bty="n")
		# # Relationship raw CPUE vs depletion
		# plot(bb[,1], CPUE_raw[,1]/(CPUE_raw[1,1]), type="l", lwd=2, col=1, xlab="Biomass depletion level", ylab="Raw standardized CPUE", xlim=c(0,1), ylim=c(0,1), )
		# lines(bb[,2], CPUE_raw[,2]/CPUE_raw[1,2], type="l", lwd=2, col="red")
		# lines(bb[,3], CPUE_raw[,3]/CPUE_raw[1,3], type="l", lwd=2, col="green")
		# lines(bb[,4], CPUE_raw[,4]/CPUE_raw[1,4], type="l", lwd=2, col="blue")
		# abline(0,1, lty=2)
	}

	#### Return data                                                                # v.names="CPUE",

	Return = list(Data= Catch_area_year_ind, Biomass=Biomass, bathym=data.bathym)
	return(Return)
}

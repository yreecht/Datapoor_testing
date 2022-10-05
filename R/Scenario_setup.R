#### Setting OM parameters

### Scenario 1:
Sim1 <- list(

  ### Number of simulation years (time step is a month)
  n_years = 15,

  ### Habitat (depth) parameters. Very important as it controls how species are distributed in space, thus affects both
  ### the population dynamics and the vessel dynamics
  Range_X = 1:40,           # the x-axis extent. the bigger, the more complex topography you can get at the cost of simulation time
  Range_Y = 1:40,           # the y-axis extent. the bigger, the more complex topography you can get at the cost of simulation time
  SD_O = 100,						    # SD of depth observation in space. the bigger, the more variable is the depth
  SpatialScale = 3,					# the spatial correlation range. The bigger, the more correlated are depth in space (patchy if small). N.B: Absolute value of high/low depends on the extent of x-axis and y-axis
  year.depth.restriction = 1000,	# If MPA, then years when MPA implemented (default > n_years)
  depth.restriction = c(50,150),	# location of MPA

  ### Species parameters
  n_species  = 4,
  # price_fish = matrix(rep(c(1.5,1,-1,2,0),15), ncol=5, byrow=T),   # price per species across time (roughly in Us west coast)
  price_fish = matrix(rep(c(1,1,-10,0),15), ncol=4, byrow=T),   # random

  ## Fish habitat preference/movement control parameters
  func_mvt_dist = c("Exp", "Exp", "Exp", "Exp"),					# shape of mvt p_function
  func_mvt_depth = c("Lognorm", "Lognorm", "Lognorm", "Lognorm"),	# shape of depth preference p_function (choice of Normal, Exp, Lognormal, Uniform)
  func_mvt_lat = c("Unif", "Unif", "Unif", "Unif"),					# shape of range preference p_function (choice of Normal, Exp, Lognormal, Uniform)

  Fish_dist_par1 = matrix(c(3, 3, 3, 10, # mvt distance mean - their mobility within month1 only values >0
                            3, 3, 3, 10, # month2   only values >0
                            3, 3, 3, 10, # month3   only values >0
                            3, 3, 3, 10, # month4   only values >0
                            3, 3, 3, 10, # month5   only values >0
                            3, 3, 3, 10, # month6   only values >0
                            3, 3, 3, 10, # month7   only values >0
                            3, 3, 3, 10, # month8   only values >0
                            3, 3, 3, 10, # month9   only values >0
                            3, 3, 3, 10, # month10  only values >0
                            3, 3, 3, 10, # month11  only values >0
                            3, 3, 3, 10),# month12  only values >0
                            nrow=12, ncol=4, byrow=T),
  Fish_dist_par2 = rep(0, 4), 										# mvt distance (not used)
  Fish_depth_par1 = matrix(c(250, 500, 300, 150, # depth preference mean per month 1
                             250, 500, 300, 150, # month2
                             220, 500, 300, 150, # month3
                             220, 500, 300, 150, # month4
                             220, 500, 300, 150, # month5
                             220, 500, 300, 150, # month6
                             220, 500, 300, 150, # month7
                             220, 500, 300, 150, # month8
                             220, 500, 300, 150, # month9
                             220, 500, 300, 150, # month10
                             220, 500, 300, 150, # month11
                             250, 500, 300, 150),# month12
                           nrow=12, ncol=4, byrow=T),
  Fish_depth_par2 = matrix(c(1, 0.5, 0.5, 1, # depth preference sd log scale per month 1
                             1, 0.5, 0.5, 1, # month2
                             1, 0.5, 0.5, 1, # month3
                             1, 0.5, 0.5, 1, # month4
                             1, 0.5, 0.5, 1, # month5
                             1, 0.5, 0.5, 1, # month6
                             1, 0.5, 0.5, 1, # month7
                             1, 0.5, 0.5, 1, # month8
                             1, 0.5, 0.5, 1, # month9
                             1, 0.5, 0.5, 1, # month10
                             1, 0.5, 0.5, 1, # month11
                             1, 0.5, 0.5, 1),# month12
                           nrow=12, ncol=4, byrow=T),
  Rangeshift_proportion = matrix(c(0, 0, 0, 0, # how much (proportion) of the population leave the fishing area by season. 0=nothing, 1=all
                                   0, 0, 0, 0, # month2
                                   0, 0, 0, 0, # month3
                                   0, 0, 0, 0, # month4
                                   0, 0, 0, 0, # month5
                                   0, 0, 0, 0.4, # month6
                                   0, 0, 0, 0.9, # month7
                                   0, 0, 0, 0.9, # month8
                                   0, 0, 0, 0.4, # month9
                                   0, 0, 0, 0, # month10
                                   0, 0, 0, 0, # month11
                                   0, 0, 0, 0),# month12          # the proportion of the population that leave the fishing groun X axis range max
                                   nrow=12, ncol=4, byrow=T),                    #
  Rangeshift_distance = c(0,0,0,50),   # 0 if the species does not perform any range shift
  ## Species Schaefer pop dyn parameters
  B0 = c(10000,10000,10000,10000)*1000,	# based on petrale, dover, darkblotched, sablefish, maybe bird/marine mammal
  r = c(0.22, 0.17, 0.08, 0.15)/12,         # based on petrale, dover, darkblotched, sablefish (mt)
  sigma_p= c(0.4, 0.4, 0.4, 0.4),            # log-normal process error = recruitment to fishery
  sigma_p_timing= c(9, 9, 9, 9),               # when does recruitment happen?
  fish.mvt = TRUE,				                        # whether animals redistribute annually based on habitat preference

  ### Parameters controlling the vessel dynamics
  Nregion = c(2,2),         # nb of fishing regions for the fishermen (equal cut of grid along X and Y)
  Nvessels = 50,					  # nb vessels
  Tot_effort = 2000,				# end year nb of effort
  CV_effort = 0.2, 					# CV of effort around the yearly mean effort
  qq_original = c(0.05,0.05,0.05,1e-5),			  # the average catchability coef by species for ALL vessels (does not mean anything. just a scaler)
  catch_trunc = c(0,0,0,1),			  # truncate value if 1 otherwise keep continuous
  CV_vessel = 0.0001,					# the between vessel variability in mean catchability
  vessel_seeds = 10,        # this creates heterogeneity in the sampled vessels
  Effort_alloc_month = c(5,4,3,2,1,1,1,2,3,4,5,5),  # it is rescaled to sum to 1
  do.tweedie = TRUE,				# include observation error in vessel catchability
  # xi = c(1.7,1.7,1.2,1.5),  # power of the tweedie distribution. reduce this to have more 0s. ==0 or >0
  xi = c(1.9,1.9,1.9,1.9),  # power of the tweedie distribution. reduce this to have more 0s. ==0 or >0
  # phi = c(0.2,0.2,0.2,0.2), # the scaler of the tweedie distribution
  phi = c(0.5,0.5,0.5,0.5), # the scaler of the tweedie distribution
  Preference = 1, 					# controls how effort concentrate into areas. The higher, the more concentrated is the effort
  Changing_preference = FALSE, 		# whether effort concentration changes over time


  ### Parameter controlling the resampling procedure from the whole fleet
  samp_prob = 0.2,             # this is the sampling probability
  samp_unit = "vessel",        # the samping unit: "vessel" or "fishing" (i.e. random sample from all fishing events)
  samp_seed = 123,             # the seed for reproducibility of the samples taken
  samp_mincutoff = 0,        # cut-off value to round values to 0
  start_year = 5,

  ### Other OM control features
  plotting = FALSE,
  parallel = FALSE,
  Interactive = FALSE
)

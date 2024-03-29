#### Setting OM parameters

Nsp <- 6

### Scenario 1:
Sim1 <- list(

  ### Number of simulation years (time step is a month)
  n_years = 15,

  ### Habitat (depth) parameters. Very important as it controls how species are distributed in space, thus affects both
  ### the population dynamics and the vessel dynamics
  Range_X = 1:20,           # the x-axis extent. the bigger, the more complex topography you can get at the cost of simulation time
  Range_Y = 1:40,           # the y-axis extent. the bigger, the more complex topography you can get at the cost of simulation time
  SD_O = 100,						    # SD of depth observation in space. the bigger, the more variable is the depth
  SpatialScale = 3,					# the spatial correlation range. The bigger, the more correlated are depth in space (patchy if small).
                            # N.B: Absolute value of high/low depends on the extent of x-axis and y-axis
  year.depth.restriction = 1000,	# If MPA, then years when MPA implemented (default > n_years)
  depth.restriction = c(50,150),	# location of MPA

  ### Species parameters
  n_species  = Nsp, # 4 targets and 2 for M and F Atlantic wolffish (different seasonality).
  sp_groups = c(1, 2, 3, 4, 5, 5), ## default seq(length.out = Nsp)
  ## ##############################
  ## col1: haddock
  ## col2: saithe
  ## col3: cod
  ## col4: golden redfish
  ## col5: Atlantic wolfish Females
  ## col6: Atlantic wolfish Males
  ## ##############################

  # price_fish = matrix(rep(c(1.5,1,-1,2,0),15), ncol=5, byrow=T),   # price per species across time (roughly in Us west coast)
  price_fish = matrix(rep(c(13.33, # 2021 average (round, fresh), https://www.rafisklaget.no/statistikk-detaljer
                            8.56,
                            23.89,
                            12.61,
                            rep(12.16, 2)),
                          15),
                      ncol = Nsp, byrow = TRUE),   # random

  ## Fish habitat preference/movement control parameters
  func_mvt_dist = c(rep("Exp", Nsp)),					# shape of mvt p_function
  func_mvt_depth = c(rep("Lognorm", Nsp)),    # shape of depth preference p_function (choice of Normal, Exp, Lognormal, Uniform)
  func_mvt_lat = c(rep("Unif", Nsp)),					# shape of range preference p_function (choice of Normal, Exp, Lognormal, Uniform)

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
                            nrow=12, ncol = Nsp, byrow = TRUE),
  Fish_dist_par2 = rep(0, Nsp), 										# mvt distance (not used)
  Fish_depth_par1 = matrix(c(87, 100, 95, 200, 140, 140, # depth preference mean per month 1
                             87, 100, 95, 200, 140, 140, # month2
                             87, 100, 95, 200, 140, 140, # month3
                             87, 100, 95, 200, 140, 140, # month4
                             87, 100, 95, 200, 140, 140, # month5
                             87, 100, 95, 200, 140, 140, # month6
                             87, 100, 95, 200, 87, 87, # month7
                             87, 100, 95, 200, 87, 87, # month8
                             87, 100, 95, 200, 87, 87, # month9
                             87, 100, 95, 200, 87, 87, # month10
                             87, 100, 95, 200, 140, 87, # month11
                             87, 100, 95, 200, 140, 87),# month12
                           nrow=12, ncol = Nsp, byrow = TRUE),
  Fish_depth_par2 = matrix(c(0.63, 0.57, 0.59, 0.35, 0.9, 0.9, # depth preference sd log scale per month 1
                             0.63, 0.57, 0.59, 0.35, 0.9, 0.9, # month2
                             0.63, 0.57, 0.59, 0.35, 0.9, 0.9, # month3
                             0.63, 0.57, 0.59, 0.35, 0.9, 0.9, # month4
                             0.63, 0.57, 0.59, 0.35, 0.9, 0.9, # month5
                             0.63, 0.57, 0.59, 0.35, 0.9, 0.9, # month6
                             0.63, 0.57, 0.59, 0.35, 0.48, 0.48, # month7
                             0.63, 0.57, 0.59, 0.35, 0.48, 0.48, # month8
                             0.63, 0.57, 0.59, 0.35, 0.48, 0.48, # month9
                             0.63, 0.57, 0.59, 0.35, 0.48, 0.48, # month10
                             0.63, 0.57, 0.59, 0.35, 0.9, 0.48, # month11
                             0.63, 0.57, 0.59, 0.35, 0.9, 0.48),# month12
                           nrow=12, ncol = Nsp, byrow = TRUE),
  Fish_range_par1 = rep(0, Nsp),                    # X-axis range min
  Fish_range_par2 = rep(50,Nsp),                    # X axis range max
  Rangeshift_catchability_adjust = matrix(c(1, 1, 1, 1, 1, 1, # Adjust the catchability so that it mimics some changes in availability of the population (migration)
                                            1, 1, 1, 1, 1, 1, # month2
                                            1, 1, 1, 1, 1, 1, # month3
                                            1, 1, 1, 1, 1, 1, # month4
                                            1, 0.8, 1, 1, 1, 1, # month5
                                            1, 0.8, 1, 1, 1, 1, # month6
                                            1, 0.8, 1, 1, 1, 1, # month7
                                            1, 0.8, 1, 1, 1, 0.4, # month8
                                            1, 0.8, 1, 1, 1, 0.4, # month9
                                            1, 0.8, 1, 1, 1, 0.4, # month10
                                            1, 0.8, 1, 1, 1, 0.4, # month11
                                            1, 0.8, 1, 1, 1, 0.4),# month12
                                            nrow=12, ncol= Nsp, byrow = TRUE),

  ## Species Schaefer pop dyn parameters
  B0 = c(80000000, 220000000, 148000000 + 460000000,	# NEA Haddock, saithe, cod (NEA+coastal),
         68600000,                                    # golden redfish,
         c(10000,10000) * 1000 / 2),                  # and atlantic wolffish (1/2)x2 (~ about 1
                                        # order of magnitude lower).
  Pop_ratio_start = rep(0.5, Nsp),          # this is the "rough" abundance of the population at the start of the simulation
  r = c(0.175, 0.175, 0.2,
        0.02625, rep(0.14, 2)) / 12, # NEA Haddock, saithe, cod (NEA+coastal), golden redfish, and atlantic wolffish (1/2)x2
  sigma_p = c(0.688, 0.45, 0.31,     # NS saithe as based on R table, unlike NEA.
              1.198, rep(0.5, 2)),   # log-normal process error = recruitment to fishery
                                        # Arbitrary for wolffish
  sigma_p_timing= c(4, 3, 2, 4, rep(8, 2)), # when does recruitment happen?
  fish.mvt = TRUE,                          # whether animals redistribute annually based on habitat preference

  ### Parameters controlling the vessel dynamics
  Nregion = c(2,2),         # nb of fishing regions for the fishermen (equal cut of grid along X and Y)
  Nvessels = 50,					  # nb vessels
  Tot_effort = 1000,				# end year nb of effort
  Tot_effort_year = seq(1000, 1000, length.out = 15),
  CV_effort = 0.2, 					# CV of effort around the yearly mean effort
  qq_original = c(0.00348, 0.0251, 0.0243,
                  0.0172, rep(0.00157 / 2, 2)) * 1e-3, # the average catchability coef by species for ALL vessels
                                        # (does not mean anything. just a scaler)
                                        # Estimated to approximately match 2021 catches with an arbitrary 1000 effort
  catch_trunc = c(rep(0, 6)),           # truncate value if 1 otherwise keep continuous
  CV_vessel = 0.1,					# the between vessel variability in mean catchability
  vessel_seeds = 10,        # this creates heterogeneity in the sampled vessels
  Effort_alloc_month = c(5,4,3,2,1,1,1,2,3,4,5,5),  # it is rescaled to sum to 1
  do.tweedie = TRUE,				# include observation error in vessel catchability
  # xi = c(1.7,1.7,1.2,1.5),  # power of the tweedie distribution. reduce this to have more 0s. ==0 or >0
  xi = c(1.9,1.9,1.9,1.9, 1.9, 1.9),  # power of the tweedie distribution. reduce this to have more 0s. ==0 or >0
  # phi = c(0.2,0.2,0.2,0.2), # the scaler of the tweedie distribution
  phi = c(0.5,0.5,0.5,0.5, 0.5, 0.5), # the scaler of the tweedie distribution
  Preference = 1, 					# controls how effort concentrate into areas. The higher, the more concentrated is the effort
  Changing_preference = FALSE, 		# whether effort concentration changes over time
  Discard_rate_beta1 = rep(0, Nsp),    # the mean discard rate (of a beta dsitribution). for a 100% discard rate write 1, for 0% discard rate write 0
  Discard_rate_beta2 = rep(0,Nsp),    # the variance of the discard rate (of a beta dsitribution). for a 100% discard rate write 0, for 0% discard rate write 0
  Discard_survival_rate1= rep(0, Nsp),   # the mean survival rate (of a beta dsitribution). for a 100% survival rate write 1, for 0% discard rate write 0
  Discard_survival_rate2= rep(0, Nsp),   # the variance  survival rate (of a beta dsitribution). for a 100% survival rate write 0, for 0% discard rate write 0

  ### Parameter controlling the resampling procedure from the whole fleet
  samp_prob = 0.05,            # this is the sampling probability
  samp_unit = "vessel",        # the samping unit: "vessel" or "fishing" (i.e. random sample from all fishing events)
  samp_seed = 123,             # the seed for reproducibility of the samples taken
  samp_mincutoff = 0,          # cut-off value to round values to 0
  start_year = 5,              # when to start sampling from the fleet
  months = c(1:12),            # the month from which you want to take the sample from

  ### Other OM control features
  plotting = FALSE,
  parallel = FALSE,
  Interactive = FALSE
)

##' This is the main configuration file that controls both the population and the vessel dynamics
##'
##' @param n_years is the number of simulation years
##' @param n_species is the catch data frame
##' @param conf  configurations
##' @param rerun_stack  force the rerun of the code to create the prediction grid
##' @details Prepare and modifies the original data files to be able to run the spatio-temporal model
##' @return Data to be provided to TMB
##' @export
##'

Sim_Settings1 <- list(

  ### Number of simulation years
  n_years = 30,

  ### Species parameters
  n_species  = 4,
  price_fish = matrix(rep(c(1,1,-1,1),30), ncol=4, byrow=T),   # price per species across time

  ## Fish habitat preference/movement control parameters
  func_mvt_dist = c("Exp", "Exp", "Exp", "Exp"),					# shape of mvt p_function
  func_mvt_depth = c("Lognorm", "Lognorm", "Lognorm", "Lognorm"),	# shape of depth preference p_function (choice of Normal, Exp, Lognormal, Uniform)
  func_mvt_lat = c("Unif", "Unif", "Unif", "Unif"),					# shape of range preference p_function (choice of Normal, Exp, Lognormal, Uniform)

  Fish_dist_par1 = c(rep(0.4, 3), 1),							# mvt distance par1
  Fish_dist_par2 = rep(0, 4), 										# mvt distance par2
  Fish_depth_par1 = c(220, 450, 250, 1200),       # depth preference par1
  Fish_depth_par2 = c(1, 1, 0.3, 0.5),            # depth preference par2
  Fish_range_par1 = rep(0, 4),                    # X-axis range par1
  Fish_range_par2 = rep(50,4),                    # X axis range par2

  ## Species Schaefer pop dyn parameters
  B0 = c(50000,820000,36000,540000),							# based on petrale, dover, darkblotched, sablefish (mt)
  r = c(0.22, 0.17, 0.08, 0.08),                  # based on petrale, dover, darkblotched, sablefish (mt)
  no.mvt = FALSE,				# whether animals redistribute annually based on habitat preference

  ### Parameters controlling the vessel dynamics
  Nvessels = 30,					  # nb vessels
  Tot_effort = 200,					# end year nb of effort
  CV_effort = 0.2, 					# CV of effort around the yearly mean effort
  qq_original = 0.2,				# the average catchability coef for ALL vessels
  CV_vessel = 0.05,					# the between vessel variability in mean catchability
  do.tweedie = TRUE,				# include observation error in catchability
  xi = c(1.9,1.9,1.2,1.2),  # power of the tweedie distribution
  phi = c(0.2,0.2,0.2,0.2), # the scaler of the tweedie distributin
  Preference = 1, 					# controls how effort concentrate into areas
  Changing_preference = FALSE, 		# whether effort concentration changes over time

  ### Habitat (depth) parameters
  Range_X = 1:30,
  Range_Y = 1:30,
  SD_O = 100,						    # SD of spatial observation
  SpatialScale = 3,					# scale of spatial variation
  year.depth.restriction = 1000,	# If MPA, then years when MPA implemented (default > n_years)
  depth.restriction = c(50,150),	# location of MPA

  ### Other OM control features
  plotting = FALSE,
  Interactive = FALSE
)										# End of simulation setting


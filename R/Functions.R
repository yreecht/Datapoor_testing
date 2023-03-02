##' Function to generate full mvt matrix of an animal: how likely they move from one cell to another
##' this is a wrapper of the "Mvt_upd1" function to calculate the movement probability between all cells
##'
##' @param Pop_param is the list of population parameters taken from the user specific simulation setting
##' @param Depth_eff logical. do we want to include the habitat preference effect. default is TRUE
##' @param Dist_eff  logical. do we want to include the distance effect (i.e. movement is restricted by distance). default is TRUE
##' @param Lat_eff  logical. do we want to include additional constraint on movement based on latitude
##' (this was initially set as potential way to set up and no-fish zone)
##' @param Depth_func  is the shape of the depth preference (if depth effect is activated): choice between "Norm", "Exp", "Lognorm", "Unif"
##' @param Dist_func  is the shape of the distance function (only "Exp" available at the moment)
##' @param Lat_func  is the shape of the depth preference (if latitude effect is activated): choice between "Norm", "Exp", "Lognorm", "Unif"
##' @param Eastings  is the extent of the grid on the x-axis: all the Eastings values
##' @param Northings  is the extent of the grid on the y-axis: all the Northings values
##' @param data.bathym  is the bathymetric map of the study area

##'
##' @return A movement probability matrix between every cell in the map
##' @details
##'
##' @export
##'
##'
##
Mvt_matrix <- function(Pop_param, Depth_eff, Dist_eff, Lat_eff, Depth_func, Dist_func="Exp", Lat_func, Eastings, Northings, data.bathym){
  Mvt_mat <- matrix(0,(length(Eastings)*length(Northings)),(length(Eastings)*length(Northings)))
  for (loc in 1:(length(Eastings)*length(Northings))){
    Mvt_mat[loc,] <- Mvt_upd1(loc, Pop_param, Depth_eff=Depth_eff, Dist_eff=Dist_eff, Lat_eff=Lat_eff, Depth_func=Depth_func, Dist_func=Dist_func, Lat_func=Lat_func, data.bathym)
  }
  return(Mvt_mat)
}


##' Function that calculates the probability of a fish to move from a specified cell "loc" to the rest of the simulation area
##' The movement probability is calculated based on the product of habitat suitability (depth) and
##' movement capability (distance between cells and more restriction if any)
##' How a fish prefers some habitat (i.e. the shape of the habitat preference function) and how far they can move are input parameters to this function
##'
##'
##' @param loc is the cell of origin from which the individual move away
##' @param Pop_param is the list of population parameters taken from the user specific simulation setting
##' @param Depth_eff logical. do we want to include the habitat preference effect. default is TRUE
##' @param Dist_eff  logical. do we want to include the distance effect (i.e. movement is restricted by distance). default is TRUE
##' @param Lat_eff  logical. do we want to include additional constraint on movement based on latitude
##' (this was initially set as potential way to set up and no-fish zone)
##' @param Depth_func  is the shape of the depth preference (if depth effect is activated): choice between "Norm", "Exp", "Lognorm", "Unif"
##' @param Dist_func  is the shape of the distance function (only "Exp" available at the moment)
##' @param Lat_func  is the shape of the depth preference (if latitude effect is activated): choice between "Norm", "Exp", "Lognorm", "Unif"
##' @param data.bathym  is the bathymetric map of the study area
##'
##' @return the vector of movement probability from the "loc" cell to the other cells
##' the order of the vector is the same as the bathymetry map (i.e. x-axis first)
##'
##' @details
##' @export
##'
##'
##
##
Mvt_upd1 <- function(loc, Pop_param, Depth_eff, Dist_eff, Lat_eff, Depth_func, Dist_func="Exp", Lat_func, data.bathym){
  Di1 <- Pop_param[1]  	# 1st value for dist function
  Di2 <- Pop_param[2]  	# 2nd value for dist function
  De1 <- Pop_param[3]  	# 1st value for depth function
  De2 <- Pop_param[4]	 	# 2nd value for depth function

  Dist <- sqrt((data.bathym$X-data.bathym$X[loc])^2+(data.bathym$Y-data.bathym$Y[loc])^2)

  Prob1 <- 1
  if(Dist_eff=="TRUE"){
    Prob1 <- dexp(Dist, rate=1/Di1)
  }

  Prob2 <- 1
  if(Depth_eff=="TRUE"){
    if(Depth_func=="Norm") Prob2 <- dnorm(data.bathym$depth, De1, De2*De1)
    if(Depth_func=="Exp") Prob2 <- dexp(data.bathym$depth, rate=1/De1)
    if(Depth_func=="Lognorm") { SD <- sqrt(log(De2^2+1)); Prob2 <- dlnorm(data.bathym$depth, (log(De1)-SD^2/2), SD) }
    if(Depth_func=="Unif") Prob2 <- dunif(data.bathym$depth, De1, De2)
  }

  Prob3 <- 1
  # if(Lat_eff=="TRUE"){
  #   if(Lat_func=="Norm") Prob3 <- dnorm(data.bathym$X, Lat1, Lat2)
  #   if(Lat_func=="Exp") Prob3 <- dexp(data.bathym$X, rate=1/Lat1)
  #   if(Lat_func=="Lognorm") Prob3 <- dlnorm(data.bathym$X, (log(Lat1)-Lat2^2/2), Lat2)
  #   if(Lat_func=="Unif") Prob3 <- dunif(data.bathym$X, Lat1, Lat2)
  # }

  Prob <- Prob1*Prob2*Prob3/sum(Prob1*Prob2*Prob3)
  return(Prob)
}



##' This function calculates the initial population distribution = "near-equilibrium" population
##' it begins with an equal distribution of fish in space then based on the movement probabilities, it redistributes fish.
##' It is an iterative process and the function recalculates the distribution over 1000iterations. In practice, this should be enough to reach
##' some sort of an equilibirum population distribution
##'
##' Note that the notion of "near-equilibrium" is not strict, especially if there is seasonality in movement. And this is ignored here.
##'
##' It is to be noted that this method works for most cases EXCEPT when the species have very short mobility. If which case, the species distribution
##' will not reach equilibirum and the distribution you get out from this function is entirely dependent on the intial distribution you provided to
##' function
##'
##'
##' @param Mvt_mat_adult is the movement probability matrices between grid cells (only provide for 1 specific season/month)
##' @param B0 This is the "fictive", random initial population you provide the model
##' @return the "near-equilibrium" population distribution. As noted in the description, it is not the equilibirum distribution in strict term.
##' The intention is to have a distribution that is sensible based on the underlying habitat map that was provided
##'
##' @details
##' @export
##'
##'
##
##
##
Stable_pop_dist <- function(Mvt_mat_adult, B0, Range_X, Range_Y){
  Pop1 <- rep(B0/(length(Range_X)*length(Range_Y)), (length(Range_X)*length(Range_Y)))
  Mvt_mat_adult1 <- t(Mvt_mat_adult)
  for (Time in 1:1000){
    Pop <- Mvt_mat_adult1%*%Pop1
    Pop1 <- Pop
  }
  Initpop_adult_mat <- matrix(Pop1, nrow=length(Range_X), ncol=length(Range_Y))
  return(Initpop_adult_mat)
}




##' Function to sample from the whole fishing fleet (from the true whole data):
##' @param data is the whole fleet data that was generated in the simulation
##' @param percent is the percent of the data to sample from
##' @param unit is the sampling unit of the data. choice of "vessel" or "unit" (fishing event)
##' @param months  is used to specify the months at which sampling occurs
##' @param seed  is the seed used for the sampling process: for reproducibility
##' @details Applies the specified sampling protocol to the original generated data from the whole fleet
##' @return the final sample data to be fed into the CPUE standardization exercise
##' @export
##'
##'

sampling_select <- function(data, percent=0.1, unit="vessel", seed = 123, months = c(1,3,6,8,10)){
  set.seed(seed)
  if (unit == "vessel") {
    vessel_sel = sample(unique(data$vessel), floor(length(unique(data$vessel))*percent))
    dat <- data %>% filter(month %in% months, vessel %in% vessel_sel)
  }
  if (unit == "fishing") {
    dat <- dat %>% filter(month %in% months)
    obs_sel = sample(1:nrow(dat), floor(nrow(dat)*percent))
    dat <- dat[obs_sel,]
  }
  return(dat)
}




##' Function to make some plot from the sampled data
##' @param data is the data you want to make the plots
##' @param years is the percent of the data to sample from
##' @details Produces the plot
##' @return a plot
##' @export
##'
##'
plot_fishing <- function(data, years, ...){
  data %>% filter(year %in% years) %>% group_by(X,Y,year) %>% summarize(total_effort = n()) %>%
    ggplot(., aes(x=X, y=Y)) +
    facet_wrap(.~year)  + theme_bw() +
    geom_tile(data = Data$bathym, aes(x=X, y=Y, col=depth, fill=depth)) +
    scale_fill_viridis_c() + scale_color_viridis_c()+ geom_point()
}


##' Function to create arbitrary fishing area to be later used in CPUE standardization
##' @param ncut_x is the number of division along the x-axis for the fishing area: equal division only
##' @param ncut_y is the number of division along the y-axis for the fishing area: equal division only
##' @param data is the data of interest
##' @details Applies the specified area division and add it to the original data
##' @return the data wih the fishing area specification
##' @export
##'
##'
equal_partition <- function(data, ncut_x=2, ncut_y=2, Range_X, Range_Y){
  out <- expand.grid(X = Range_X, Y = Range_Y)
  if (ncut_x>1) out <- out %>% mutate(X_cut = factor(cut(X, ncut_x), labels=c(1:ncut_x))) else out <- out %>% mutate(X_cut = rep(1, length(X)))
  if (ncut_y>1) out <- out %>% mutate(Y_cut = factor(cut(Y, ncut_y), labels=c(1:ncut_y))) else out <- out %>% mutate(Y_cut = rep(1, length(Y)))
  area <- as.factor(apply(out[,3:4], 1, function(x) paste(x, collapse="_")))
  area <- factor(area, labels = 1:(ncut_x*ncut_y))
  area <- data.frame(X=out$X, Y=out$Y, area)

  data <- data %>% left_join(area)
  return(data)
}


##' Convenience function to move from gridID (cell ID) to X, Y coordinates and vice versa
##'
##' @param gridID is the gridID (write this is you want to )
##' @param X is the X coordinate of the gridID
##' @param Y is the Y coordinate of the gridID
##'
##' @details
##' @example
##'
Convert_grid_coordinates <- function(GridID=NULL, X=NULL, Y=NULL, Settings){
  if (!is.null(GridID)) {
    X = ((GridID-1) %% length(Settings$Range_X)) + 1
    Y = trunc((GridID-1) / length(Settings$Range_X)) + 1
    return(c(X,Y))
  } else {
    GridID = (Y-1)*length(Settings$Range_X) + X
    return(GridID)
  }
}

##' Convenience function to aggregate parameters/status (carrying capacity, biomass,...)
##' for species simulated through several groups.
##'
##' Internally used to calculate common growth rates for species simulated as several groups
##' (typically for separate simulations of males and females).
##' @title Function \code{group_cols()}
##' @param mat A matrix with groups in column.
##' @param groups Vector of groups (same value for grouped columns), with as many elements as columns in \code{mat}.
##' @param FUN The aggregation function (sum by default).
##' @return A matrix of the same dimensions as \code{mat}, where columns of common groups are summed (by default).
##' @author Yves Reecht
group_cols <- function(mat, groups, FUN = sum)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  2 Mar 2023, 12:38
    require(dplyr)

    res <- sapply(seq(length.out = ncol(mat)),
                  function(i, mat, groups, fun)
           {
               apply(mat[ , groups %in% groups[i], drop = FALSE],
                     1, FUN = fun)
           }, mat = mat, groups = groups, fun = FUN, simplify = FALSE) %>%
        bind_cols(.name_repair = "unique") %>% as.matrix()

    colnames(res) <- colnames(mat)

    return(res)
}

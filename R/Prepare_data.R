##' Prepare/filter data for the TMB model run
##' @param data is the catch data frame
##' @param Pred_data is the prediction data frame
##' @param conf  configurations
##' @details Prepare and modifies the original data files to be able to run the spatio-temporal model
##' @return Data to be provided to TMB
##' @export
##'


Prepare_data = function(data, conf, Pred_data) {

  ### Step 1: extracting the variables included in the formula and divide between fixed and random effects
  data$yy <- as.numeric(data[,'Sp1'])

  fixed_effect_formula = nobars(conf$model_formula)
  RE_effects = barnames(findbars(conf$model_formula))
  mgcv_mod <- mgcv::gam(fixed_effect_formula, data = data)
  X <- model.matrix(mgcv_mod)


  RE_indexes = matrix(rep(0L,4),2,2)
  nobs_RE = 0L
  ln_tau_G_index= 0L
  if (length(RE_effects)>0)
  {
    RE_indexes <- sapply(RE_effects, function(x) as.numeric(factor(data[[x]], labels=1:length(unique(data[[x]]))))) - 1L
    RE_indexes <- as.matrix(RE_indexes)
    nobs_RE <- unname(apply(RE_indexes, 2L, max)) + 1L
    if (length(nobs_RE) == 0L) nobs_RE <- 0L
    ln_tau_G_index <- unlist(lapply(seq_along(nobs_RE), function(i) rep(i, each = nobs_RE[i]))) - 1L
  }

  ### Extract general feature from the data
  yobs= data %>% select(starts_with("sp"))
  Nspecies = ncol(yobs)
  Nyear <- length(unique(data$year))

  ### Creating the projection data i.e. combination of different variables. Need to choose a reference category
  X_pred <- mgcv::predict.gam(mgcv_mod, type = "lpmatrix", newdata = Pred_data)
  RE_indexes_proj = matrix(rep(0L,4),2,2)
  if (length(RE_effects)>0)
  {
    temp <- data.frame(val=unique(data[, RE_effects]), new_val = 1:length(unique(data[, RE_effects])) - 1L)
    RE_indexes_proj <- as.matrix(temp$new_val[match(Pred_data[, RE_effects], temp$val)])
    RE_indexes_proj[is.na(RE_indexes_proj)] <- 999L
    # RE_indexes_proj = as.matrix(RE_indexes_proj)
  }

  ### Stuff needed for spatiotemporal A matrix.
  if (conf$spatial_model == 0){
    # adding fake coordinates (as a placeholder)
    data$X <- 1:nrow(data)
    data$Y <- 1:nrow(data)
    spde <- make_mesh(data, xy_cols = c("X", "Y"), cutoff = conf$INLAmesh_cutoff)
  } else {
    # Construct our mesh:
    spde <- make_mesh(data, xy_cols = c("X", "Y"), cutoff = conf$INLAmesh_cutoff)
  }
  plot(spde)

  if (conf$barrier == 1){
    # create the barrier mesh
    Barrier_range_fraction = 0.1

    spde <- add_barrier_mesh(
      spde, Atlantic_proj, range_fraction = Barrier_range_fraction,
      proj_scaling = 1, plot = TRUE)

  }

  data$sdm_orig_id <- seq(1, nrow(data))
  data$sdm_x <- spde$loc_xy[,1,drop=TRUE]
  data$sdm_y <- spde$loc_xy[,2,drop=TRUE]
  fake_data <- unique(data.frame(sdm_x = data$sdm_x, sdm_y = data$sdm_y))
  fake_data[["sdm_spatial_id"]] <- seq(1, nrow(fake_data))
  data <- base::merge(data, fake_data, by = c("sdm_x", "sdm_y"),
                      all.x = TRUE, all.y = FALSE)
  data <- data[order(data$sdm_orig_id),, drop=FALSE]

  A <- spde$A
  A_st <- INLA::inla.spde.make.A(spde$mesh,
                                 loc = as.matrix(fake_data[, c("sdm_x", "sdm_y"), drop = FALSE]))

  n_s <- nrow(spde$mesh$loc)
  Nmesh <- n_s


  ### Now prepare the input data for the prediction

  # design matrix for prediction for the positive model
  if (length(grep( "k =", fixed_effect_formula)) == 0 )  X_proj <- model.matrix(fixed_effect_formula, Pred_data)
  if (length(grep( "k =", fixed_effect_formula)) > 0 )   X_proj = mgcv::predict.gam(mgcv_mod, type = "lpmatrix", newdata = Pred_data)

  proj_data <- subset(Pred_data, year==Sim2$start_year)

  # create the A matrix: it is always the same because the underlying mesh does not change
  A_proj<- INLA::inla.spde.make.A(spde$mesh,
                                  loc = as.matrix(proj_data[, c("X", "Y"), drop = FALSE]))

  Npred = nrow(A_proj)

  Pred_data$A_spatial_index_proj <- as.factor(apply(Pred_data[,c('X','Y')],1,function(x) paste(x, collapse="_")))
  Pred_data$A_spatial_index_proj <- as.integer(factor(Pred_data$A_spatial_index_proj, labels=(1:length(unique(Pred_data$A_spatial_index_proj)))))

  ### Starting model building

  ## tmb data
  tmb_data <- list(
    Nfactor               = conf$Nfactor,
    X                     = as.matrix(X),
    yobs                  = as.matrix(yobs),
    RE_indexes            = RE_indexes,
    nobs_RE               = nobs_RE,
    ln_tau_G_index        = ln_tau_G_index,
    family                = 0L,
    link                  = 1L,
    sim                   = conf$sim,
    incl_target           = conf$incl_target,
    incl_sp_int           = conf$incl_sp_int,
    spatial_model         = conf$spatial_model,
    include_st            = conf$include_st,
    barrier               = conf$barrier,
    do_spatialDFA         = conf$do_spatialDFA,
    spde                  = spde$spde$param.inla[c("M0","M1","M2")],
    # spde_barrier        = make_barrier_spde(spde),
    # barrier_scaling     = if(conf$barrier==1) spde$barrier_scaling else 1,
    Aobs                  = A,
    Ast                   = A_st,
    A_spatial_index       = data$sdm_spatial_id - 1L,
    year_i                = make_year_i(data$year),
    Nyear                 = Nyear,
    Nmesh                 = Nmesh,
    Npred                 = Npred,
    Do_predict            = conf$Do_predict,
    se_fit                = conf$se_fit,
    X_proj                = as.matrix(X_pred),
    A_proj					      = A_proj,
    A_spatial_index_proj  = Pred_data$A_spatial_index_proj-1L,
    RE_indexes_proj       = RE_indexes_proj
  )


  # tmb param specs
  vals = rnorm(Nspecies*(tmb_data$Nfactor))
  vals = rnorm(Nspecies*(tmb_data$Nfactor) - (tmb_data$Nfactor)*(tmb_data$Nfactor-1)/2)

  Nsize = ifelse(conf$do_spatialDFA == 0, Nspecies, conf$Nfactor)

  tmb_params <- list(
    beta                = matrix(0, nrow=ncol(X), ncol=Nspecies),
    omega               = matrix(0, nrow=Nmesh, ncol=Nsize),
    epsilon_st          = array(0, dim=c(Nmesh, Nyear, Nsize)),
    transf_rho          = rep(0,Nsize),
    logKappa            = rep(0.1,Nsize),
    logTauO             = rep(0.1, Nsize),
    logTauE             = rep(0.1,Nsize),
    Species_target      = rep(1,Nspecies),
    L_val_spatial       = vals,
    L_val_target        = vals,
    L_val_sp            = rnorm(Nspecies*(Nspecies-1)/2) ,
    RE_species_latent   = matrix(0, nrow=nrow(yobs), ncol=tmb_data$Nfactor) ,
    RE_species_int      = matrix(0, nrow=nrow(yobs), ncol=Nspecies) ,
    logsds              = rep(0, Nspecies),
    thetaf              = rep(0.1, Nspecies) ,
    ln_phi              = rep(0.1, Nspecies) ,
    ln_tau_G            = if(sum(nobs_RE)>0) matrix(0.1, nrow = length(nobs_RE), ncol=Nspecies) else matrix(0.1,nrow=1,ncol=Nspecies) ,
    RE                  = if(sum(nobs_RE)>0) matrix(0, nrow=sum(nobs_RE), ncol=Nspecies) else matrix(0,nrow=1,ncol=Nspecies)
  )

  ### Output list from the Prepare_data function
  output <- list()
  output$tmb_data <- tmb_data
  output$tmb_params <- tmb_params
  output$Pred_data <- Pred_data
  output$mgcv_mod <- mgcv_mod
  output$data <- data
  output$spde <- spde

  return(output)

}


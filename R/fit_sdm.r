##' Run the model with the specified configurations
##' @param data$data$tmb_data data input for the TMB model run
##' @param data$data$tmb_paramss parameter input for the TMB model run
##' @param conf  configurations file for the model and data
##' @details the model is run in two steps to make model parameters estimable
##' @return different outputs needed for summarizing and plotting the results at later stage
##' @useDynLib mackerel_ST
##' @export fit_mackerelST
##'


fit_sdm <- function(data, conf) {

	### model to use for the computation
  version <- paste0(getwd(), "/src/mult_species"); dllversion = "mult_species"

  dyn.load(dynlib(version))
  dyn.unload(dynlib(version))
  if (.Platform$OS.type == "unix") compile(paste0(version, ".cpp"),"-O0 -g")
  if (.Platform$OS.type == "windows") compile(paste0(version, ".cpp"), "-O1 -g", DLLFLAGS="")
  dyn.load(dynlib(version))

  ### Now changing configuration depending on the specification put in the configuration files
  Map_phase1 <- list()

		# configuration depending on the how the spatial info is treated: non spatial model, spatial model, spatio-temporal model
			if (conf$spatial_model == 0) {
				Map_phase1$omega <- factor(rep(NA, length(data$tmb_params$omega)))
				Map_phase1$epsilon_st <- factor(rep(NA, length(data$tmb_params$epsilon_st)))
				Map_phase1$transf_rho <- factor(rep(NA, length(data$tmb_params$transf_rho)))
				Map_phase1$logKappa <- factor(rep(NA, length(data$tmb_params$logKappa)))
				Map_phase1$logTauO <- factor(rep(NA, length(data$tmb_params$logTauO)))
				Map_phase1$logTauE <- factor(rep(NA, length(data$tmb_params$logTauE)))
			}
			if (conf$spatial_model == 1 & conf$include_st == 0) {
				Map_phase1$epsilon_st <- factor(rep(NA, length(data$tmb_params$epsilon_st)))
				Map_phase1$transf_rho <- factor(rep(NA, length(data$tmb_params$transf_rho)))
				Map_phase1$logTauE <- factor(rep(NA, length(data$tmb_params$logTauE)))
			}
			if (conf$spatial_model == 1 & conf$include_st == 1) {
				Map_phase1$transf_rho <- factor(rep(NA, length(data$tmb_params$transf_rho)))
			}
		# if not including the fisher targeting behavior
		if (conf$incl_target == 0) {
			Map_phase1$L_val_target <- factor(rep(NA, length(data$tmb_params$L_val_target)))
			Map_phase1$RE_species_latent <- factor(rep(NA, length(data$tmb_params$RE_species_latent)))
		}
		# if not including the spatial DFA
		if (conf$do_spatialDFA == 0) {
			Map_phase1$L_val_spatial <- factor(rep(NA, length(data$tmb_params$L_val_spatial)))
		}
		# if including the spatial DFA (for identifiability. the scale of the spatial variation is controlled by the L_val_spatial values)
		if (conf$do_spatialDFA == 1) {
		  Map_phase1$logTauO <- factor(rep(NA, length(data$tmb_params$logTauO))) # did not work well. tauE is important
		  Map_phase1$logTauE <- factor(rep(NA, length(data$tmb_params$logTauE))) # did not work well. tauE is important
		}
		# if not including species interaction
		if (conf$incl_sp_int == 0) {
			Map_phase1$logsds <- factor(rep(NA, length(data$tmb_params$logsds)))
			Map_phase1$L_val_sp <- factor(rep(NA, length(data$tmb_params$L_val_sp)))
			Map_phase1$RE_species_int <- factor(rep(NA, length(data$tmb_params$RE_species_int)))
		}
		if (conf$incl_sp_int == 1) {
			# Map_phase1$logTauO <- factor(rep(NA, length(data$tmb_params$logTauO))) # did not work well. tauE is important
			# Map_phase1$logTauE <- factor(rep(NA, length(data$tmb_params$logTauE))) # did not work well. tauE is important
			# Map_phase1$logsds <- factor(rep(NA, length(data$tmb_params$logsds)))
		}
		# configuring the random effect estimation in the model
		if(sum(data$tmb_data$nobs_RE)==0) {
		  Map_phase1$RE <- factor(rep(NA, length(data$tmb_params$RE)))
		  Map_phase1$ln_tau_G <- factor(rep(NA, length(data$tmb_params$ln_tau_G)))
		}
		random_list <- NULL
		if(sum(data$tmb_data$nobs_RE)>0)  random_list <- c(random_list, "RE")
		if(data$tmb_data$incl_target == 1)  random_list <- c(random_list, "RE_species_latent")
		if(data$tmb_data$incl_sp_int == 1)  random_list <- c(random_list, "RE_species_int")
		if(conf$spatial_model == 1 & conf$include_st >= 0)  random_list <- c(random_list, "omega")
		if(conf$spatial_model == 1 & conf$include_st >0)  random_list <- c(random_list, "epsilon_st")

  ### compile the model
  tmb_obj_phase1 <- MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase1,
                              random = random_list, DLL = dllversion, silent = TRUE)

	### Run the model
  opt_phase1 <- fit_tmb(tmb_obj_phase1, lower=-15, upper=15, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))


  opt_phase1$diagnostics


  # sdrep <- NA #sdreport(tmb_obj_phase2)
  # #
  # dens_pred <- NA
  # epsilon_pred <- NA
  #
  # if (conf$Do_predict == 1) {
  #   if (conf$bias_correct == TRUE) sdrep <- sdreport(tmb_obj_phase1, bias.correct =TRUE, bias.correct.control = list(sd= TRUE))
  #   if (conf$bias_correct == FALSE) sdrep <- sdreport(tmb_obj_phase1, bias.correct =FALSE)
  #   # extract the density estimates on the prediction grid
  #   pred_grid = data$pred_grid
  #   mu_proj2_new <- tmb_obj_phase1$report()$mu_proj
  #   dens_pred <- data.table::as.data.table(mu_proj2_new)
  #   colnames(dens_pred) <- c("LOC", "YEAR", "AGE", "Pred")
  #   dens_pred$YEAR <- as.factor(dens_pred$YEAR)
  #   dens_pred$YEAR <- as.factor(as.numeric(as.character(factor(dens_pred$YEAR, labels=sort(unique(data$data$YEAR))))))
  #   pred_grid$YEAR <- as.factor(as.numeric(as.character(factor(pred_grid$YEAR, labels=sort(unique(data$data$YEAR))))))
  #   pred_grid$LOC <- rep(1:(nrow(pred_grid)/conf$Nyear), conf$Nyear)
  #   dens_pred <- dens_pred %>% left_join(pred_grid)
  #   dens_pred$logPred <- log(dens_pred$Pred)
  #
  #   # extract estimated spatio-temporal field on the prediction grid
  #   epsilon2_new <- tmb_obj_phase1$report()$epsilon_st_A_proj
  #   epsilon_pred <- data.table::as.data.table(epsilon2_new)
  #   colnames(epsilon_pred) <- c("LOC", "YEAR", "AGE", "Pred")
  #   epsilon_pred$YEAR <- as.factor(epsilon_pred$YEAR)
  #   epsilon_pred$YEAR <- as.factor(as.numeric(as.character(factor(epsilon_pred$YEAR, labels=sort(unique(data$data$YEAR))))))
  #   pred_grid$LOC <- rep(1:(nrow(pred_grid)/conf$Nyear), conf$Nyear)
  #   epsilon_pred <- epsilon_pred %>% left_join(pred_grid)
  # }
  #
  #
  output_run <- list(obj = tmb_obj_phase1, opt = opt_phase1)

  return(output_run)
}

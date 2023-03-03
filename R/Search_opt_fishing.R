##' The below function is an utility functions that tries to change the configuration of the scenario to
##' match as much as possible the catch data structure (i.e. amount of zeros, its distribution)
##' it runs a non-linear optimizer to find the best configuration thus takes some time to run
##'
##' @param x this is the parameters to be optimized to match the simulation model catch distribution to the real catch data
##' @param weight this is the weight to put between minimizing the difference in catch quantile and the amount of zeros in the data
##' @param Sim_setting  the simulation configuration setting
##' @param focus_species indicating whether the interest is on the "landed" species or the "discarded" species
##' @param catch_quantile the quantile of the catch distribution to compare between the real and simulated data
##' @param catch_data the real catch data to match with (make sure the column order match the simulated data)
##' @details
##' @return the optimized parameter values
##'

Obj_func <- function(x, weight=1, Sim_setting=Sim1, focus_species="landed", catch_quantile = seq(0.05,0.95,by=0.05),
                     catch_data, ...) {

  ### Changes to the scenario input list to optimize:
    Sim_setting$qq_original = x[1:6] * 1e-3 # the average catchability coef by species for ALL vessels to be adjusted
    Sim_setting$xi = x[7:12]                # power of the tweedie distribution to be adjusted

    Data <- Generate_scenario_data(Sim_Settings = Sim_setting, seed_input=12)
    Data$Data <- as.data.frame(Data$Data)

### optimize based on matching the distribution of the catch data to the simulated data
### In practice, it tries to match the quantiles (inputs) of the distribution by minimizing the sum across species of the absolute relative error
    if (focus_species=="landed") sel_cols = grep("Sp", colnames(Data$Data))[!grep("Sp", colnames(Data$Data)) %in% grep("_disc", colnames(Data$Data))]
    if (focus_species=="discarded") sel_cols = grep("_disc", colnames(Data$Data))

    aaa <- apply(Data$Data[,sel_cols], 2, function(x) quantile(x, catch_quantile))
    aaa_true <- apply(catch_data, 2, function(x) quantile(x, catch_quantile))

    LL <- abs(aaa-aaa_true)
    LL[which(aaa_true != 0)] <- LL[which(aaa_true != 0)]/aaa_true[which(aaa_true != 0)]
    LL <- sum(LL)

### Then adding the constraint on the amount of zero in the catch
### In practice, this constraint requires the user to use some "weight"
### so that it becomes part of the optimization procedure.
### again, we minimize here the sum across species of the absolute relative error
  bbb <- apply(Data$Data[,sel_cols], 2, function(x) sum(x==0)/length(x))
  bbb_true <- apply(catch_data, 2, function(x) sum(x==0)/length(x))

  LL1 <- abs(bbb-bbb_true)
  LL1[which(bbb_true != 0)] <- LL1[which(bbb_true != 0)]/bbb_true[which(bbb_true != 0)]
  LL1 <- sum(LL1)
  LL <-  LL + weight * LL1

  return(LL)

}


### Running now the optimization routine

catch_data= datdat %>% dplyr::select(hyse,sei,torsk,`vanlig uer`,gråsteinbit)
catch_data$gråsteinbit_m = catch_data$gråsteinbit

# Nloptr require only parameters to be optimized (or fixed). So we need to configure before hand
Obj_func1 <- function(x, weight=1){
  return(Obj_func(x, weight, Sim_setting=Sim1, focus_species="landed", catch_quantile = seq(0.05,0.95,by=0.05),
                      catch_data=catch_data))
}


library(nloptr)
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-1 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-1,
              "maxeval"= 160000,
              "local_opts" = local_opts,
              "print_level" = 0 )

res <- nloptr ( x0 = c(0.06, 0.9, 0.3, 0.2, rep(0.003 / 2, 2), 1.84,1.84,1.87,1.71,1.55,1.55),
                eval_f = Obj_func1,
                lb = c(rep(0,6), rep(1.2,6)),
                ub = c(rep(10,6),rep(2,6)),
                opts = list("algorithm"="NLOPT_LN_COBYLA",
                            "xtol_rel"=1.0e-1),
                weight = 10)

print(res)



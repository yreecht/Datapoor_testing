###### Example code to test the performance of the different models


samp_list <- seq(0.1, 1, by=0.1)
nsim = 50

  testing <- function(samps, sim){

    Final_data <- sampling_select(data = Data$Data %>% as.data.frame(), percent = samp_list[samps], unit=Sim2$samp_unit, seed=sim, months = c(11:12))

    Year_adj <- 1   # if the data is taken towards the end of the year, add the year adjustement factor

    nb_x = 3
    nb_y = 3

    Final_data <- equal_partition(Final_data, nb_x, nb_y)
    Final_data$depth_scl <- scale(Final_data$depth)

    truncate_fn <- function(x) ifelse(x < Sim2$samp_mincutoff, 0, x)
    Final_data <- Final_data %>% filter(year >= Sim2$start_year) %>% # to remove the initial year effect of the simulator ("burn-in" time)
      mutate_at(vars(starts_with("Sp")), truncate_fn) # truncate variables as in real world
    Final_data <- Final_data %>% mutate(year_fct = as.factor(year),
                                        area_fct = as.factor(area),
                                        month_fct = as.factor(month),
                                        vessel_fct = as.factor(vessel))


    Final_data_df <- Final_data %>% pivot_longer(cols = starts_with("Sp"), names_to = "Species", values_to = "CPUE_trunc")


    #### projection grid
    qcs_grid <- do.call(rbind, replicate(Sim2$n_years, Data$bathym, simplify=FALSE))
    qcs_grid$year <- as.numeric(rep(1:Sim2$n_years, each=nrow(Data$bathym)))
    qcs_grid$vessel <- unique(Final_data_df$vessel)[1]
    qcs_grid$vessel_fct <- unique(Final_data_df$vessel)[1]
    qcs_grid <-  qcs_grid %>% mutate(X = as.numeric(X), Y= as.numeric(Y), depth_scl = (depth-mean(qcs_grid$depth))/sd(qcs_grid$depth))
    qcs_grid$CPUE_trunc <- 1
    qcs_grid <- equal_partition(qcs_grid, nb_x, nb_y)
    qcs_grid$year_area_fct <- as.factor(apply(qcs_grid[, c('year','area')], 1, function(x) paste(x, collapse="_")))
    qcs_grid <- qcs_grid %>% filter(year >= Sim2$start_year)

    # now adding the average species catch in each year, area combination
    qcs_grid <- qcs_grid %>% left_join(
      Final_data %>% pivot_longer(cols = starts_with("Sp"), names_to = "Species", values_to = "CPUE") %>%
        group_by(year, area, Species) %>% summarize(Mean = mean(CPUE)) %>% ungroup() %>%
        pivot_wider(names_from = "Species", values_from = "Mean"))
    qcs_grid <- qcs_grid %>% mutate(area_fct = factor(area), year_fct = as.factor(year), yy=1)



    ##### generating single species data (for single species models)
    Final_data_df$vessel_fct <- as.factor(Final_data_df$vessel)

    out <- c()
    for (Which_sp in c("Sp1","Sp2","Sp3","Sp4","Sp5")){

      Final_data_bycatch = Final_data_df %>% filter(Species == Which_sp)
      Final_data_bycatch$year_fct <- as.factor(Final_data_bycatch$year)
      Final_data_bycatch <- Final_data_bycatch %>% mutate(depth_scl = (depth-mean(qcs_grid$depth))/sd(qcs_grid$depth))
      Final_data_bycatch$year_area_fct <- as.factor(apply(Final_data_bycatch[, c('year','area')], 1, function(x) paste(x, collapse="_")))
      Final_data_bycatch <- Final_data_bycatch %>% mutate(area = droplevels(area))

      ####### gam:
      library(mgcv)
      gam1 <- try(gam(CPUE_trunc ~ 0 + as.factor(year) + as.factor(area) +  s(depth_scl) +
                    s(vessel_fct, bs="re"), data=Final_data_bycatch, family = tw))
      if (TRUE %in% class(gam1) == "try-error" | !(all(levels(qcs_grid$area) %in% levels(Final_data_bycatch$area)))) {
        IA_gam <- c()
      } else {
        predicted <- predict(gam1, newdata = qcs_grid, type="response", se.fit=F)
        pred <- qcs_grid
        pred$pred = predicted
        pred <- pred %>% mutate(year_or = year, year = year_or + Year_adj)
        IA_gam <- pred %>% group_by(year) %>% summarize(IA = sum(pred, na.rm=T)) %>% mutate(Method="GAM", IA = IA/IA[1])

        IA_true <- data.frame(IA_true=true_index[-c(1:(Sim2$start_year-1)),Which_sp]/true_index[Sim2$start_year+Year_adj,Which_sp], year = Sim2$start_year:Sim2$n_years)
        IA_gam <- IA_gam %>% left_join(IA_true) %>% mutate(RE = IA- IA_true)

        IA_gam$samps <- samp_list[samps]
        IA_gam$iter <- sim
        IA_gam$Species = Which_sp
      }
      out <- rbind(out, IA_gam)
    }

    return(out)
  }


  scen <- expand.grid(samps = 1:length(samp_list), sim = 1:nsim)

  library(future)
  library(furrr)
  plan(multisession, workers = 2)
  OUT <- c()
  OUT <- furrr::future_pmap_dfr(scen, testing, .options = furrr::furrr_options(seed = TRUE))

  failure <- OUT %>% group_by(samps, Species) %>% distinct(iter) %>% summarize(n=n())
  MRE_out <- OUT %>% group_by(samps, year, Species) %>% summarize(mean_RE = mean(RE), lwr95=quantile(RE, 0.025, na.rm=T), upp95=quantile(RE, 0.975, na.rm=T))

  ggplot(MRE_out, aes(x=year, y=mean_RE, col=as.factor(samps))) + facet_grid(~Species) + geom_line(size=1.5) + theme_bw() +
    geom_ribbon(aes(ymin=lwr95, ymax=upp95, fill=as.factor(samps))) +
    geom_hline(yintercept=0) + scale_color_viridis_d(option = "magma") + scale_fill_viridis_d(option = "magma")

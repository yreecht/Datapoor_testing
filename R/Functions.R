##' Function to sample from the whole fishing fleet:
##' @param data is the whole fleet data that was generated
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
##' @return the final sample data to be fed into the CPUE standardization exercise
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
equal_partition <- function(data, ncut_x=2, ncut_y=2){
  out <- expand.grid(X = Sim2$Range_X, Y = Sim2$Range_Y)
  out <- out %>% mutate(X_cut = factor(cut(X, ncut_x), labels=c(1:ncut_x)),
                        Y_cut = factor(cut(Y, ncut_y), labels=c(1:ncut_y)))
  area <- as.factor(apply(out[,3:4], 1, function(x) paste(x, collapse="_")))
  area <- factor(area, labels = 1:(ncut_x*ncut_y))
  area <- data.frame(X=out$X, Y=out$Y, area)

  data <- data %>% left_join(area)
  return(data)
}





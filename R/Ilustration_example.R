set.seed(12)
model_bathym <- RMgauss(var=150^2, scale=25)
map_grid <- expand.grid(X=1:40, Y=1:40)
Bathym <- RFsimulate(model_bathym, x=map_grid)
data.bathym1 <- data.frame(ID = 1:nrow(map_grid), map_grid , depth=Bathym@data[,1])
data.bathym$depth = data.bathym$depth+abs(min(data.bathym$depth))	# to avoid getting negative depth here depth takes positive real number
image.plot(Sim_Settings$Range_X, Sim_Settings$Range_Y, matrix(data.bathym$depth,nrow=length(Sim_Settings$Range_X), ncol=length(Sim_Settings$Range_Y)), ylab="", xlab="")

est_model <- RMstable(alpha=2, var=NA, scale=NA)
Bathym <- RFsimulate(model_bathym, x=map_grid)
data.bathym <- data.frame(map_grid , depth=Bathym@data[,1])
hist(data.bathym$depth)
bla <- RFfit(est_model, data=data.bathym)


image.plot(1:40, 1:40, matrix(Bathym@data[,1],nrow=40, ncol=40), ylab="", xlab="")

coordinates(data.bathym) <- ~ X + Y
library(automap)
m <- autofitVariogram(depth~1, data.bathym)
plot(m)
library(gstat)
v <- variogram(depth~1, data.bathym)#create a variogram of the sorting data
m <- fit.variogram(v, vgm(psill=m$var_model[2,2], model=as.character(m$var_model[2,1]), range=m$var_model[2,3], nugget =m$var_model[1,2], kappa=m$var_model[2,4]))    #fit a model to the variogram
plot(v, model= m)
g <- gstat(id = "depth", formula = depth~1, data=data.bathym, model = m, nmax=5)
library(raster)
map_raster <- rasterFromXYZ(data.bathym, crs="")
vals <- raster::interpolate(object=map_raster, model=g, xyOnly=TRUE, progress="text") #Interpolate the object to a raster
xyz <- rasterToPoints(vals)


### Tweedie parameters effects
  par(mfrow=c(2,2))
  set.seed(1)
  hist(rtweedie(10000, mu= 0.05, power=1.7, phi=0.2), main = "power=1.7, phi=0.2")
  set.seed(1)
  hist(rtweedie(10000, mu= 0.05, power=1.2, phi=0.2), main = "power=1.3, phi=0.2")
  set.seed(1)
  hist(rtweedie(10000, mu= 0.05, power=1.5, phi=0.2), main = "power=1.5, phi=0.2")
  set.seed(1)
  hist(rtweedie(10000, mu= 0.05, power=1.7, phi=0.4), main = "power=1.7, phi=0.4")
  set.seed(1)
  hist(rtweedie(10000, mu= 0.05, power=1.9, phi=0.2), main = "power=1.7, phi=0.4")
  set.seed(1)
  hist(rtweedie(10000, mu= 0.05, power=1.9, phi=0.5), main = "power=1.9, phi=0.5")


#### Illustration movement probability
  Par_mvt_adult[[1]][1,1] = 1
  bla <- Mvt_upd1(loc = 40*40-39, Pop_param=Par_mvt_adult[[1]][,1], Depth_eff="TRUE", Dist_eff="TRUE", Lat_eff="TRUE", Dist_func=Sim_Settings$func_mvt_dist[xxx], Depth_func=Sim_Settings$func_mvt_depth[1], Lat_func=Sim_Settings$func_mvt_lat[1], data.bathym)
  image.plot(matrix(bla, nrow=40, ncol=40))
  Par_mvt_adult[[1]][1,1] = 5
  bla <- Mvt_upd1(loc = 10, Pop_param=Par_mvt_adult[[1]][,1], Depth_eff="TRUE", Dist_eff="TRUE", Lat_eff="TRUE", Dist_func=Sim_Settings$func_mvt_dist[xxx], Depth_func=Sim_Settings$func_mvt_depth[1], Lat_func=Sim_Settings$func_mvt_lat[1], data.bathym)
  image.plot(matrix(bla, nrow=40, ncol=40))
  Par_mvt_adult[[1]][1,1] = 100
  bla <- Mvt_upd1(loc = 10, Pop_param=Par_mvt_adult[[1]][,1], Depth_eff="TRUE", Dist_eff="TRUE", Lat_eff="TRUE", Dist_func=Sim_Settings$func_mvt_dist[xxx], Depth_func=Sim_Settings$func_mvt_depth[1], Lat_func=Sim_Settings$func_mvt_lat[1], data.bathym)
  image.plot(matrix(bla, nrow=40, ncol=40))

  plot(1:600, dlnorm(1:600, log(250)-1^2/2, 1), type="l", main = "log normal", ylab="", xlab="", col="red",lwd=2)
  plot(1:50, dexp(1:50, rate = 1/1), type="l", main = "Exponential", ylab="", xlab="", col="red",lwd=2)
  plot(1:50, dexp(1:50, rate = 1/5), type="l", main = "Exponential", ylab="", xlab="", col="red",lwd=2)
  plot(1:50, dexp(1:50, rate = 1/100), type="l", main = "Exponential", ylab="", xlab="", col="red",lwd=2)

#### Habitat preference function

  image.plot(Sim_Settings$Range_X, Sim_Settings$Range_Y, matrix(data.bathym$depth,nrow=length(Sim_Settings$Range_X), ncol=length(Sim_Settings$Range_Y)))

  par(mfrow=c(2,2), mar=c(2,3,2,1), oma=c(1,1,1,1))
  plot(1:600, dnorm(1:600, 250, 50), type="l", main="Gaussian", ylab="", xlab="", col="red",lwd=2)
  plot(1:600, dlnorm(1:600, log(250)-0.5^2/2, 0.5), type="l", main = "log normal", ylab="", xlab="", col="red",lwd=2)
  plot(1:600, dexp(1:600, rate = 1/200), type="l", main = "Exponential", ylab="", xlab="", col="red",lwd=2)
  plot(1:600, dunif(1:600, 100, 400), type="l", main = "Uniform", ylab="", xlab="", col="red",lwd=2)

#### Fishing location
  ggplot(data = Data$Data %>% filter(year == 1)) + geom_raster(data=Data$bathym, aes(x=X, y=Y, fill=depth)) + scale_fill_viridis_c() +
    geom_point(aes(x=X, y=Y)) + facet_wrap(.~month)






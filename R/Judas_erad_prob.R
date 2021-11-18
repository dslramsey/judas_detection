
library(terra)
library(simecol)
library(sf)
library(ggspatial)
library(gridExtra)
library(tidyverse)
library(lubridate)
library(posterior)
source("R/judas functions.r")

goat_ids<- read_csv("data/Santiago_judas_ID.csv")
goat_ids<- goat_ids %>% mutate(ID = factor(as.character(ID)))

out<- readRDS("out/goat_surv.rds") # generate from discrete_survival_analysis.r 
fit<- as_draws_rvars(out$draws(c("beta0","beta")))

beta0<- tibble(mean=E(fit$beta0), sd=sd(fit$beta0))
beta1<- tibble(mean=E(fit$beta)[1], sd=sd(fit$beta)[1])
beta2<- tibble(mean=E(fit$beta)[-1], sd=sd(fit$beta)[-1])
sigma<- tibble(mean=goat_chr$sigma, sd=goat_chr$se)
X<- model.matrix(~Sex + Trans, data=goat_ids)[,-1]

#----------------------------------------------
# Goats
# Santiago Island - Galapagos 

temp <- tempfile()
unzip("Data/santiago_island.zip", exdir = temp)
shp_files<-list.files(temp, pattern = ".shp$",full.names=TRUE)

santiago<- read_sf(shp_files) #CRS is UTM Zone 16 (units are km)

hr_center<- goat_ids %>% mutate(Time = as.numeric(floor((date_max - date_min)/28) + 1))
nperiod<- hr_center$Time

timeseq<- seq(as.Date("2004-6-01"),as.Date("2005-12-01"),length.out=4)

results<- matrix(NA,nrow=length(timeseq)-1,ncol=13)

for(i in 1:(length(timeseq)-1)) {
  inds<- which(hr_center$date_min < timeseq[i+1] & hr_center$date_max > timeseq[i])
  dat<- slice(hr_center, inds)
  
  tottime<- rep(NA,length(inds))
  for(z in 1:length(inds)) {
    tottime[z]<- intersect.period(c(timeseq[i],timeseq[i+1]),
                                c(hr_center$date_min[inds[z]],hr_center$date_max[inds[z]]))
  }
  tottime<- tottime/28
  tottime<- floor(tottime)+1

  nsims<- 50
  cellsize <- 0.5
  tmpsims<- matrix(NA,nsims,4)
  n<- nrow(dat)

  if(i==1) prior<- runif(nsims, 0, 1)
  else prior<- prior_updated
  
  for(j in 1:nsims) {
    cat("doing bootstrap sample ",j," for period ",i,"\n")
    b1<- rnorm(n, beta0$mean[inds], beta0$sd[inds])
    b2<- rnorm(n, beta1$mean, beta1$sd)
    ce<- rcov_effects(beta2, X[inds,])
    params<- data.frame(b1=b1,b2=b2,ce=ce,sigma=sigma$mean[inds])
    bvn.surf<- make.surface(dat,params,santiago,nperiod=tottime,cellsize=cellsize,CE=0.95,Pu=1,prior=prior[j])
    tmpsims[j,]<- unlist(bvn.surf$Table)
  }
  results[i,1:3]<- c(mean(tmpsims[,1]),quantile(tmpsims[,1],c(0.05,0.95)))
  results[i,4:6]<- c(mean(tmpsims[,2]),quantile(tmpsims[,2],c(0.05,0.95)))
  results[i,7:9]<- c(mean(tmpsims[,3]),quantile(tmpsims[,3],c(0.05,0.95)))
  results[i,10:12]<- c(mean(tmpsims[,4]),quantile(tmpsims[,4],c(0.05,0.95)))
  results[i,13]<- n
  prior_updated<- as.vector(tmpsims[,4])
}
results<- data.frame(results)
names(results)<- c("cellsens","cellsens.l","cellsens.u","Cov","Cov.l","Cov.u","SSe","SSe.l","SSe.u",
                   "PE","PE.l","PE.u","n")

results
##------------------------
# Plot single realisation (takes a while due to smaller cell resolution (100m))

goat_locs<- read_csv("data/santiago_judas_locs.csv")
goat_locs<- goat_locs %>% mutate(ID = factor(as.character(ID)))

rast_map<- list()
Pd<- list()

cellsize<- 0.5 # decrease cellsize for better resolution (takes a while to process at cell resolution 0.1)
#cellsize<- 0.1

for(i in 1:(length(timeseq)-1)) {
  cat("doing period ",i,"\n")
  inds<- which(hr_center$date_min < timeseq[i+1] & hr_center$date_max > timeseq[i])
  dat<- slice(hr_center, inds)
  
  tottime<- rep(NA,length(inds))
  for(z in 1:length(inds)) {
    tottime[z]<- intersect.period(c(timeseq[i],timeseq[i+1]),
                                  c(hr_center$date_min[inds[z]],hr_center$date_max[inds[z]]))
  }
  tottime<- tottime/28
  tottime<- floor(tottime)+1
  
  if(i==1) prior<- runif(1)
  else prior<- prior_updated
  
  n<- nrow(dat)
  b1<- rnorm(n, beta0$mean[inds], beta0$sd[inds])
  b2<- rnorm(n, beta1$mean, beta1$sd)
  ce<- rcov_effects(beta2, X[inds,])
  params<- data.frame(b1=b1,b2=b2,ce=ce,sigma=sigma$mean[inds])
  bvn.surf<- make.surface(dat,params,santiago, nperiod=tottime, cellsize=cellsize, CE=0.95, Pu=1, prior=prior)
  rast_map[[i]]<- bvn.surf$Raster
  Pd[[i]]<- bvn.surf$Table
  prior_updated<- bvn.surf$Table$PE
}


# Now plot
plt<- list()
labs<- c("Jun-Dec 2004","Jan-Jun 2005","Jul-Dec 2005")
for(i in 1:(length(timeseq)-1)) {
  inds<- which(hr_center$date_min < timeseq[i+1] & hr_center$date_max > timeseq[i])
  dat<- slice(hr_center, inds)
  locs<- goat_locs %>% filter(ID %in% dat$ID)
  pp<- as.data.frame(rast_map[[i]], xy=TRUE)

plt[[i]]<- santiago %>% ggplot() +
  annotation_scale(location = "bl", width_hint = 0.25,height = unit(0.5,"cm"),text_cex=1) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.15, "in"), pad_y = unit(0.15, "in"),
                         style = north_arrow_fancy_orienteering) +
  geom_raster(aes(x=x, y=y, fill=lyr.1), data=pp, interpolate=F) + 
  scale_fill_distiller(palette = "Spectral", na.value = "transparent",limits=c(0,1)) +
  geom_point(aes(xbar, ybar), data=dat, shape=16, size=1, color=adjustcolor("black", alpha.f = 0.9)) +
  geom_point(aes(East,North), data=locs, size=0.5, color="grey70") +
  labs(x="Easting",y="Northing",fill="Probability") +
  geom_sf(fill=NA) +
  labs(title=paste0("Period ",i,": ",labs[i])) +
  theme_bw()
}

win.graph(7,12)
grid.arrange(plt[[1]],plt[[2]],plt[[3]], ncol=1)

#-----------------------------------------
# Santa Cruz pigs
#-----------------------------------------

pig_ids<- read_csv("data/Santacruz_judas_ID.csv")
pig_ids<- pig_ids %>% mutate(ID = factor(as.character(ID)))


out<- readRDS("out/pig_surv.rds") # from discrete_survival_analysis
fit<- as_draws_rvars(out$draws(c("beta0","beta")))

beta0<- tibble(mean=E(fit$beta0), sd=sd(fit$beta0))
beta1<- tibble(mean=E(fit$beta)[1], sd=sd(fit$beta)[1])
beta2<- tibble(mean=E(fit$beta)[-1], sd=sd(fit$beta)[-1])
sigma<- tibble(mean=pig_chr$sigma, sd=pig_chr$se)
X<- model.matrix(~ Sex, pig_ids)[,-1]

# Read in Island shapefile

temp <- tempfile()
unzip("Data/santacruz_island_zones.zip", exdir = temp)
shp_files<-list.files(temp, pattern = ".shp$",full.names=TRUE)

santacruz<- read_sf(shp_files) #CRS is California Albers, epsg:3310 (units are km)
zone1<- filter(santacruz, Zone==1) # Extract zone 1

hr_center<- pig_ids %>% mutate(Time = as.numeric(floor((date_max - date_min)/28) + 1))
nperiod<- hr_center$Time

timeseq<- seq(as.Date("2005-05-17"),as.Date("2006-08-01"),length.out=4)

results<- matrix(NA,nrow=length(timeseq)-1,ncol=13)

for(i in 1:(length(timeseq)-1)) {
  inds<- which(hr_center$date_min < timeseq[i+1] & hr_center$date_max > timeseq[i])
  dat<- slice(hr_center, inds)
  
  tottime<- rep(NA,length(inds))
  for(z in 1:length(inds)) {
    tottime[z]<- intersect.period(c(timeseq[i],timeseq[i+1]),
                                  c(hr_center$date_min[inds[z]],hr_center$date_max[inds[z]]))
  }
  tottime<- tottime/28
  tottime<- floor(tottime)+1
  
  nsims<- 50
  cellsize <- 0.1
  tmpsims<- matrix(NA,nsims,4)
  n<- nrow(dat)
  
  if(i==1) prior<- runif(nsims, 0, 1)
  else prior<- prior_updated
  
  for(j in 1:nsims) {
    cat(paste("doing bootstrap sample ",j," for period ",i,sep=""),"\n")
    b1<- rnorm(n, beta0$mean[inds], beta0$sd[inds])
    b2<- rnorm(n, beta1$mean, beta1$sd)
    params<- data.frame(b1=b1,b2=b2,sigma=sigma$mean[inds])
    bvn.surf<- make.surface(dat,params,zone1,nperiod=tottime,cellsize=cellsize,CE=0.95,Pu=1,prior=prior[j])
    tmpsims[j,]<- unlist(bvn.surf$Table)
  }
  results[i,1:3]<- c(mean(tmpsims[,1]),quantile(tmpsims[,1],c(0.05,0.95)))
  results[i,4:6]<- c(mean(tmpsims[,2]),quantile(tmpsims[,2],c(0.05,0.95)))
  results[i,7:9]<- c(mean(tmpsims[,3]),quantile(tmpsims[,3],c(0.05,0.95)))
  results[i,10:12]<- c(mean(tmpsims[,4]),quantile(tmpsims[,4],c(0.05,0.95)))
  results[i,13]<- n
  prior_updated<- as.vector(tmpsims[,4])
}

results<- data.frame(results)
names(results)<- c("cellsens","cellsens.l","cellsens.u","Cov","Cov.l","Cov.u","SSe","SSe.l","SSe.u",
                   "PE","PE.l","PE.u","n")


results

##------------------------
# Single realisation for plotting

pig_locs<- read_csv("data/santacruz_judas_locs.csv")
pig_locs<- pig_locs %>% mutate(ID = factor(as.character(ID)),East=East/1000, North=North/1000)


rast_map<- list()
Pd<- list()

cellsize<- 0.1 # 100m

for(i in 1:(length(timeseq)-1)) {
  inds<- which(hr_center$date_min < timeseq[i+1] & hr_center$date_max > timeseq[i])
  dat<- slice(hr_center, inds)
  
  tottime<- rep(NA,length(inds))
  for(z in 1:length(inds)) {
    tottime[z]<- intersect.period(c(timeseq[i],timeseq[i+1]),
                                  c(hr_center$date_min[inds[z]],hr_center$date_max[inds[z]]))
  }
  tottime<- tottime/28
  tottime<- floor(tottime)+1
  
  if(i==1) prior<- runif(1)
  else prior<- prior_updated
  
  n<- nrow(dat)
  b1<- rnorm(n, beta0$mean[inds], beta0$sd[inds])
  b2<- rnorm(n, beta1$mean, beta1$sd)
  params<- data.frame(b1=b1,b2=b2,sigma=sigma$mean[inds])
  bvn.surf<- make.surface(dat,params,zone1, nperiod=tottime, cellsize=cellsize, CE=0.95, Pu=1, prior=prior)
  rast_map[[i]]<- bvn.surf$Raster
  Pd[[i]]<- bvn.surf$Table
  prior_updated<- bvn.surf$Table$PE
}

# Now plot
plt<- list()
labs<- c("Jun-Oct 2005","Nov 2005-Mar 2006","Apr-Aug 2006")
for(i in 1:(length(timeseq)-1)) {
  inds<- which(hr_center$date_min < timeseq[i+1] & hr_center$date_max > timeseq[i])
  dat<- slice(hr_center, inds)
  locs<- pig_locs %>% filter(ID %in% dat$ID)
  pp<- as.data.frame(rast_map[[i]], xy=TRUE)
  
  plt[[i]]<- zone1 %>% ggplot() +
    annotation_scale(location = "bl", width_hint = 0.25,height = unit(0.5,"cm"),text_cex=1) +
    annotation_north_arrow(location = "tr", which_north = "true", 
                           pad_x = unit(0.05, "in"), pad_y = unit(0.05, "in"),
                           style = north_arrow_fancy_orienteering) +
    geom_raster(aes(x=x, y=y, fill=lyr.1), data=pp, interpolate=F) + 
    scale_fill_distiller(palette = "Spectral", na.value = "transparent",limits=c(0,1)) +
    geom_point(aes(xbar, ybar), data=dat, shape=16, size=1, color=adjustcolor("black", alpha.f = 0.9)) +
    geom_point(aes(East,North), data=locs, size=0.5, color="grey70") +
    labs(x="Easting",y="Northing",fill="Probability") +
    geom_sf(fill=NA) +
    labs(title=paste0("Period ",i,": ",labs[i])) +
    theme_bw()
  
}

win.graph(7,12)
grid.arrange(plt[[1]],plt[[2]],plt[[3]], ncol=1)

#---------------------------------
library(cpp11)

cpp_source("src/neighbor_matrix.cpp")

tmpmat<- matrix(0, 30,30)
tmpmat[15,20]<- 1

wkern<- matrix(0.5, 5,5)

neighbourhood(tmpmat, wdist=wkern, state=1)


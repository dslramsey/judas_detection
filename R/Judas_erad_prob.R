
#----------------------------------------------------------------------
# Calculate probability of eradication in area (shapefile), given 
# sample of Judas animals detect no 'wild' individuals
#
# SSe - Surveillance sensitivity (sensu Anderson et al 2013) aka 
# probability of detection in area (adjusted for coverage)
#
# Cov - coverage of area by combined Judas animals (proportion of cells
# with > 0 SSe)
#
# PE - probability of eradication, given SSE and prior probability
#
# D. Ramsey (22/11/2021)
#----------------------------------------------------------------------

library(terra)
library(simecol)
library(sf)
library(ggspatial)
library(gridExtra)
library(tidyverse)
library(purrr)
library(lubridate)
library(posterior)
source("R/judas functions.r")

#----------------------------------------------
# Goats
# Santiago Island - Galapagos 
#----------------------------------------------

# read in Judas IDs containing info on deployment period (date_min, date_max)
goat_ids<- read_csv("data/Santiago_judas_ID.csv")
goat_ids<- goat_ids %>% mutate(ID = factor(as.character(ID)))

# read in parameter estimates (generate from discrete_survival_analysis.r) 
out<- readRDS("out/goat_surv.rds") 
fit<- as_draws_rvars(out$draws(c("beta0","beta")))

# Calculate home range scale (sigma) from Judas locations 
goat_locs<- read_csv("data/santiago_judas_locs.csv")
goat_locs<- goat_locs %>% mutate(ID = factor(as.character(ID)))
goat_chr<- goat_locs %>% group_by(ID) %>% nest() %>% arrange(ID)
goat_chr<- goat_chr %>% mutate(CNorm=map(data, calc_hr))
goat_chr<- goat_chr %>% unnest(CNorm) %>% select(-data)

# Assemble parameter estimates
beta0<- tibble(mean=E(fit$beta0), sd=sd(fit$beta0)) # Intercept
beta1<- tibble(mean=E(fit$beta)[1], sd=sd(fit$beta)[1]) # distance effect
beta2<- tibble(mean=E(fit$beta)[-1], sd=sd(fit$beta)[-1]) #sex/status
sigma<- tibble(mean=goat_chr$sigma, sd=goat_chr$se) # Home range scale
X<- model.matrix(~Sex + Trans, data=goat_ids)[,-1] # sex/status covariates
pars<- list(beta0=beta0,beta1=beta1,beta2=beta2,sigma=sigma,X=X)

# Read in shapefile of the island
temp <- tempfile()
unzip("Data/santiago_island.zip", exdir = temp)
shp_files<-list.files(temp, pattern = ".shp$",full.names=TRUE)

santiago<- read_sf(shp_files) #CRS is UTM Zone 16 (units are km)

# Function to simulate posterior distributions of SSe, Cov and PE
posterior_sims<- function(prior, ii, df, pars, shape, np, cs, PU) {
  n<- nrow(df)
  b0<- rnorm(n, pars$beta0$mean[ii], pars$beta0$sd[ii])
  b1<- rnorm(n, pars$beta1$mean, pars$beta1$sd)
  ce<- rcov_effects(pars$beta2, pars$X[ii,])
  sigma<- pars$sigma$mean[ii]
  params<- data.frame(b0=b0,b1=b1,ce=ce,sigma=sigma)
  bvn.surf<- make.surface(df,params,shape,nperiod=np,cellsize=cs,Pu=PU,prior=prior)
  vals<- bvn.surf$Table
  vals
}


# Create periods
timeseq<- seq(as.Date("2004-6-01"),as.Date("2005-12-01"),length.out=4)

nsims<- 50  # No. of samples. increase for serious inference

cellsize <- 0.5 # grid cell resolution (km). Smaller cells take longer to run

sim_summary<- list()  # for collecting intermediate results

# The following code finds Judas that were present at some point in each period 
for(i in 1:(length(timeseq)-1)) {
  cat("Doing period ",i,"\n")
  inds<- which(goat_ids$date_min < timeseq[i+1] & goat_ids$date_max > timeseq[i])
  dat<- slice(goat_ids, inds)
  
  tottime<- rep(NA,length(inds))
  for(z in 1:length(inds)) {
    tottime[z]<- intersect.period(c(timeseq[i],timeseq[i+1]),
                                c(goat_ids$date_min[inds[z]],goat_ids$date_max[inds[z]]))
  }
  tottime<- tottime/28  # Months
  tottime<- floor(tottime)+1
  
  # Initialise or update prior probability of eradication
  if(i==1) prior<- runif(nsims, 0, 1) else prior<- prior_updated
  
  sim_out<- prior %>% map_dfr(posterior_sims, ii=inds, df=dat, pars=pars, shape=santiago,
                              np=tottime, cs=cellsize, PU=1)
   
  sim_summary[[i]]<- imap_dfr(sim_out, ~summarise_sims(.y, .x, period=i))
  
  prior_updated<- as.vector(sim_out$PE)
}

sim_summary<- map_dfr(sim_summary, rbind)
sim_summary


##------------------------
# Plot single realisation 

plt<- list() # Collect plots

labs<- c("Jun-Dec 2004","Jan-Jun 2005","Jul-Dec 2005")

cellsize<- 0.5 # decrease cellsize for better resolution (takes a while to process at cell resolution 0.1)
#cellsize<- 0.1

for(i in 1:(length(timeseq)-1)) {
  cat("doing period ",i,"\n")
  inds<- which(goat_ids$date_min < timeseq[i+1] & goat_ids$date_max > timeseq[i])
  dat<- slice(goat_ids, inds)
  locs<- goat_locs %>% filter(ID %in% dat$ID)
  
  tottime<- rep(NA,length(inds))
  for(z in 1:length(inds)) {
    tottime[z]<- intersect.period(c(timeseq[i],timeseq[i+1]),
                                  c(goat_ids$date_min[inds[z]],goat_ids$date_max[inds[z]]))
  }
  tottime<- tottime/28
  tottime<- floor(tottime)+1
  
  n<- nrow(dat)
  b0<- rnorm(n, beta0$mean[inds], beta0$sd[inds])
  b1<- rnorm(n, beta1$mean, beta1$sd)
  ce<- rcov_effects(beta2, X[inds,])
  params<- data.frame(b0=b0,b1=b1,ce=ce,sigma=sigma$mean[inds])
  bvn.surf<- make.surface(dat,params,santiago, nperiod=tottime, cellsize=cellsize, Pu=1, prior=0.5)
  
  pp<- as.data.frame(bvn.surf$Raster, xy=TRUE)

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

# Calculate home range scale (sigma) from Judas locations 
pig_locs<- read_csv("data/santacruz_judas_locs.csv")
pig_locs<- pig_locs %>% mutate(ID = factor(as.character(ID)),East=East/1000, North=North/1000)
pig_chr<- pig_locs %>% group_by(ID) %>% nest() %>% arrange(ID)
pig_chr<- pig_chr %>% mutate(CNorm=map(data, calc_hr))
pig_chr<- pig_chr %>% unnest(CNorm) %>% select(-data)

# Assemble parameter estimates
beta0<- tibble(mean=E(fit$beta0), sd=sd(fit$beta0))
beta1<- tibble(mean=E(fit$beta)[1], sd=sd(fit$beta)[1])
beta2<- tibble(mean=E(fit$beta)[-1], sd=sd(fit$beta)[-1])
sigma<- tibble(mean=pig_chr$sigma, sd=pig_chr$se)
X<- as.matrix(model.matrix(~ Sex, pig_ids)[,-1])
pars<- list(beta0=beta0,beta1=beta1,beta2=beta2,sigma=sigma,X=X)

# Read in Island shapefile

temp <- tempfile()
unzip("Data/santacruz_island_zones.zip", exdir = temp)
shp_files<-list.files(temp, pattern = ".shp$",full.names=TRUE)

santacruz<- read_sf(shp_files) #CRS is California Albers, epsg:3310 (units are km)
zone1<- filter(santacruz, Zone==1) # Extract zone 1


timeseq<- seq(as.Date("2005-05-17"),as.Date("2006-08-01"),length.out=4)

nsims<- 50  # No. of samples. increase for serious inference

cellsize <- 0.5 # grid cell resolution (km). Smaller cells take longer to run

sim_summary<- list()  # for collecting intermediate results

# The following code finds Judas that were present at some point in each period 
for(i in 1:(length(timeseq)-1)) {
  cat("Doing period ",i,"\n")
  inds<- which(pig_ids$date_min < timeseq[i+1] & pig_ids$date_max > timeseq[i])
  dat<- slice(pig_ids, inds)
  
  tottime<- rep(NA,length(inds))
  for(z in 1:length(inds)) {
    tottime[z]<- intersect.period(c(timeseq[i],timeseq[i+1]),
                                  c(pig_ids$date_min[inds[z]],pig_ids$date_max[inds[z]]))
  }
  tottime<- tottime/28  # Months
  tottime<- floor(tottime)+1
  
  # Initialise or update prior probability of eradication
  if(i==1) prior<- runif(nsims, 0, 1) else prior<- prior_updated
  
  sim_out<- prior %>% map_dfr(posterior_sims, ii=inds, df=dat, pars=pars, shape=zone1,
                              np=tottime, cs=cellsize, PU=1)
  
  sim_summary[[i]]<- imap_dfr(sim_out, ~summarise_sims(.y, .x, period=i))
  
  prior_updated<- as.vector(sim_out$PE)
}

sim_summary<- map_dfr(sim_summary, rbind)
sim_summary


##------------------------
# Single realisation for plotting


plt<- list()

labs<- c("Jun-Oct 2005","Nov 2005-Mar 2006","Apr-Aug 2006")

cellsize<- 0.1 # 100m

for(i in 1:(length(timeseq)-1)) {
  cat("doing period ",i,"\n")
  inds<- which(pig_ids$date_min < timeseq[i+1] & pig_ids$date_max > timeseq[i])
  dat<- slice(pig_ids, inds)
  locs<- pig_locs %>% filter(ID %in% dat$ID)
  
  tottime<- rep(NA,length(inds))
  for(z in 1:length(inds)) {
    tottime[z]<- intersect.period(c(timeseq[i],timeseq[i+1]),
                                  c(pig_ids$date_min[inds[z]],pig_ids$date_max[inds[z]]))
  }
  tottime<- tottime/28
  tottime<- floor(tottime)+1
  
  n<- nrow(dat)
  b0<- rnorm(n, beta0$mean[inds], beta0$sd[inds])
  b1<- rnorm(n, beta1$mean, beta1$sd)
  ce<- rcov_effects(beta2, X[inds,])
  params<- data.frame(b0=b0,b1=b1,ce=ce,sigma=sigma$mean[inds])
  bvn.surf<- make.surface(dat, params, zone1, nperiod=tottime, cellsize=cellsize, Pu=1, prior=0.5)
  
  pp<- as.data.frame(bvn.surf$Raster, xy=TRUE)
  
  plt[[i]]<- zone1 %>% ggplot() +
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




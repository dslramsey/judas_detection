#-------------------------------------------------------------
# Functions for Judas interaction analysis
#-------------------------------------------------------------


# Circular normal HR function
calc_hr<- function(df) {
  bnorm.obj<- function(par, df) {
    # circular normal log likelihood function
    e<- dnorm(df$xdev, 0, par, log=TRUE) + dnorm(df$ydev, 0, par, log=TRUE)
    -sum(e)
  }
  df<- df %>% mutate(X= East, Y = North,
                     xdev = X-mean(X), ydev = Y-mean(Y)) 
  par.init<- mean(c(abs(df$xdev),abs(df$ydev)))
  par.max<- max(c(abs(df$xdev),abs(df$ydev)))
  mle<- optim(par.init, bnorm.obj, lower=0, upper=par.max, method="Brent", hessian=TRUE, df=df)
  se<- as.numeric(sqrt(solve(mle$hessian)))
  tibble(sigma=mle$par,se=se)
}

calc.latlong.dist<- function(xy1,xy2)
{
  # uses spherical law of cosines to calculate distance between two lat/long
  # coordinates in decimal degrees
  R<- 6371 # Earths radius
  xy1<- (pi * xy1)/180 # radians
  xy2<- (pi * xy2)/180 
  D<- acos(sin(xy1[,1])*sin(xy2[,1]) + cos(xy1[,1])*cos(xy2[,1])*cos(xy2[,2]-xy1[,2]))  
  return(R*D)
}
#--------------------------------------------------
calc.dist<-  function(xy1, xy2)
{ # calculate Euclidean distance from xy1 to xy2
  deltax<- (xy1[,1] - xy2[,1])^2
  deltay<- (xy1[,2] - xy2[,2])^2
  return(sqrt(deltax + deltay))
}
#----------------------------------------------------

rescale_shp<- function(shp, units="m") {
  # rescale sf object to different units (e.g. m -> km)
  p4str<- st_crs(shp)$proj4string
  bits<- strsplit(p4str,"+", fixed=TRUE)[[1]]
  ind<- pmatch("units", bits)
  bits[ind]<- paste0("units=",units," ")
  newstr<- paste(bits, collapse="+")
  shp<- st_transform(shp, crs=newstr)
  shp
}

#--------------------------------------------
kernel2D<- function(parms, maxdim, np, eps) {
  # 2D kernel function fully vectorised
  calc.dist<- function(x,y){
    return(sqrt((x^2 + y^2))) 
  }  
  max.disp<- parms$sigma * 2.45
  mx<-  ceiling(max.disp/eps)
  while(length(-mx:mx) >= maxdim) {mx<- mx - 1}  #make sure final size of kernel is compatible
  gxy <- -mx:mx * eps
  dmat<- outer(gxy, gxy, FUN=calc.dist)
  if(!is.null(parms$ce)) 
    {lam<- exp(parms$b0 + parms$b1 * dmat + parms$ce)} #ce are other covariate effects
  else 
    {lam<- exp(parms$b1 + parms$b2 * dmat)}
  lam[dmat > max.disp]<- 0 #truncate beyond max.disp
  lam<- lam * np  # total risk is risk per period * number of periods (np)
  return(lam)
}


#--------------------------------------------
make.surface<- function(dat, parms, shape, nperiod, cellsize, Pu=1, prior, verbose=FALSE) {
  n<- nrow(dat)
  rast_map<- rast(vect(shape), resolution=cellsize)
  rast_map<- rasterize(vect(shape), rast_map, 0, background = 0)
  rast.mat<- as.matrix(rast_map, wide=TRUE)
  dims<- dim(rast.mat)
  for(i in 1:n){
    if(verbose) cat(paste("doing judas ",i," of ",n,sep=""),"\n")
    tmprast<- rast_map
    cells<- cellFromXY(tmprast, cbind(dat$xbar[i],dat$ybar[i]))
    tmprast[cells]<- 1
    wkern<- kernel2D(parms=parms[i,], maxdim=min(dims), np=nperiod[i], eps=cellsize)
    tmpmat<- as.matrix(tmprast, wide=TRUE)
    tmpmat<- simecol::neighbours(tmpmat, state=1, wdist=wkern)
    rast.mat<- rast.mat + tmpmat
  }
  values(rast_map)<- rast.mat
  den <- mask(rast_map, vect(shape))
  den<- app(den, function(x){1-exp(-x)}) #Probability scale
  
  N<- length(which(values(which.lyr(den >= 0))==1)) # Total cells
  nn<- length(which(values(which.lyr(den >= 0.0001))==1)) # Covered cells
  
  seu_avg<- global(den, 'mean', na.rm=TRUE)$mean
  SSe<- 1 - (1 - seu_avg * nn/N)^Pu
  PE <- prior/(1-SSe*(1-prior))
  result<- tibble(Cov=round(nn/N,3),SSe=round(SSe,3),PE=round(PE,3))
  list(Table=result, Raster=den)
}

#-----------------------------

rcov_effects<- function(beta, X) {
  X<- as.matrix(X)
  n<- nrow(X)
  nr<- nrow(as.matrix(beta))
  if(nr==1)
    rbetas<- as.matrix((replicate(n, rnorm(nr,beta$mean, beta$sd))))
  else
    rbetas<- t(replicate(n, rnorm(nr,beta$mean, beta$sd)))
  cov_eff<- rep(NA,n)
  for(i in 1:n) {
    cov_eff[i]<- X[i,] %*% rbetas[i,]
  }
  as.vector(cov_eff)
}

#-----------------------------

intersect.period<- function(x, y) {
  # find overlap between two date ranges x and y
  period<- case_when(
    (x[2] < y[1]) ~ 0,
    (y[2] < x[1]) ~ 0,
    ((x[1] <= y[1]) & (y[1] <= x[2]) & (x[2] <= y[2])) ~ as.numeric(x[2]-y[1]), 
    ((y[1] <= x[1]) & (x[1] <= y[2]) & (y[2] <= x[2])) ~ as.numeric(y[2]-x[1]),
    ((x[1] <= y[1]) & (y[1] <= y[2]) & (y[2] <= x[2])) ~ as.numeric(y[2]-y[1]),
    ((y[1] <= x[1]) & (x[1] <= x[2]) & (x[2] <= y[2])) ~ as.numeric(x[2]-x[1]))
  return(period)
}

#-----------------------------
summarise_sims<- function(.y, .x, period, ndec=3) {
  xbar<- mean(.x, na.rm=TRUE)
  xbar<- round(xbar, ndec)
  cl<- quantile(.x, c(0.05, 0.95))
  cl<- round(cl, ndec)
  tibble(Period=period, stat=.y, Mean=xbar,lcl=cl[1],ucl=cl[2])
}


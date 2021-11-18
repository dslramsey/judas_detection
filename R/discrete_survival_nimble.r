
#===================================
#
# Detection probability - survival analysis
# -----
library(tidyverse)
library(lubridate)
library(sf)
library(purrr)
library(nimble)
library(posterior)
source("R/judas functions.r")

goat_dyads<- read_csv("data/Santiago_judas_dyads.csv")
goat_dyads<- goat_dyads %>% mutate(id1 = factor(as.character(id1)), id2 = factor(as.character(id2)))
goat_ids<- goat_dyads %>% select(ID=id1,Sex=Sex1,Trans=Trans1) %>% distinct() %>% arrange(ID)

id1<- as.vector(unclass(goat_dyads$id1))
id2<- as.vector(unclass(goat_dyads$id2))

dyad<- as.vector(unclass(factor(goat_dyads$Dyad)))
ndyads<- max(dyad)
N<- nrow(goat_dyads)
M<- nrow(goat_ids)
X<- model.matrix(~ distance + Sex1 + Trans1, data=goat_dyads)[,-1]
K<- ncol(X)
sex<- model.matrix(~Sex + Trans, goat_ids)[,-1]

# Calculate home range scale (sigma) from Judas locations 

goat_locs<- read_csv("data/santiago_judas_locs.csv")
goat_locs<- goat_locs %>% mutate(ID = factor(as.character(ID)))
goat_locs<- goat_locs %>% group_by(ID) %>% nest() %>% arrange(ID)


goat_chr<- goat_locs %>% mutate(CNorm=map(data, calc_hr))
goat_chr<- goat_chr %>% unnest(CNorm) %>% select(-data)

maxD<- goat_chr$sigma * 2.45

#------------------------------------------------------
# Discrete survival analysis

# Nimble function to undertake estimation of Marginal detection probability (Eqn. 2)

integrateP<- function(b0, b1, beta, XX, maxD) {
  integrand<- function(x, b0, b1, beta, XX, maxD) {
    e<- 2*x/maxD^2 
    sex.effects<- as.vector(beta %*% XX)
    int<- e * (1-exp(-exp(b0 + b1*x + sex.effects)))
    return(int)
  }
  return(integrate(integrand,lower=0,upper=maxD, b0=b0, b1=b1, beta=beta, XX=XX, maxD=maxD,
                   rel.tol = .Machine$double.eps^0.5)[[1]])
}

RintegrateP <- nimbleRcall(function(b0=double(0), b1=double(0), beta=double(1),
                                    XX=double(1), maxD=double(0)){}, Rfun='integrateP',
                           returnType = double(0))

# Nimble code -  Note: since no time varying parameters/covariates we can group intervals
# and use a binomial distribution rather than bernoulli.

code<- nimbleCode({
  
  for (i in 1:N) {
    ev[i] ~ dbin(q[i], time[i])
    cloglog(q[i])<- b1 + inprod(X[i,1:K],beta[1:K]) + eps_id1[id1[i]] + eps_id2[id2[i]]
  }
  
  for(k in 1:M) {
    eps_id1[k] ~ dnorm(0, sd=sigma_id1)
    eps_id2[k] ~ dnorm(0, sd=sigma_id2)
    beta0[k] <- b1 + eps_id1[k]
    Pave[k]<- RintegrateP(beta0[k], beta[1], beta[2:K], sex[k,1:(K-1)], maxD[k])
  }
  
  for(j in 1:K) {
    beta[j] ~ dnorm(0, sd=5)
  }
  b1 ~ dnorm(0, sd=5)
  
  sigma_id1 ~ T(dt(0, 1, 4),0,)
  sigma_id2 ~ T(dt(0, 1, 4),0,)
 
})  	
#----------------------

constants<- list(N=N, M=M, K=K, id1=id1, id2=id2)
data<- list(time=goat_dyads$time,ev=goat_dyads$event, X=X, sex=sex, maxD=maxD)


inits1<- list(b1=-4, beta= c(-1,1,1,1), sigma_id1=runif(1),sigma_id2=runif(1),
              q=rep(0.01,N),eps_id1=runif(M,-2,0),eps_id2=runif(M, -2,0))

## create the model object
Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits1)
Rmcmc<- compileNimble(Rmodel)

ModSpec <- configureMCMC(Rmodel)
#ModSpec$removeSamplers(c("eps_id1"), print = FALSE)
#ModSpec$removeSamplers(c("eps_id2"), print = FALSE)
#ModSpec$addSampler(target =c("eps_id1"), type = 'AF_slice', print=FALSE)
#ModSpec$addSampler(target =c("eps_id2"), type = 'AF_slice', print=FALSE)
ModSpec$resetMonitors()
ModSpec$addMonitors(c("b1","beta","sigma_id1","sigma_id2","beta0","Pave"))
                      

Cmcmc <- buildMCMC(ModSpec)
Cmodel <- compileNimble(Cmcmc, project = Rmodel, resetFunctions = T)

n.iter<- 2000
n.burnin<- 1000
n.thin<- 1
n.chains<- 3

inits<- function(){inits1}

samp<- runMCMC(Cmodel, niter=n.iter, nburnin=n.burnin, thin=n.thin, nchains=n.chains,inits=inits,
               samplesAsCodaMCMC = TRUE)


tmp<- subset_draws(as_draws(samp), variable=c("b1","beta","sigma_id1","sigma_id2"))
summarise_draws(tmp)

tmp<- subset_draws(as_draws(samp), variable=c("Pave"))
summarise_draws(tmp)
#--------------------------------------------------------------------
# Plots
#-------------------------------------------------------------------- 

# Predicted probability of detection at activity centres (i.e. distance = 0)

fit<- as_draws_rvars(samp)

beta0<- draws_of(fit$beta0)
beta<- draws_of(fit$beta)
beta1<- beta[,1]
beta2<- beta[,-1]

n<- nrow(goat_ids)
XX<- model.matrix(~Sex + Trans, goat_ids)[,-1]

plam<- matrix(NA,nrow=dim(beta)[1],ncol=n)
for(i in 1:n){
  plam[,i]<- 1-exp(-exp(beta0[,i] + beta2 %*% XX[i,]))
}

plam<- as.data.frame(plam)
names(plam)<- goat_ids$ID
plam$samp<- 1:nrow(plam)

plam<- plam %>% pivot_longer(-samp, names_to = "ID", values_to = "p")
plam<- left_join(plam, goat_ids)
plam<- plam %>% mutate(Sex = factor(Sex, levels=c("Female","Super","Male")),
                       Trans = factor(Trans, labels = c("Not Translocated","Translocated")))

win.graph(10,5)
plam %>% ggplot(aes(ID, p, color=interaction(Sex,Trans))) +
  geom_boxplot(fill="grey90",outlier.color = NA) +
  labs(x="Judas ID", y="Detection probability at HR centre") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())



#--------------------------------------------------------
# Predicted detection curves with distance from HR centre
#--------------------------------------------------------
fit<- as_draws_rvars(samp)

beta0<- E(fit$beta0)
beta<- E(fit$beta)
beta1<- beta[1]
beta2<- beta[-1]

n<- nrow(goat_ids)
distance<- seq(0,20,0.1)
XX<- model.matrix(~Sex + Trans, goat_ids)[,-1]

plam<- matrix(NA,nrow=length(distance), ncol=n)

for(i in 1:n){
  plam[,i]<- 1-exp(-exp(beta0[i] + beta1*distance + as.vector(XX[i,] %*% beta2)))
}
plam<- as.data.frame(plam)
names(plam)<- goat_ids$ID
plam<- plam %>% mutate(Distance=distance)

plam<- plam %>% pivot_longer(-Distance, names_to = "ID", values_to = "p")
plam<- left_join(plam, goat_ids)
plam<- plam %>% mutate(Sex = factor(Sex, levels=c("Female","Super","Male")),
                       Trans = factor(Trans, labels = c("Not Translocated","Translocated")))


win.graph(8,7)
plam %>% ggplot(aes(Distance, p, group=ID)) +
  geom_line(color="grey50") +
  facet_grid(rows=vars(Sex), cols=vars(Trans)) +
  labs(x="Distance (km)",y="Probability of association") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(face="bold"))

#--------------------------------------------------------------------
# Unconditional detection probability (Eqn. 2)
#-------------------------------------------------------------------- 

fit<- subset_draws(as_draws_rvars(samp), variable="Pave")

plam<- goat_ids %>% mutate(pave= E(fit$Pave), lcl = quantile2(fit$Pave)[1,], ucl=quantile2(fit$Pave)[2,])
plam<- plam %>% mutate(Sex = factor(Sex, levels=c("Female","Super","Male")),
                       Trans = factor(Trans, labels = c("Not Translocated","Translocated")))

win.graph(10,5)
plam<- plam %>% arrange(pave)
plam %>%  mutate(ID = fct_reorder(ID, pave), Sex) %>%
  ggplot(aes(x=Trans, y=pave, color=Sex)) +
  geom_pointrange(aes(ymin=lcl,ymax=ucl), position=position_dodge2(width=0.5)) +
  facet_wrap(~Sex) +
  labs(x="",y="Probability of association",color="Status") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(face="bold"),
        axis.title = element_text(face="bold", size=12),
        axis.text.x = element_text(size=11))

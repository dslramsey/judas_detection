
#=========================================================
#
# Judas detection probability - discrete survival analysis
#
# Dave Ramsey (22/11/2021)
#=========================================================

library(tidyverse)
library(lubridate)
library(sf)
library(purrr)
library(cmdstanr)
library(posterior)
source("R/judas functions.r")

#==================================================
# Santiago Judas goats
#==================================================


goat_dyads<- read_csv("data/Santiago_judas_dyads.csv")
goat_dyads<- goat_dyads %>% mutate(id1 = factor(as.character(id1)), id2 = factor(as.character(id2)))
goat_ids<- read_csv("data/Santiago_judas_ID.csv")
goat_ids<- goat_ids %>% mutate(ID = factor(as.character(ID)))

id1<- as.vector(unclass(goat_dyads$id1))
id2<- as.vector(unclass(goat_dyads$id2))

dyad<- as.vector(unclass(factor(goat_dyads$Dyad)))
ndyads<- max(dyad)
N<- nrow(goat_dyads)
M<- nrow(goat_ids)
X<- model.matrix(~ distance + Sex1 + Trans1, data=goat_dyads)[,-1]
K<- ncol(X)
sex<- model.matrix(~Sex + Trans, goat_ids)[,-1]

#---- Calculate home range scale (sigma) from Judas locations ----

goat_locs<- read_csv("data/santiago_judas_locs.csv")
goat_locs<- goat_locs %>% mutate(ID = factor(as.character(ID)))
goat_locs<- goat_locs %>% group_by(ID) %>% nest() %>% arrange(ID)


goat_chr<- goat_locs %>% mutate(CNorm=map(data, calc_hr))
goat_chr<- goat_chr %>% unnest(CNorm) %>% select(-data)

maxD<- goat_chr$sigma * 2.45

#---- Discrete survival analysis -----

set.seed(123)

mod<- cmdstan_model("stan/discrete_survival.stan")

data<- list(N=N,M=M,K=K,X=X,id1=id1,id2=id2,
            time=as.integer(goat_dyads$time),
            ev=as.integer(goat_dyads$event),sex=sex,maxD=maxD)
            

## MCMC settings (for illustration - increase for serious inference)
ni <- 500
nt <- 1
nb <- 500
nc <- 2

 
inits <- lapply(1:nc, function(i)
  list(b1= -3, beta=c(-0.5,1,1,1), sigma_id1=runif(1), sigma_id2=runif(1), q=runif(N, 0, 0.1)))

out<- mod$sample(
  data=data,
  init = inits,
  iter_warmup = nb,
  iter_sampling = ni,
  chains=nc,
  parallel_chains = nc,
  refresh = 50)


out$summary(c("b1","beta","sigma_id1","sigma_id2"))

out$summary("pave")

out$save_object(file="out/goat_surv.rds")

#==================================================
# Santa Cruz Island Judas pigs
#==================================================

pig_dyads<- read_csv("data/Santacruz_judas_dyads.csv")
pig_dyads<- pig_dyads %>% mutate(id1 = factor(as.character(id1)), id2 = factor(as.character(id2)))
pig_ids<- pig_dyads %>% select(ID=id1,Sex=Sex1) %>% distinct() %>% arrange(ID)

id1<- as.vector(unclass(pig_dyads$id1))
id2<- as.vector(unclass(pig_dyads$id2))

dyad<- as.vector(unclass(factor(pig_dyads$Dyad)))
ndyads<- max(dyad)
N<- nrow(pig_dyads)
M<- nrow(pig_ids)
X<- model.matrix(~ distance + Sex1, data=pig_dyads)[,-1]
K<- ncol(X)
sex<- as.matrix(ifelse(pig_ids$Sex=="B", 0, 1))

#--- Calculate home range scale (sigma) from Judas locations
pig_locs<- read_csv("data/santacruz_judas_locs.csv")
# rescale coordinates to km
pig_locs<- pig_locs %>% mutate(ID = factor(as.character(ID)),East=East/1000, North=North/1000)
pig_locs<- pig_locs %>% group_by(ID) %>% nest() %>% arrange(ID)

pig_chr<- pig_locs %>% mutate(CNorm=map(data, calc_hr))
pig_chr<- pig_chr %>% unnest(CNorm) %>% select(-data)

maxD<- pig_chr$sigma * 2.45

#------------------------------------------------------
# Discrete survival analysis
#-----------------------------------------------------------

mod<- cmdstan_model("stan/discrete_survival.stan")

data<- list(N=N,M=M,K=K,X=X,id1=id1,id2=id2,
            time=as.integer(pig_dyads$time),
            ev=as.integer(pig_dyads$event),sex=sex,maxD=maxD)


## MCMC settings
ni <- 500
nt <- 1
nb <- 500
nc <- 2


inits <- lapply(1:nc, function(i)
  list(b1= -3, beta=c(-0.5,1), sigma_id1=runif(1), sigma_id2=runif(1), q=runif(N, 0, 0.1)))

out<- mod$sample(
  data=data,
  init = inits,
  iter_warmup = nb,
  iter_sampling = ni,
  chains=nc,
  parallel_chains = nc,
  refresh = 50)


out$summary(c("b1","beta","sigma_id1","sigma_id2"))

out$summary("pave")

out$save_object(file="out/pig_surv.rds")


#---------------------------------------------------
# Plots 
#---------------------------------------------------

# Predicted probability of detection at activity centres (i.e. distance = 0)

out<- readRDS("out/goat_surv.rds")

fit<- as_draws_rvars(out$draws(c("beta0","beta")))

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
# Predicted detection curves with distance from activity centre
#--------------------------------------------------------
fit<- as_draws_rvars(out$draws(c("beta0","beta")))

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

fit<- summarise_draws(out$draws("pave"), mean, ~quantile(.x, c(0.05,0.95)))
 
plam<- goat_ids %>% mutate(pave= fit$mean, lcl = fit$`5%`, ucl=fit$`95%`)
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


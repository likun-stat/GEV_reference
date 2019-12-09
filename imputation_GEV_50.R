source("~/Desktop/Research/GEV_reference/generic_samplers.R")
source("~/Desktop/Research/GEV_reference/ReferencePrior_utils.R")
source("~/Desktop/Research/GEV_reference/ReferencePrior_sampler.R")



######################################################################################
### -------------------------- Precipitation (xi=0.15) --------------------------------
######################################################################################

# library(extRemes)
xi_prior_type <- 'Beta'
mu <- 0
xi <- 0.15
tau <- 1
Time <- 50
mu_T_true <- qevd(1-1/Time, loc=0, shape=xi, scale=1)

## Generate fake data
n <- 50
n.experiments <- 100
Y_all <- matrix(NA, nrow=n.experiments, ncol=n)
set.seed(111)
for (iter in 1:n.experiments){
  Y_all[iter,] <- revd(n, loc=mu, scale=tau, shape=xi)
}



## Initial setup
Tau <- seq(0.1,10,length.out = 35)
if(xi>0.9) Mu_T <- seq(30,60,length.out = 35) else Mu_T <- seq(-10,10,length.out = 35)
Xi <- seq(-0.5,0.5,length.out = 35)
Par=expand.grid(tau=Tau, mu_T=Mu_T, xi=Xi)
Par<-cbind(mu=Par$mu, xi=Par$xi, tau=Par$tau)
lik_Par <- rep(NA,nrow(Par))

n.updates <- 10000
thin <- 10
echo.interval <- 1000


## -----------------------   Run through 100 experiments  -----------------------------
for (iter in 1:n.experiments){
  
  Y <- Y_all[iter,]
  
  ## Search initial values 
  for(i in 1:nrow(Par)){
    mu_T_tmp <- Par[i,1]; xi_tmp <- Par[i,2]; tau_tmp <- Par[i,3]
    tmp <- Lik3(Y=Y, mu_T=mu_T_tmp, xi=xi_tmp, tau=tau_tmp, Time=Time)
    if(xi_prior_type=="Beta") lik_Par[i]<-tmp+beta_prior(xi = xi_tmp)-log(tau_tmp)
    if(xi_prior_type=="h11") lik_Par[i]<-tmp+h11(xi_tmp)-log(tau_tmp)
    if(xi_prior_type=="h22_1") lik_Par[i]<--tmp+h22_1(xi_tmp)-log(tau_tmp)
    if(xi_prior_type=="h22_2") lik_Par[i]<--tmp+h22_2(xi_tmp)-log(tau_tmp)
    if(xi_prior_type=="MDI") lik_Par[i]<-tmp+MDI(xi = xi_tmp)-log(tau_tmp)
  }
  
  mu_T_init <- Par[which.max(lik_Par),1]
  xi_init <- Par[which.max(lik_Par),2]
  tau_init <- Par[which.max(lik_Par),3]
  
  initial.values <- list(mu_T=mu_T_init, xi=xi_init, tau=tau_init, Time=Time)
  true.params <- list(mu_T=mu_T_true, xi=xi, tau=tau)
  
  res <- GEV.sampler.02_T(Y=Y, iter.num = iter,
                        xi_prior_type=xi_prior_type,
                        initial.values=initial.values,
                        n.updates=n.updates, thin=thin,
                        experiment.name="50Beta-prior",
                        echo.interval=echo.interval,
                        true.params=true.params, sd.ratio=NULL, lower.prob.lim=0.5)
  
}





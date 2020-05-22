source("~/Desktop/Research/GEV_reference/generic_samplers.R")
source("~/Desktop/Research/GEV_reference/ReferencePrior_utils.R")
source("~/Desktop/Research/GEV_reference/ReferencePrior_sampler.R")

args <- commandArgs(trailingOnly=TRUE)

n <- as.numeric(args[1])
n.experiments <- as.numeric(args[2])


######################################################################################
### -------------------------- Precipitation (xi=0.15) --------------------------------
######################################################################################

# library(extRemes)
xi_prior_type <- 'Beta'
mu <- 0
xi <- 0.15
tau <- 1

## Generate fake data
# n <- 50
# n.experiments <- 100
Y_all <- matrix(NA, nrow=n.experiments, ncol=n)
set.seed(111)
for (iter in 1:n.experiments){
  Y_all[iter,] <- revd(n, loc=mu, scale=tau, shape=xi)
}



## Initial setup
Tau <- seq(0.1,10,length.out = 35)
Mu <- seq(-10,10,length.out = 35)
Xi <- seq(-0.5,0.5,length.out = 35)
Par=expand.grid(tau=Tau,mu=Mu,xi=Xi)
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
    tmp <- Lik(Par[i,],Y)
    xi_tmp <- Par[i,2];tau_tmp<-Par[i,3]
    if(xi_prior_type=="Beta") lik_Par[i]<-tmp+beta_prior(xi = xi_tmp)-log(tau_tmp)
    if(xi_prior_type=="h11") lik_Par[i]<-tmp+h11(xi_tmp)-log(tau_tmp)
    if(xi_prior_type=="h22_1") lik_Par[i]<--tmp+h22_1(xi_tmp)-log(tau_tmp)
    if(xi_prior_type=="h22_2") lik_Par[i]<--tmp+h22_2(xi_tmp)-log(tau_tmp)
    if(xi_prior_type=="MDI") lik_Par[i]<-tmp+MDI(xi = xi_tmp)-log(tau_tmp)
  }
  
  mu_init <- Par[which.max(lik_Par),1]
  xi_init <- Par[which.max(lik_Par),2]
  tau_init <- Par[which.max(lik_Par),3]
  
  initial.values <- list(mu=mu_init, xi=xi_init, tau=tau_init)
  true.params <- list(mu=mu, xi=xi, tau=tau)
  
  res <- GEV.sampler.02(Y=Y, iter.num = iter,
                        xi_prior_type=xi_prior_type,
                        initial.values=initial.values,
                        n.updates=n.updates, thin=thin,
                        experiment.name="Beta-prior",
                        echo.interval=echo.interval,
                        true.params=true.params, sd.ratio=NULL, lower.prob.lim=0.5)
  
}





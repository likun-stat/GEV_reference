
######################################################################################
### -------------------------- Precipitation (xi=0.15) --------------------------------
######################################################################################

library(extRemes)

mu <- 0
xi <- 0.15
tau <- 1

n <- 500
N <- 100
delta<-0.1

Tau <- seq(0.1,10,length.out = 35)
Mu <- seq(-10,10,length.out = 35)
Xi <- seq(-0.5,0.5,length.out = 35)
Par=expand.grid(tau=Tau,mu=Mu,xi=Xi)
Par<-cbind(mu=Par$mu, xi=Par$xi, tau=Par$tau)



lik_beta_Par <- rep(NA,nrow(Par))
lik_h11_Par <- rep(NA,nrow(Par))
lik_h22_1_Par <- rep(NA,nrow(Par))
lik_h22_2_Par <- rep(NA,nrow(Par))
lik_MDI_Par <- rep(NA,nrow(Par))

lik_beta_GML <- matrix(NA,nrow=N,ncol=3)
lik_h11_GML <- matrix(NA,nrow=N,ncol=3)
lik_h22_1_GML <- matrix(NA,nrow=N,ncol=3)
lik_h22_2_GML <- matrix(NA,nrow=N,ncol=3)
lik_MDI_GML <- matrix(NA,nrow=N,ncol=3)

set.seed(1234)
for(iter in 1:N){
  Y <- revd(n, loc=mu, scale=tau, shape=xi)
  
  
  for(i in 1:nrow(Par)){
    tmp <- Lik(Par[i,],Y)
    xi_tmp <- Par[i,2];tau_tmp<-Par[i,3]
    lik_beta_Par[i]<-tmp+beta_prior(xi = xi_tmp)
    lik_h11_Par[i]<-tmp+h11(xi_tmp)-log(tau_tmp)
    lik_h22_1_Par[i]<--tmp+h22_1(xi_tmp)-log(tau_tmp)
    lik_h22_2_Par[i]<--tmp+h22_2(xi_tmp)-log(tau_tmp)
    lik_MDI_Par[i]<-tmp+MDI(xi = xi_tmp)-log(tau_tmp)
  }
  
  
  ## For beta
  GML <- c(Par[which.max(lik_beta_Par),1],Par[which.max(lik_beta_Par),2], Par[which.max(lik_beta_Par),3])
  start <- c(GML[1],GML[2],GML[3])
  lower <- c(GML[1]-delta, GML[2]-delta, GML[3]-10*delta)
  upper <- c(GML[1]+delta, GML[2]+delta, GML[3]+10*delta)
  Res <- optim(start, Lik_beta, Y=Y, 
               # method = "L-BFGS-B", lower = lower, upper = upper, 
               control=list(trace=0,fnscale=-1))
  lik_beta_GML[iter,]<-Res$par
  
  
  ## For h11
  GML <- c(Par[which.max(lik_h11_Par),1],Par[which.max(lik_h11_Par),2], Par[which.max(lik_h11_Par),3])
  start <- c(GML[1],GML[2],GML[3])
  lower <- c(GML[1]-delta, GML[2]-delta, GML[3]-10*delta)
  upper <- c(GML[1]+delta, GML[2]+delta, GML[3]+10*delta)
  Res <- optim(start, Lik_h11, Y=Y, 
               # method = "L-BFGS-B", lower = lower, upper = upper, 
               control=list(trace=0,fnscale=-1))
  lik_h11_GML[iter,]<-Res$par

  ## For h22_1
  GML <- c(Par[which.max(lik_h22_1_Par),1],Par[which.max(lik_h22_1_Par),2], Par[which.max(lik_h22_1_Par),3])
  start <- c(GML[1],GML[2],GML[3])
  lower <- c(GML[1]-delta, GML[2]-delta, GML[3]-10*delta)
  upper <- c(GML[1]+delta, GML[2]+delta, GML[3]+10*delta)
  Res <- optim(start, Lik_h22_1, Y=Y, 
               # method = "L-BFGS-B", lower = lower, upper = upper, 
               control=list(trace=0,fnscale=-1))
  lik_h22_1_GML[iter,]<-Res$par

  
  ## For h22_2
  GML <- c(Par[which.max(lik_h22_2_Par),1],Par[which.max(lik_h22_2_Par),2], Par[which.max(lik_h22_2_Par),3])
  start <- c(GML[1],GML[2],GML[3])
  lower <- c(GML[1]-delta, GML[2]-delta, GML[3]-10*delta)
  upper <- c(GML[1]+delta, GML[2]+delta, GML[3]+10*delta)
  Res <- optim(start, Lik_h22_2, Y=Y, 
               # method = "L-BFGS-B", lower = lower, upper = upper, 
               control=list(trace=0,fnscale=-1))
  lik_h22_2_GML[iter,]<-Res$par
  
  
  ## For MDI
  GML <- c(Par[which.max(lik_MDI_Par),1],Par[which.max(lik_MDI_Par),2], Par[which.max(lik_MDI_Par),3])
  start <- c(GML[1],GML[2],GML[3])
  lower <- c(GML[1]-delta, GML[2]-delta, GML[3]-10*delta)
  upper <- c(GML[1]+delta, GML[2]+delta, GML[3]+10*delta)
  Res <- optim(start, Lik_MDI, Y=Y, 
               # method = "L-BFGS-B", lower = lower, upper = upper, 
               control=list(trace=0,fnscale=-1))
  lik_MDI_GML[iter,]<-Res$par
  
  cat('Done with',iter,'\n')
}


cat(round(mean(lik_beta_GML[,2] - xi),4),'&',
round(mean(lik_h11_GML[,2] - xi),4),'&',
round(mean(lik_h22_1_GML[,2] - xi),4),'&',
round(mean(lik_h22_2_GML[,2] - xi),4),'&',
round(mean(lik_MDI_GML[,2] - xi),4),'&','\n')

cat(round(sd(lik_beta_GML[,2]),4),'&',
round(sd(lik_h11_GML[,2]),4),'&',
round(sd(lik_h22_1_GML[,2]),4),'&',
round(sd(lik_h22_2_GML[,2]),4),'&',
round(sd(lik_MDI_GML[,2]),4),'&','\n')



lik_beta_50 <- rep(NA,N)
lik_h11_50 <- rep(NA,N)
lik_h22_1_50 <-rep(NA,N)
lik_h22_2_50 <- rep(NA,N)
lik_MDI_50 <- rep(NA,N)
lik_beta_100 <- rep(NA,N)
lik_h11_100 <- rep(NA,N)
lik_h22_1_100 <-rep(NA,N)
lik_h22_2_100 <- rep(NA,N)
lik_MDI_100 <- rep(NA,N)

for(iter in 1:N){
  lik_beta_50[iter] <- qevd(0.98, loc=lik_beta_GML[iter,1], shape=lik_beta_GML[iter,2], scale=lik_beta_GML[iter,3])
  lik_beta_100[iter] <- qevd(0.99, loc=lik_beta_GML[iter,1], shape=lik_beta_GML[iter,2], scale=lik_beta_GML[iter,3])
  
  lik_h11_50[iter] <- qevd(0.98, loc=lik_h11_GML[iter,1], shape=lik_h11_GML[iter,2], scale=lik_h11_GML[iter,3])
  lik_h11_100[iter] <- qevd(0.99, loc=lik_h11_GML[iter,1], shape=lik_h11_GML[iter,2], scale=lik_h11_GML[iter,3])
  
  lik_h22_1_50[iter] <- qevd(0.98, loc=lik_h22_1_GML[iter,1], shape=lik_h22_1_GML[iter,2], scale=lik_h22_1_GML[iter,3])
  lik_h22_1_100[iter] <- qevd(0.99, loc=lik_h22_1_GML[iter,1], shape=lik_h22_1_GML[iter,2], scale=lik_h22_1_GML[iter,3])
  
  lik_h22_2_50[iter] <- qevd(0.98, loc=lik_h22_2_GML[iter,1], shape=lik_h22_2_GML[iter,2], scale=lik_h22_2_GML[iter,3])
  lik_h22_2_100[iter] <- qevd(0.99, loc=lik_h22_2_GML[iter,1], shape=lik_h22_2_GML[iter,2], scale=lik_h22_2_GML[iter,3])
  
  lik_MDI_50[iter] <- qevd(0.98, loc=lik_MDI_GML[iter,1], shape=lik_MDI_GML[iter,2], scale=lik_MDI_GML[iter,3])
  lik_MDI_100[iter] <- qevd(0.99, loc=lik_MDI_GML[iter,1], shape=lik_MDI_GML[iter,2], scale=lik_MDI_GML[iter,3])
}

cat(round(mean(lik_beta_50 - qevd(0.98, loc=mu, shape=xi, scale=tau)),4),'&',
round(mean(lik_h11_50 - qevd(0.98, loc=mu, shape=xi, scale=tau)),4),'&',
round(mean(lik_h22_1_50 - qevd(0.98, loc=mu, shape=xi, scale=tau)),4),'&',
round(mean(lik_h22_2_50 - qevd(0.98, loc=mu, shape=xi, scale=tau)),4),'&',
round(mean(lik_MDI_50 - qevd(0.98, loc=mu, shape=xi, scale=tau)),4),'&','\n')

cat(round(sd(lik_beta_50),4),'&',
round(sd(lik_h11_50),4),'&',
round(sd(lik_h22_1_50),4),'&',
round(sd(lik_h22_2_50),4),'&',
round(sd(lik_MDI_50),4),'&','\n')


cat(round(mean(lik_beta_100 - qevd(0.99, loc=mu, shape=xi, scale=tau)),4),'&',
round(mean(lik_h11_100 - qevd(0.99, loc=mu, shape=xi, scale=tau)),4),'&',
round(mean(lik_h22_1_100 - qevd(0.99, loc=mu, shape=xi, scale=tau)),4),'&',
round(mean(lik_h22_2_100 - qevd(0.99, loc=mu, shape=xi, scale=tau)),4),'&',
round(mean(lik_MDI_100 - qevd(0.99, loc=mu, shape=xi, scale=tau)),4),'&','\n')

cat(round(sd(lik_beta_100),4),'&',
round(sd(lik_h11_100),4),'&',
round(sd(lik_h22_1_100),4),'&',
round(sd(lik_h22_2_100),4),'&',
round(sd(lik_MDI_100),4),'&','\n')




######################################################################################
### -------------------------- Temperature (xi= -0.3) --------------------------------
######################################################################################

library(extRemes)

mu <- 0
xi <- -0.3
tau <- 1

n <- 50
N <- 100
delta<-0.1

Tau <- seq(0.1,10,length.out = 35)
Mu <- seq(-10,10,length.out = 35)
Xi <- seq(-0.5,0.5,length.out = 35)
Par=expand.grid(tau=Tau,mu=Mu,xi=Xi)
Par<-cbind(mu=Par$mu, xi=Par$xi, tau=Par$tau)



lik_beta_Par <- rep(NA,nrow(Par))
lik_h11_Par <- rep(NA,nrow(Par))
lik_h22_1_Par <- rep(NA,nrow(Par))
lik_h22_2_Par <- rep(NA,nrow(Par))
lik_MDI_Par <- rep(NA,nrow(Par))

lik_beta_GML <- matrix(NA,nrow=N,ncol=3)
lik_h11_GML <- matrix(NA,nrow=N,ncol=3)
lik_h22_1_GML <- matrix(NA,nrow=N,ncol=3)
lik_h22_2_GML <- matrix(NA,nrow=N,ncol=3)
lik_MDI_GML <- matrix(NA,nrow=N,ncol=3)

set.seed(1234)
for(iter in 1:N){
  Y <- revd(n, loc=mu, scale=tau, shape=xi)
  
  
  for(i in 1:nrow(Par)){
    tmp <- Lik(Par[i,],Y)
    xi_tmp <- Par[i,2];tau_tmp<-Par[i,3]
    lik_beta_Par[i]<-tmp+beta_prior(xi = xi_tmp)
    lik_h11_Par[i]<-tmp+h11(xi_tmp)-log(tau_tmp)
    lik_h22_1_Par[i]<--tmp+h22_1(xi_tmp)-log(tau_tmp)
    lik_h22_2_Par[i]<--tmp+h22_2(xi_tmp)-log(tau_tmp)
    lik_MDI_Par[i]<-tmp+MDI(xi = xi_tmp)-log(tau_tmp)
  }
  
  
  ## For beta
  GML <- c(Par[which.max(lik_beta_Par),1],Par[which.max(lik_beta_Par),2], Par[which.max(lik_beta_Par),3])
  start <- c(GML[1],GML[2],GML[3])
  lower <- c(GML[1]-delta, GML[2]-delta, GML[3]-10*delta)
  upper <- c(GML[1]+delta, GML[2]+delta, GML[3]+10*delta)
  Res <- optim(start, Lik_beta, Y=Y, 
               # method = "L-BFGS-B", lower = lower, upper = upper, 
               control=list(trace=0,fnscale=-1))
  lik_beta_GML[iter,]<-Res$par
  
  
  ## For h11
  GML <- c(Par[which.max(lik_h11_Par),1],Par[which.max(lik_h11_Par),2], Par[which.max(lik_h11_Par),3])
  start <- c(GML[1],GML[2],GML[3])
  lower <- c(GML[1]-delta, GML[2]-delta, GML[3]-10*delta)
  upper <- c(GML[1]+delta, GML[2]+delta, GML[3]+10*delta)
  Res <- optim(start, Lik_h11, Y=Y, 
               # method = "L-BFGS-B", lower = lower, upper = upper, 
               control=list(trace=0,fnscale=-1))
  lik_h11_GML[iter,]<-Res$par
  
  ## For h22_1
  GML <- c(Par[which.max(lik_h22_1_Par),1],Par[which.max(lik_h22_1_Par),2], Par[which.max(lik_h22_1_Par),3])
  start <- c(GML[1],GML[2],GML[3])
  lower <- c(GML[1]-delta, GML[2]-delta, GML[3]-10*delta)
  upper <- c(GML[1]+delta, GML[2]+delta, GML[3]+10*delta)
  Res <- optim(start, Lik_h22_1, Y=Y, 
               # method = "L-BFGS-B", lower = lower, upper = upper, 
               control=list(trace=0,fnscale=-1))
  lik_h22_1_GML[iter,]<-Res$par
  
  
  ## For h22_2
  GML <- c(Par[which.max(lik_h22_2_Par),1],Par[which.max(lik_h22_2_Par),2], Par[which.max(lik_h22_2_Par),3])
  start <- c(GML[1],GML[2],GML[3])
  lower <- c(GML[1]-delta, GML[2]-delta, GML[3]-10*delta)
  upper <- c(GML[1]+delta, GML[2]+delta, GML[3]+10*delta)
  Res <- optim(start, Lik_h22_2, Y=Y, 
               # method = "L-BFGS-B", lower = lower, upper = upper, 
               control=list(trace=0,fnscale=-1))
  lik_h22_2_GML[iter,]<-Res$par
  
  
  ## For MDI
  GML <- c(Par[which.max(lik_MDI_Par),1],Par[which.max(lik_MDI_Par),2], Par[which.max(lik_MDI_Par),3])
  start <- c(GML[1],GML[2],GML[3])
  lower <- c(GML[1]-delta, GML[2]-delta, GML[3]-10*delta)
  upper <- c(GML[1]+delta, GML[2]+delta, GML[3]+10*delta)
  Res <- optim(start, Lik_MDI, Y=Y, 
               # method = "L-BFGS-B", lower = lower, upper = upper, 
               control=list(trace=0,fnscale=-1))
  lik_MDI_GML[iter,]<-Res$par
  
  cat('Done with',iter,'\n')
}


cat(round(mean(lik_beta_GML[,2] - xi),4),'&',
    round(mean(lik_h11_GML[,2] - xi),4),'&',
    round(mean(lik_h22_1_GML[,2] - xi),4),'&',
    round(mean(lik_h22_2_GML[,2] - xi),4),'&',
    round(mean(lik_MDI_GML[,2] - xi),4),'&','\n')

cat(round(sd(lik_beta_GML[,2]),4),'&',
    round(sd(lik_h11_GML[,2]),4),'&',
    round(sd(lik_h22_1_GML[,2]),4),'&',
    round(sd(lik_h22_2_GML[,2]),4),'&',
    round(sd(lik_MDI_GML[,2]),4),'&','\n')



lik_beta_50 <- rep(NA,N)
lik_h11_50 <- rep(NA,N)
lik_h22_1_50 <-rep(NA,N)
lik_h22_2_50 <- rep(NA,N)
lik_MDI_50 <- rep(NA,N)
lik_beta_100 <- rep(NA,N)
lik_h11_100 <- rep(NA,N)
lik_h22_1_100 <-rep(NA,N)
lik_h22_2_100 <- rep(NA,N)
lik_MDI_100 <- rep(NA,N)

for(iter in 1:N){
  lik_beta_50[iter] <- qevd(0.98, loc=lik_beta_GML[iter,1], shape=lik_beta_GML[iter,2], scale=lik_beta_GML[iter,3])
  lik_beta_100[iter] <- qevd(0.99, loc=lik_beta_GML[iter,1], shape=lik_beta_GML[iter,2], scale=lik_beta_GML[iter,3])
  
  lik_h11_50[iter] <- qevd(0.98, loc=lik_h11_GML[iter,1], shape=lik_h11_GML[iter,2], scale=lik_h11_GML[iter,3])
  lik_h11_100[iter] <- qevd(0.99, loc=lik_h11_GML[iter,1], shape=lik_h11_GML[iter,2], scale=lik_h11_GML[iter,3])
  
  lik_h22_1_50[iter] <- qevd(0.98, loc=lik_h22_1_GML[iter,1], shape=lik_h22_1_GML[iter,2], scale=lik_h22_1_GML[iter,3])
  lik_h22_1_100[iter] <- qevd(0.99, loc=lik_h22_1_GML[iter,1], shape=lik_h22_1_GML[iter,2], scale=lik_h22_1_GML[iter,3])
  
  lik_h22_2_50[iter] <- qevd(0.98, loc=lik_h22_2_GML[iter,1], shape=lik_h22_2_GML[iter,2], scale=lik_h22_2_GML[iter,3])
  lik_h22_2_100[iter] <- qevd(0.99, loc=lik_h22_2_GML[iter,1], shape=lik_h22_2_GML[iter,2], scale=lik_h22_2_GML[iter,3])
  
  lik_MDI_50[iter] <- qevd(0.98, loc=lik_MDI_GML[iter,1], shape=lik_MDI_GML[iter,2], scale=lik_MDI_GML[iter,3])
  lik_MDI_100[iter] <- qevd(0.99, loc=lik_MDI_GML[iter,1], shape=lik_MDI_GML[iter,2], scale=lik_MDI_GML[iter,3])
}

cat(round(mean(lik_beta_50 - qevd(0.98, loc=mu, shape=xi, scale=tau)),4),'&',
    round(mean(lik_h11_50 - qevd(0.98, loc=mu, shape=xi, scale=tau)),4),'&',
    round(mean(lik_h22_1_50 - qevd(0.98, loc=mu, shape=xi, scale=tau)),4),'&',
    round(mean(lik_h22_2_50 - qevd(0.98, loc=mu, shape=xi, scale=tau)),4),'&',
    round(mean(lik_MDI_50 - qevd(0.98, loc=mu, shape=xi, scale=tau)),4),'&','\n')

cat(round(sd(lik_beta_50),4),'&',
    round(sd(lik_h11_50),4),'&',
    round(sd(lik_h22_1_50),4),'&',
    round(sd(lik_h22_2_50),4),'&',
    round(sd(lik_MDI_50),4),'&','\n')


cat(round(mean(lik_beta_100 - qevd(0.99, loc=mu, shape=xi, scale=tau)),4),'&',
    round(mean(lik_h11_100 - qevd(0.99, loc=mu, shape=xi, scale=tau)),4),'&',
    round(mean(lik_h22_1_100 - qevd(0.99, loc=mu, shape=xi, scale=tau)),4),'&',
    round(mean(lik_h22_2_100 - qevd(0.99, loc=mu, shape=xi, scale=tau)),4),'&',
    round(mean(lik_MDI_100 - qevd(0.99, loc=mu, shape=xi, scale=tau)),4),'&','\n')

cat(round(sd(lik_beta_100),4),'&',
    round(sd(lik_h11_100),4),'&',
    round(sd(lik_h22_1_100),4),'&',
    round(sd(lik_h22_2_100),4),'&',
    round(sd(lik_MDI_100),4),'&','\n')


library(extRemes)
euler <- -digamma(1)
mu <- 20
tau <- 0.5
xi <- 0.2
n <- 100000
z <- revd(n, loc=mu, scale=tau, shape=xi)

W <- 1+ (xi+0.01)*(z-mu-0.01)/(tau+0.01)
sum(log(W))/n

W0 <- 1+ (xi)*(z-mu)/(tau)
sum(log(W0))/n
(xi+1)*euler
((xi+1)/xi)*sum(log(W0))/n


fu <- function(xi) 2^(2*xi)*gamma(xi+1/2)*(1+xi)/sqrt(pi)
Fu <- Vectorize(fu)
curve(Fu, from = -0.49, to = 1, ylim=c(0,10))
abline(h=1)

a=-2; delta=1
fu <- function(r) (1+r)^a-1-a*r^{a-1}-a*(a-1)*((1+r)^{a-2}+1)*r^2/2
Fu <- Vectorize(fu)
curve(Fu, from =-1, to = 10, ylim=c(-10,0))

a=2; delta=1
fu <- function(x) a^x- 1-x*log(a) - (a^{delta}+a^{delta})*(log(a))^2*x^2
Fu <- Vectorize(fu)
curve(Fu, from =-delta, to = delta, ylim=c(-10,0))


Lik<-function(par,Y){
  mu<-par[1]; xi<-par[2]; tau<-par[3]
  W <- 1+ (xi)*(Y-mu)/(tau)
  if(any(W<0)) return(-10^15) else return(-length(Y)*log(tau)-(xi+1)*sum(log(W))/xi-sum(W^{-1/xi}))
}

mu <- 20
tau <- 0.5
xi <- 0.2
n <- 1000000
Y <- revd(n, loc=mu, scale=tau, shape=xi)

Lik(c(tau/xi+min(Y)+0.1,xi,tau),Y)

max(Y)

delta<-0.02
mu+delta-(tau+delta)/(xi-delta)
Res <- optim(c(mu+delta,xi+delta,tau+delta), Lik, Y=Y, method = "L-BFGS-B", lower = c(mu-delta,xi-delta,tau-delta), upper = c(mu+delta,xi+delta,tau+delta), control=list(fnscale=-1))

Lik(Res$par,Y)
beta_hat<-Res$par[1]-Res$par[3]/Res$par[2]
beta_hat
mu-tau/xi
W_hat<-1+ (Res$par[2])*(Y-Res$par[1])/(Res$par[3])


tmp <- function(a){
  -n*sum((W_hat+a)^{-1/Res$par[2]-1})/sum((W_hat+a)^{-1/Res$par[2]})+(Res$par[2]+1)*sum((W_hat+a)^{-1})
}
tmp <- Vectorize(tmp)
curve(tmp, from=Res$par[2]*(beta_hat-min(Y)), to=0,xlab=expression(beta),ylab=expression(paste(L(theta)-L(hat(theta)))))
abline(v=Res$par[1]-Res$par[3]/Res$par[2])


## Check W0^{-1-1/xi} and W0^{-1-1/hat(xi)}
W0 <- 1+ (xi)*(Y-mu)/(tau)
W_hat<-1+ (Res$par[2])*(Y-Res$par[1])/(Res$par[3])
sum(W0^{-2-1/xi})/n
sum(W0^{-2-1/Res$par[2]})/n
sum(W_hat^{-2-1/Res$par[2]})/n
gamma(2*xi+2)

-xi*gamma(xi+2)*digamma(xi+2)
sum(W_hat^{-1-1/Res$par[2]}*log(W_hat))/n


sum(log(W_hat))/n





rdist<-function(a,b){sqrt(sum((a-b)^2))}

(Lik(par=c(18,10,1000),Y)-Lik(par=c(mu,tau,xi),Y))/(rdist(c(18,10,1000),c(mu,tau,xi)))^5

(Lik(par=c(18,tau,xi),Y)-Lik(par=c(mu,tau,xi),Y))/(rdist(c(18,tau,xi),c(mu,tau,xi)))^5

(Lik(par=c(min(Y),tau,xi),Y)-Lik(par=c(mu,tau,xi),Y))/(rdist(c(min(Y),tau,xi),c(mu,tau,xi)))^5





##########################################################################################
## ------------------------------- Generic Illustrations ---------------------------------
##########################################################################################

# calculate the likelihhood
Lik<-function(par,Y){
  mu<-par[1]; xi<-par[2]; tau<-par[3]
  if(tau<=0) return(NA)
  if(xi==0) return(-length(Y)*log(tau)-sum((Y-mu)/tau)-sum(exp(-(Y-mu)/tau)))
  W <- 1+ (xi)*(Y-mu)/(tau)
  if(any(W<0)) return(NA) else return(-length(Y)*log(tau)-(xi+1)*sum(log(W))/xi-sum(W^{-1/xi}))
}


bound <- function(tau,xi,alpha=0.1) tau/xi+min(Y)
bound_max <- function(tau,xi,alpha=0.1) tau/xi+max(Y)
beta_hat_line <- function(beta_hat,tau, xi) beta_hat + tau/xi
beta_bottom_line <- function(tau, xi) 16.4 + tau/xi

## ---------------- For positive shape parameter -------------------
set.seed(123)
xi_try<-0.2
mu <- 20
tau <- 0.5
n <- 1000
Y <- extRemes::revd(n, loc=mu, scale=tau, shape=xi_try)
delta<-0.02
mu+delta-(tau+delta)/(xi_try-delta)
Res <- optim(c(mu+delta,xi_try+delta,tau+delta), Lik, Y=Y, method = "L-BFGS-B", lower = c(mu-delta,xi_try-delta,tau-delta), upper = c(mu+delta,xi_try+delta,tau+delta), control=list(fnscale=-1))
beta_hat <- Res$par[1]-Res$par[3]/Res$par[2]



### 1. Generic support illustrations
pdf("~/Desktop/Illustration/Figures/supp_pos.pdf",width=5, height=5)
x<-cbind(tau=c(1000,0,0,1000,1000),mu=c(-1000,-1000,min(Y),bound(1000,0.1),-1000))
plot(x,type='n',xlab=expression(tau), ylab='',xlim=c(-0.3,1.5),ylim=c(0.3*min(Y),2.3*min(Y)),
     cex.main=1.6, cex.lab=1.5,
     yaxt='n', xaxt='n',main=expression(paste(xi,' > ',0)))
title(ylab=expression(mu),line=0, cex.lab=1.5)
axis(side=1, at=c(0,1.565), labels=c('0',expression(infinity)),cex.axis=1.5)
abline(v=0,col="grey",lty=2)
text(-0.12,min(Y)-0.05*min(Y),expression(Y[(1)]),pos=3,cex=1.5)

alpha=0.7
lines(x=c(0,-1000),y=c(min(Y),bound(-1000,0.1)),col="grey",lty=2)
polygon(x,col=scales::alpha('steelblue',alpha))
dev.off()


pdf("~/Desktop/Illustration/Figures/supp_zero.pdf",width=5, height=5)
x<-cbind(tau=c(1000,0,0,1000,1000),mu=c(-1000,-1000,1000,1000,-1000))
plot(x,type='n',xlab=expression(tau),ylab='',xlim=c(-0.3,1.5),ylim=c(0.3*min(Y),2.3*min(Y)),
     cex.main=1.6,cex.lab=1.5,
     yaxt='n', xaxt='n', main=expression(paste(xi==0)))
title(ylab=expression(mu),line=0, cex.lab=1.5)
axis(side=1, at=c(0,1.565), labels=c('0',expression(infinity)),cex.axis=1.5)
abline(v=0,col="grey",lty=2)

alpha=0.7
polygon(x,col=scales::alpha('steelblue',alpha))
dev.off()

pdf("~/Desktop/Illustration/Figures/supp_neg.pdf",width=5, height=5)
x<-cbind(tau=c(1000,0,0,1000,1000),mu=c(1000,1000,max(Y),bound_max(1000,-xi),1000))
plot(x,type='n',xlab=expression(tau),ylab='',xlim=c(-0.3,1.5),ylim=c(0.3*min(Y),2.3*min(Y)),
     cex.main=1.6,cex.lab=1.5,
     yaxt='n', xaxt='n',main=expression(paste(xi,' < ',0)))
title(ylab=expression(mu),line=0, cex.lab=1.5)
axis(side=1, at=c(0,1.565), labels=c('0',expression(infinity)),cex.axis=1.5)
abline(v=0,col="grey",lty=2)
text(-0.12,max(Y)-0.05*max(Y),expression(Y[(n)]),pos=3,cex=1.5)

alpha=0.7
lines(x=c(0,-1000),y=c(max(Y),bound_max(-1000,-xi)),col="grey",lty=2)
polygon(x,col=scales::alpha('steelblue',alpha))
dev.off()




### 2. Distance between Y(1) and the boundary
pdf("~/Desktop/Illustration/Figures/supp_pos_boundary.pdf",width=5, height=5)
x<-cbind(tau=c(1000,0,0,1000,1000),mu=c(-1000,-1000,min(Y),bound(1000,0.8*xi),-1000))
plot(x,type='n',xlab=expression(tau), ylab='',xlim=c(-0.3,1.5),ylim=c(0.3*min(Y),2.3*min(Y)),
     cex.main=1.6, cex.lab=1.5,
     yaxt='n', xaxt='n',main=expression(paste(xi[0],' > ',0)))
title(ylab=expression(mu),line=0, cex.lab=1.5)
axis(side=1, at=c(0,1.565), labels=c('0',expression(infinity)),cex.axis=1.5)
abline(v=0,col="grey",lty=2)
text(-0.12,min(Y)-0.05*min(Y),expression(Y[(1)]),pos=3,cex=1.5)

alpha=0.7
lines(x=c(0,-1000),y=c(min(Y),bound(-1000,0.8*xi)),col="grey",lty=2)
polygon(x,col=scales::alpha('steelblue',alpha))

lines(x=c(10,-1000),y=c(bound(10,0.8*xi)-5,bound(-1000,0.8*xi)-5),col="grey80",lty=2)
text(-0.08,min(Y)-5.1-0.05*min(Y),expression(beta[0]),pos=3,cex=1.5)
points(0.7,bound(0.7,0.8*xi)-5,pch=20,col='red', cex=1.8)
text(0.7,bound(0.7,0.8*xi)-5,expression(paste('(',tau[0],', ',mu[0],')')),pos=4,cex=1.5)
dev.off()


### 3. Convergence rate comparison
pdf("~/Desktop/Illustration/Figures/supp_rate_comp.pdf",width=5, height=5)
x<-cbind(tau=c(1000,0,0,1000,1000),mu=c(-1000,-1000,min(Y),bound(1000,0.8*xi),-1000))
plot(x,type='n',xlab=expression(tau), ylab='',xlim=c(-0.3,1.5),ylim=c(0.3*min(Y),2.3*min(Y)),
     cex.main=1.6, cex.lab=1.5,
     yaxt='n', xaxt='n',main=expression(paste(xi[0],' > ',0)))
title(ylab=expression(mu),line=0, cex.lab=1.5)
axis(side=1, at=c(0,1.565), labels=c('0',expression(infinity)),cex.axis=1.5)
abline(v=0,col="grey",lty=2)
text(-0.12,min(Y)-0.05*min(Y),expression(Y[(1)]),pos=3,cex=1.5)

alpha=0.7
lines(x=c(0,-1000),y=c(min(Y),bound(-1000,0.8*xi)),col="grey",lty=2)
polygon(x,col=alpha('steelblue',alpha))

lines(x=c(10,-1000),y=c(bound(10,0.8*xi)-5,bound(-1000,0.8*xi)-5),col="grey80",lty=2)
text(-0.08,min(Y)-5.1-0.05*min(Y),expression(beta[0]),pos=3,cex=1.5)
points(0.7,bound(0.7,0.8*xi)-5,pch=20,col='red', cex=1.2)
text(0.8,bound(0.7,0.8*xi)-5,expression(theta[0]),pos=1,cex=1.3)

Theta<-seq(0,2*pi,length.out=500)
circle<-cbind(0.7+0.05*cos(Theta),bound(0.7,0.8*xi)-5+0.05*sin(Theta)*23)
points(circle,type='l',lwd=2,col='red')
dev.off()


### 3(b). K_tilde
pdf("~/Desktop/Illustration/Figures/K_tilde.pdf",width=5, height=5)
x<-cbind(tau=c(1000,0,0,1000,1000),mu=c(-1000,-1000,min(Y),bound(1000,0.8*xi),-1000))
plot(x,type='n',xlab=expression(tau), ylab='',xlim=c(-0.3,1.5),ylim=c(0.3*min(Y),2.3*min(Y)),
     cex.main=1.6, cex.lab=1.5,
     yaxt='n', xaxt='n',main=expression(paste(xi[0],' > ',0)))
title(ylab=expression(mu),line=0, cex.lab=1.5)
axis(side=1, at=c(0,1.565), labels=c('0',expression(infinity)),cex.axis=1.5)
abline(v=0,col="grey",lty=2)
text(-0.12,min(Y)-0.05*min(Y),expression(Y[(1)]),pos=3,cex=1.5)

alpha=0.3
lines(x=c(0,-1000),y=c(min(Y),bound(-1000,0.8*xi)),col="grey",lty=2)
polygon(x,col=scales::alpha('steelblue',alpha),border="grey")

beta0<-mu-tau/xi_try
k<-cbind(tau=c(0.7-0.34,0.7+0.34,0.7+0.34,0.7-0.34,0.7-0.34),
         mu=c(beta_hat_line(beta0-13,0.7-0.34,xi = 0.8*xi),beta_hat_line(beta0-13,0.7+0.34,xi = 0.8*xi),
              beta_hat_line(beta0+6.7,0.7+0.34,xi = 0.8*xi),beta_hat_line(beta0+6.7,0.7-0.34,xi = 0.8*xi),
              beta_hat_line(beta0-13,0.7-0.34,xi = 0.8*xi)))
polygon(k,col=scales::alpha('yellow',alpha))

lines(x=c(10,-1000),y=c(bound(10,0.8*xi)-5,bound(-1000,0.8*xi)-5),col="grey80",lty=2)
lines(x=c(0.7,0.7),y=c(bound(0.7,0.8*xi)-5,-1000),col="grey80",lty=2)
lines(x=c(0.7-0.34,0.7-0.34),y=c(beta_hat_line(beta0-13,0.7-0.34,xi = 0.8*xi),-1000),col="grey80",lty=2)
text(-0.08,min(Y)-5.1-0.05*min(Y),expression(beta[0]),pos=3,cex=1.5)
points(0.7,bound(0.7,0.8*xi)-5,pch=20,col='red', cex=1.8)
text(0.7,bound(0.7,0.8*xi)-5,expression(paste('(',tau[0],', ',mu[0],')')),pos=4,cex=1.5)
dev.off()



### 4. Evaluate the log-likelihood at different cross section
pdf("~/Desktop/Illustration/Figures/sim_cross_section_pos.pdf",width=5, height=5)
m<-1.5  # Controls the cross section: xi_cross = m*xi_MLE
resolution <- 250
xi_MLE <- Res$par[2]
Tau=seq(0,1.7,length.out = resolution)
Mu=seq(16,24,length.out = resolution)
# Mu=seq(0.7*max(Y),1.4*max(Y),length.out = 200)
Par=expand.grid(tau=Tau,mu=Mu)
Par<-cbind(mu=Par$mu, xi=m*xi_MLE, tau=Par$tau)
lik_Par<-rep(NA,nrow(Par))
for(i in 1:nrow(Par)){
  tmp<-Lik(Par[i,],Y)
  lik_Par[i]<-tmp
}
range(lik_Par,na.rm = TRUE)
lik_Par[which(lik_Par<quantile(lik_Par,0.25,na.rm = TRUE))]<-quantile(lik_Par,0.25,na.rm = TRUE)

xlim=c(-0.3,1.5)
ylim=c(range(Mu)[1]+0.05*diff(range(Mu)),range(Mu)[2]-0.05*diff(range(Mu)))
x<-cbind(tau=c(0,0,1000,1000,0),mu=c(-1000,min(Y),bound(1000,m*xi_MLE),-1000,-1000))
# polygon(x=c(xlim[1],xlim[1],xlim[2],xlim[2],xlim[1]),y=c(ylim[1],ylim[2],ylim[2],ylim[1],ylim[1]),lty=2,col=alpha('steelblue',0.05))
plot(x,type='n',xlab='',ylab='',xlim=xlim, ylim=ylim,
     main=expression(paste(xi,' > ',hat(xi)[n])), cex.main=1.6)
title(xlab=expression(tau), ylab=expression(mu), cex.lab=1.5)
abline(v=0,col="grey",lty=2)
text(0.02,min(Y),expression(Y[(1)]),pos=2,cex=1.5)
image(x=Tau,y=Mu,
      z=matrix(lik_Par,length(Tau)),col=terrain.colors(40),add=TRUE)

alpha=0
lines(x=c(0,-1000),y=c(min(Y),bound(-1000,m*xi_MLE)),col="grey",lty=2)
polygon(x,col=scales::alpha('steelblue',alpha))
lines(x=c(0,-1000),y=c(beta_hat,beta_hat_line(beta_hat,-1000,xi=m*xi_MLE)),col="grey",lty=2)
lines(x=c(0,1000),y=c(beta_hat_line(beta_hat,0,xi=m*xi_MLE),beta_hat_line(beta_hat,1000,xi=m*xi_MLE)),lty=2)  # Beta hat line
text(0.02,beta_hat,expression(hat(beta)[n]),pos=2,cex=1.5)
points(Par[which.max(lik_Par),][3], Par[which.max(lik_Par),][1], col="blue", pch=20, cex=1.4)
dev.off()


pdf("~/Desktop/Illustration/Figures/sim_cross_section_zero.pdf",width=5, height=5)
m<-1  # Controls the cross section: xi_cross = m*xi_MLE
xi_MLE <- Res$par[2]
Tau=seq(0,1.7,length.out = resolution)
Mu=seq(16,24,length.out = resolution)
# Mu=seq(0.7*max(Y),1.4*max(Y),length.out = 200)
Par=expand.grid(tau=Tau,mu=Mu)
Par<-cbind(mu=Par$mu, xi=m*xi_MLE, tau=Par$tau)
lik_Par<-rep(NA,nrow(Par))
for(i in 1:nrow(Par)){
  tmp<-Lik(Par[i,],Y)
  lik_Par[i]<-tmp
}
range(lik_Par,na.rm = TRUE)
lik_Par[which(lik_Par<quantile(lik_Par,0.25,na.rm = TRUE))]<-quantile(lik_Par,0.25,na.rm = TRUE)

xlim=c(-0.3,1.5)
ylim=c(range(Mu)[1]+0.05*diff(range(Mu)),range(Mu)[2]-0.05*diff(range(Mu)))
x<-cbind(tau=c(0,0,1000,1000,0),mu=c(-1000,min(Y),bound(1000,m*xi_MLE),-1000,-1000))
# polygon(x=c(xlim[1],xlim[1],xlim[2],xlim[2],xlim[1]),y=c(ylim[1],ylim[2],ylim[2],ylim[1],ylim[1]),lty=2,col=alpha('steelblue',0.05))
plot(x,type='n',xlab='',ylab='',xlim=xlim, ylim=ylim,
     main=expression(paste(xi==hat(xi)[n])), cex.main=1.6)
title(xlab=expression(tau), ylab=expression(mu), cex.lab=1.5)
abline(v=0,col="grey",lty=2)
text(0.02,min(Y),expression(Y[(1)]),pos=2,cex=1.5)
image(x=Tau,y=Mu,
      z=matrix(lik_Par,length(Tau)),col=terrain.colors(40),add=TRUE)

alpha=0
lines(x=c(0,-1000),y=c(min(Y),bound(-1000,m*xi_MLE)),col="grey",lty=2)
polygon(x,col=scales::alpha('steelblue',alpha))
lines(x=c(0,-1000),y=c(beta_hat,beta_hat_line(beta_hat,-1000,xi=m*xi_MLE)),col="grey",lty=2)
lines(x=c(0,1000),y=c(beta_hat_line(beta_hat,0,xi=m*xi_MLE),beta_hat_line(beta_hat,1000,xi=m*xi_MLE)),lty=2)  # Beta hat line
text(0.02,beta_hat,expression(hat(beta)[n]),pos=2,cex=1.5)
points(Par[which.max(lik_Par),][3], Par[which.max(lik_Par),][1], col="blue", pch=20, cex=1.4)
dev.off()


pdf("~/Desktop/Illustration/Figures/sim_cross_section_neg.pdf",width=5, height=5)
m<-0.5  # Controls the cross section: xi_cross = m*xi_MLE
xi_MLE <- Res$par[2]
Tau=seq(0,1.7,length.out = resolution)
Mu=seq(16,24,length.out = resolution)
# Mu=seq(0.7*max(Y),1.4*max(Y),length.out = 200)
Par=expand.grid(tau=Tau,mu=Mu)
Par<-cbind(mu=Par$mu, xi=m*xi_MLE, tau=Par$tau)
lik_Par<-rep(NA,nrow(Par))
for(i in 1:nrow(Par)){
  tmp<-Lik(Par[i,],Y)
  lik_Par[i]<-tmp
}
range(lik_Par,na.rm = TRUE)
lik_Par[which(lik_Par<quantile(lik_Par,0.25,na.rm = TRUE))]<-quantile(lik_Par,0.25,na.rm = TRUE)

xlim=c(-0.3,1.5)
ylim=c(range(Mu)[1]+0.05*diff(range(Mu)),range(Mu)[2]-0.05*diff(range(Mu)))
x<-cbind(tau=c(0,0,1000,1000,0),mu=c(-1000,min(Y),bound(1000,m*xi_MLE),-1000,-1000))
# polygon(x=c(xlim[1],xlim[1],xlim[2],xlim[2],xlim[1]),y=c(ylim[1],ylim[2],ylim[2],ylim[1],ylim[1]),lty=2,col=alpha('steelblue',0.05))
plot(x,type='n',xlab='',ylab='',xlim=xlim, ylim=ylim,
     main=expression(paste(xi,' < ',hat(xi)[n])), cex.main=1.6)
title(xlab=expression(tau), ylab=expression(mu), cex.lab=1.5)
abline(v=0,col="grey",lty=2)
text(0.02,min(Y),expression(Y[(1)]),pos=2,cex=1.5)
image(x=Tau,y=Mu,
      z=matrix(lik_Par,length(Tau)),col=terrain.colors(40),add=TRUE)

alpha=0
lines(x=c(0,-1000),y=c(min(Y),bound(-1000,m*xi_MLE)),col="grey",lty=2)
polygon(x,col=scales::alpha('steelblue',alpha))
lines(x=c(0,-1000),y=c(beta_hat,beta_hat_line(beta_hat,-1000,xi=m*xi_MLE)),col="grey",lty=2)
lines(x=c(0,1000),y=c(beta_hat_line(beta_hat,0,xi=m*xi_MLE),beta_hat_line(beta_hat,1000,xi=m*xi_MLE)),lty=2)  # Beta hat line
text(0.02,beta_hat,expression(hat(beta)[n]),pos=2,cex=1.5)
points(Par[which.max(lik_Par),][3], Par[which.max(lik_Par),][1], col="blue", pch=20, cex=1.4)
dev.off()

lik_Par[which.max(lik_Par)]
Res$value




pdf("~/Desktop/Illustration/Figures/sim_cross_section_shape_zero.pdf",width=5, height=5)
m<-0  # Controls the cross section: xi_cross = m*xi_MLE
xi_MLE <- Res$par[2]
Tau=seq(0,1.7,length.out = resolution)
Mu=seq(16,24,length.out = resolution)
# Mu=seq(0.7*max(Y),1.4*max(Y),length.out = 200)
Par=expand.grid(tau=Tau,mu=Mu)
Par<-cbind(mu=Par$mu, xi=m*xi_MLE, tau=Par$tau)
lik_Par<-rep(NA,nrow(Par))
for(i in 1:nrow(Par)){
  tmp<-Lik(Par[i,],Y)
  lik_Par[i]<-tmp
}
range(lik_Par,na.rm = TRUE)
lik_Par[which(lik_Par<quantile(lik_Par,0.25,na.rm = TRUE))]<-quantile(lik_Par,0.25,na.rm = TRUE)

xlim=c(-0.3,1.5)
ylim=c(range(Mu)[1]+0.05*diff(range(Mu)),range(Mu)[2]-0.05*diff(range(Mu)))
x<-cbind(tau=c(0,0,1000,1000,0),mu=c(-1000,min(Y),bound(1000,m*xi_MLE),-1000,-1000))
# polygon(x=c(xlim[1],xlim[1],xlim[2],xlim[2],xlim[1]),y=c(ylim[1],ylim[2],ylim[2],ylim[1],ylim[1]),lty=2,col=alpha('steelblue',0.05))
plot(x,type='n',xlab='',ylab='',xlim=xlim, ylim=ylim,
     main=expression(paste(xi,' = ',0)), cex.main=1.6)
title(xlab=expression(tau), ylab=expression(mu), cex.lab=1.5)
abline(v=0,col="grey",lty=2)
text(0.02,min(Y),expression(Y[(1)]),pos=2,cex=1.5)
image(x=Tau,y=Mu,
      z=matrix(lik_Par,length(Tau)),col=terrain.colors(40),add=TRUE)

alpha=0
lines(x=c(0,-1000),y=c(min(Y),bound(-1000,m*xi_MLE)),col="grey",lty=2)
polygon(x,col=scales::alpha('steelblue',alpha))
lines(x=c(0,-1000),y=c(beta_hat,beta_hat_line(beta_hat,-1000,xi=m*xi_MLE)),col="grey",lty=2)
lines(x=c(0,1000),y=c(beta_hat_line(beta_hat,0,xi=m*xi_MLE),beta_hat_line(beta_hat,1000,xi=m*xi_MLE)),lty=2)  # Beta hat line
text(0.02,beta_hat,expression(hat(beta)[n]),pos=2,cex=1.5)
points(Par[which.max(lik_Par),][3], Par[which.max(lik_Par),][1], col="blue", pch=20, cex=1.4)
dev.off()

# Compared to global MLE
lik_Par[which.max(lik_Par)]
Res$value



## ---------------- For negative shape parameter -------------------
mu <- 20
tau <- 0.5
xi <- -0.2
n <- 1000
Y <- revd(n, loc=mu, scale=tau, shape=xi)

library(scales)
bound <- function(tau,xi) tau/xi+max(Y)

pdf("~/Desktop/Illustration/illus1.pdf",width=5, height=5)
plot(x,type='n',xlab=expression(tau),ylab='',xlim=c(-0.1,0.7),ylim=c(0.9*max(Y),1.03*max(Y)),yaxt='n',main=expression(xi==hat(xi)[n]))
title(ylab=expression(mu),line=0, cex.lab=1)
abline(v=0,col="grey",lty=2)
text(-0.055,max(Y)-0.011*max(Y),"max(Y)",pos=3)
text(-0.02,mu-tau/xi-0.008*max(Y),expression(hat(beta)[n]),pos=3)

alpha=0.1
x<-cbind(tau=c(0,0,1000,1000,0),mu=c(1000,max(Y),bound(1000,xi),1000,1000))
lines(x=c(0,-1000),y=c(max(Y),bound(-1000,xi)),col="grey",lty=2)
polygon(x,col=alpha('steelblue',alpha))

x<-cbind(tau=c(0,0,1000,1000,0),mu=c(1000,mu-tau/xi,bound(1000,xi)-max(Y)+mu-tau/xi,1000,1000))
lines(x=c(0,-1000),y=c(mu-tau/xi,bound(-1000,xi)-max(Y)+mu-tau/xi),col="grey",lty=2)
polygon(x,col=alpha('steelblue',alpha+0.3))
text(0.75,min(Y)-0.3*min(Y),expression(beta[0]),pos=3)

points(tau,mu,pch=20)
r=0.04
Theta<-seq(0,2*pi,length.out=500)
circle<-cbind(tau+r*cos(Theta),mu+r*sin(Theta)/0.24)
points(circle,type='l',lwd=2,col='red')
dev.off()




##########################################################################################
## ---------------------------------- Generic beta(xi) -----------------------------------
##########################################################################################

# calculate the likelihhood
Lik<-function(par,Y){
  mu<-par[1]; xi<-par[2]; tau<-par[3]
  if(tau<=0) return(NA)
  W <- 1+ (xi)*(Y-mu)/(tau)
  if(any(W<0)) return(NA) else return(-length(Y)*log(tau)-(xi+1)*sum(log(W))/xi-sum(W^{-1/xi}))
}


bound <- function(tau,xi,alpha=0.1) tau/xi+min(Y)
bound_max <- function(tau,xi,alpha=0.1) tau/xi+max(Y)
beta_hat_line <- function(beta_hat,tau, xi) beta_hat + tau/xi
beta_bottom_line <- function(tau, xi) 16.4 + tau/xi

set.seed(123)
xi_try<-0.2
mu <- 20
tau <- 0.5
n <- 1000
Y <- extRemes::revd(n, loc=mu, scale=tau, shape=xi_try)
delta<-0.02
mu+delta-(tau+delta)/(xi_try-delta)
Res <- optim(c(mu+delta,xi_try+delta,tau+delta), Lik, Y=Y, method = "L-BFGS-B", lower = c(mu-delta,xi_try-delta,tau-delta), upper = c(mu+delta,xi_try+delta,tau+delta), control=list(fnscale=-1))
beta_hat <- Res$par[1]-Res$par[3]/Res$par[2]

## Use simulations to see what the graph looks like
xi_cross<- seq(-0.999999,0.9,length=50) # Controls the cross section: xi_cross = m*xi_MLE
Tau=seq(0,12,length.out = 200)
Mu=seq(16,50,length.out = 250)
max_Lik<-rep(NA,length(xi_cross))
max_mu<-rep(NA,length(xi_cross))
max_tau<-rep(NA,length(xi_cross))
for(j in 1:length(xi_cross)){
  Par=expand.grid(tau=Tau,mu=Mu)
  Par<-cbind(mu=Par$mu, xi=xi_cross[j], tau=Par$tau)
  lik_Par<-rep(NA,nrow(Par))
  for(i in 1:nrow(Par)){
    tmp<-Lik(Par[i,],Y)
    lik_Par[i]<-tmp
  }
  max_Lik[j]<-max(lik_Par,na.rm=TRUE)
  max_mu[j]<-Par[which.max(lik_Par),1]
  max_tau[j]<-Par[which.max(lik_Par),3]
  cat('j=',j,'\n')
}

plot(xi_cross,max_mu-max_tau/xi_cross, type='l',ylim=c(-20,40))
fun1<- function(x) max(Y)-0.23/x
curve(fun,from=-1,to=-0.001,add=TRUE,col='red')
fun2<- function(x) min(Y)-0.5/x
curve(fun,from=0.001,to=1,add=TRUE,col='blue')

fun1<- function(x) max(Y)+5-0.23/x
fun2<- function(x) min(Y)-5-0.5/x
pdf("~/Desktop/Illustration/beta_xi.pdf",width=5, height=5)
curve(fun2,from=0.0126, to=1.8, n=500, lwd=1.3, xlab=expression(xi),ylab='',xlim=c(-1,1.8), ylim=c(-25,65),xaxt='n',yaxt='n')
curve(fun1,from=-1, to=-0.007, n=500, lwd=1.3, add=TRUE)
title(ylab=expression(beta[n](xi)),line=0.5, cex.lab=1.5)
abline(v=0,col="grey65",lty=2)
abline(v=-1,col="grey65",lty=3)
abline(v=1.8,col="grey65",lty=3)
points(x=c(-1,0),y=c(max(Y)+5,max(Y)+5),type='l',lty=2,col='grey50')
points(x=c(0,1.8),y=c(min(Y)-5,min(Y)-5),type='l',lty=2,col='grey50')
text(-0.05, max(Y)+5-0.011*max(Y),expression(Y[(n)]),pos=4,cex=1.5)
text(0.05,min(Y)-5-0.011*max(Y),expression(Y[(1)]),pos=2,cex=1.5)
# text(-0.05, 59, expression(infinity),pos=4)
# text(0.05,-19,expression(-infinity),pos=2)
axis(side=1, at=c(-1,0,1.8), labels=c('-1','0','n-1'))
dev.off()


##########################################################################################
## ---------------------------------- 3d support plot-------------------------------------
##########################################################################################

hessian <- function(par, Y, Res){
  mu<-par[1]; xi<-par[2]; tau<-par[3]
  beta <- mu-tau/xi
  n<-length(Y)
  mu_hat<-Res$par[1]; xi_hat<-Res$par[2];tau_hat<-Res$par[3]
  beta_hat <- mu_hat-tau_hat/xi_hat
  W_hat <- 1+ (xi_hat)*(Y-mu_hat)/(tau_hat)
  pTth<-W_hat+xi_hat*(beta_hat-beta)/tau_hat
  sum2<-sum(pTth^{-2})
  sum2xi<-sum(pTth^{-2-1/xi_hat})
  sum1xi<-sum(pTth^{-1-1/xi_hat})
  sumxi<-sum(pTth^{-1/xi_hat})
  # -n*(1+1/xi_hat)*sum2+n*(tau_hat/tau)^{-1/xi_hat}*(1/xi_hat)*(1/xi_hat+1)*sum2xi-
  #   (tau_hat/tau)^{-1/xi_hat}*(1/xi_hat^2-1)*sumxi*sum2+
  #   (tau_hat/tau)^{-2/xi_hat}*(1/xi_hat)*(1/xi_hat^2-1)*sumxi*sum2xi-
  #   (tau_hat/tau)^{-2/xi_hat}*(1/xi_hat^3)*(sum1xi)^2
  xi_hat*sum2-sum2xi
}


## 3-dimensional plot
beta0<-mu-tau/xi
tau=seq(0,1.5,length.out =100)
mu=seq(beta0+0.1,beta0+10,length.out =100)
volcano<-expand.grid(list(tau, mu))
xi<-t(matrix(volcano[,1]/(volcano[,2]-beta0),100,100))

library(plotly)
a <- list(text = "Support", xref = "paper", yref = "paper", 
          yanchor = "bottom", xanchor = "center", align = "center", 
          x = 0.5, y = 0.92, showarrow = FALSE)

scene1 <- list(aspectmode = "manual", aspectratio = list(x = 1, y = 1, z = 0.95),
               zaxis = list(showline = TRUE), 
               camera = list(eye = list(x = 0.9, y = -2.7, z = 0.3)))

# volcano is a numeric matrix that ships with R
p <- plot_ly(z = ~xi,x=~tau,y=~mu) %>% add_surface(colorscale = list(c(0,'#BA52ED'), c(1,'#FCB040')), showscale = FALSE) %>% 
  layout(annotations = a, scene = scene1)
p

Sys.setenv("plotly_username"="lfz5044")
Sys.setenv("plotly_api_key"="nn8hSl9HkAYX9RMvHw3L")
chart_link = api_create(p, filename="surface-1")
chart_link






##########################################################################################
## -------------------------------- Plot the h functions ---------------------------------
##########################################################################################

## h11
zeta=1.20205
gam=-digamma(1)
11*pi^6/2160-zeta^2
h11<-11*pi^4/360-6*zeta^2/pi^2
sqrt(h11)

h<- function(xi){
  if(xi+1/2<1e-4) return(sqrt(2*pi^2/3))
  p<- (1+xi)^2*gamma(2*xi+1)
  q<- gamma(2+xi)*(digamma(1+xi)+1/xi+1)
  s<-1-0.5772157+1/xi
  # r<-q-p/xi
  tmp <- (pi^2/6)-(q-s*gamma(xi+2))^2/(p-gamma(xi+2)^2)
  h11<-sqrt(tmp)/abs(xi)
  return(h11)
}
h<-Vectorize(h)
curve(h, from=-0.49,to=20,xlim=c(-0.5,20))
points(0,sqrt(h11))

curve(h, from=-0.5,to=30,xlab=expression(xi), ylab=expression(h[11]^{1/2}))
tmp<-function(xi) if(xi<0.45) return(NA) else sqrt(pi^2/(6*xi^2))
tmp<-Vectorize(tmp)
curve(tmp,from=-0.5,to=30,add=TRUE,col='red',lty=2)

library(ggplot2)
dat<-data.frame(xi=c(-1/2,0),y=c(sqrt(2*pi^2/3),sqrt(h11)))
# dat1<-data.frame(xi=c(0,0),y=c(-1,sqrt(h11)))
# dat2<-data.frame(xi=c(-1,0),y=c(sqrt(h11),sqrt(h11)))
x<-c(seq(-0.5,-0.05,by=0.01),seq(0.05,10,by=0.01))
h_approx<-approxfun(x=x,y=h(x))
H11<-h_approx

ggplot(data.frame(xi=c(-0.5, 10)), aes(xi)) + 
  stat_function(fun=h_approx, col='#036180',alpha=0.5, size=1.1) +
  stat_function(fun=tmp, linetype="dashed", col='#0a0001',alpha=0.9) +
  # scale_colour_manual("Function", values=c("blue","red"), breaks=c("square","exp")) +
  # geom_line(data=dat1,linetype="dashed", aes(y=y)) +
  # geom_line(data=dat2,linetype="dashed", aes(y=y)) +
  geom_point(data=dat, size=2, colour = "black", fill = "white", aes(y=y)) +
  xlab(expression(xi))+ylab(expression(sqrt(h[11])))+
  coord_cartesian(xlim=c(-0.51, 10), ylim=c(0.1,2.6)) +
  scale_y_continuous(breaks =c(0,0.5,1,sqrt(h11),1.5,2,2.5,sqrt(2*pi^2/3)),labels = c('0.0','0.5','1.0','','1.5','2.0','2.5','')) +
  scale_x_continuous(breaks =c(-0.5,0,2.5,5,7.5,10),labels = c(-0.5,'  0',2.5,5,7.5,10)) +
  theme(panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major.x = element_line(colour = "grey95"),
        panel.grid.major.y = element_line(colour = "grey95"),
        legend.position = 'none', legend.justification = c("right","top"), legend.title=element_blank(), 
        legend.text=element_text(size=10), legend.spacing.x = unit(0.2,'cm'),
        legend.box.background = element_rect(colour = "black",size=0.7),
        legend.text.align=0.5,
        plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"))
ggsave("~/Desktop/Illustration/h11_1.pdf",width = 4, height = 4)
# ggsave("~/Desktop/Illustration/h11_1.eps",width = 4, height = 4, units="in", device="eps")


ggplot(data.frame(xi=c(10, 20)), aes(xi)) + 
  stat_function(fun=function(xi) sqrt(pi^2/(6*xi^2)), linetype="dashed", col='#0a0001') +
  stat_function(fun=h, alpha=0.3, size=1.1, col = '#036180') +
  scale_colour_manual("Function", values=c("blue","red"), breaks=c("square","exp")) +
  xlab(expression(xi))+ylab(expression(sqrt(h[11])))+
  coord_cartesian(ylim=c(0.1,2.6)) +
  scale_y_continuous(breaks =c(0,0.5,1,1.5,2,2.5)) +
  theme(panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major.x = element_line(colour = "grey95"),
        panel.grid.major.y = element_line(colour = "grey95"),
        legend.position = c(.27, .94), legend.justification = c("right","top"), legend.title=element_blank(), 
        legend.text=element_text(size=10), legend.spacing.x = unit(0.2,'cm'),
        legend.box.background = element_rect(colour = "black",size=0.7),
        legend.text.align=0.5,
        plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"))
ggsave("~/Desktop/Illustration/h11_2.pdf",width = 4, height = 4)


## h33
zeta=1.20205
gam=-digamma(1)
dev2<-gam^2+pi^2/6
dev3<--gam^3-3*gam*pi^2/6-2*zeta
dev4<-gam^4+6*gam^2*pi^2/6+8*gam*zeta+27*pi^4/(90*2)
h33<- dev2+dev3+dev4/4
sqrt(h33)

h<- function(xi){
  if(xi+1/2<1e-4) return(sqrt(2*pi^2/3))
  p<- (1+xi)^2*gamma(2*xi+1)
  q<- gamma(2+xi)*(digamma(1+xi)+1/xi+1)
  s<-1-0.5772157+1/xi
  # r<-q-p/xi
  tmp <- (pi^2/6)+s^2-2*q/xi+p/xi^2
  h33<-sqrt(tmp)/abs(xi)
  return(h33)
}
h<-Vectorize(h)
curve(h, from=-0.49,to=1,ylim=c(0,10))
points(0,sqrt(h33))

curve(h, from=-0.5,to=30,xlab=expression(xi), ylab=expression(h[11]^{1/2}))
tmp<-function(xi) sqrt(pi^2/(6*xi^2))
curve(tmp,from=-0.5,to=30,add=TRUE,col='red',lty=2)

library(ggplot2)
dat<-data.frame(xi=c(0),y=c(sqrt(h33)))
# dat1<-data.frame(xi=c(0,0),y=c(-1,sqrt(h11)))
# dat2<-data.frame(xi=c(-1,0),y=c(sqrt(h11),sqrt(h11)))
x<-c(seq(-0.499,-0.05,by=0.001),seq(0.05,10,by=0.01))
h_approx<-approxfun(x=x,y=h(x))
H33<-h_approx

fun1<-function(xi) {if(xi<1) return(-1);sqrt((gamma(2*xi+1)-2*gamma(xi+1)*(1+digamma(xi+1))+pi^2/6+(1-gam)^2)/xi^2)}
fun1<-Vectorize(fun1)

ggplot(data.frame(xi=c(-0.49, 3.15)), aes(xi)) + 
  stat_function(fun=h_approx, col='#d113be',alpha=0.7, size=1.1) +
  stat_function(fun=function(xi) sqrt(4/(2*xi+1)), linetype="dashed", col='#0a0001',alpha=0.9) +
  # stat_function(fun=fun1, linetype="dashed", col='#0a0001',alpha=0.9) +
  geom_point(data=dat, size=2, colour = "black", fill = "white", aes(y=y)) +
  xlab(expression(xi))+ylab(expression(sqrt(h[33])))+
  coord_cartesian(xlim=c(-0.49, 3.15), ylim=c(0.5,15)) +
  scale_y_continuous(breaks=c(0,sqrt(h33),3,6,9,12,15),labels = c(0,'',3,6,9,12,15)) +
  scale_x_continuous(breaks =c(-0.52,0,1,2,3),labels = c(-0.5,0,1,2,3)) +
  theme(panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major.x = element_line(colour = "grey95"),
        panel.grid.major.y = element_line(colour = "grey95"),
        legend.position = 'none', legend.justification = c("right","top"), legend.title=element_blank(), 
        legend.text=element_text(size=10), legend.spacing.x = unit(0.2,'cm'),
        legend.box.background = element_rect(colour = "black",size=0.7),
        legend.text.align=0.5,
        plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"))
ggsave("~/Desktop/Illustration/h33_1.pdf",width = 4, height = 4)
# ggsave("~/Desktop/Illustration/h11_1.eps",width = 4, height = 4, units="in", device="eps")



## h22_1
zeta=1.20205
gam=-digamma(1)
dev2<-gam^2+pi^2/6
dev3<--gam^3-3*gam*pi^2/6-2*zeta
dev4<-gam^4+6*gam^2*pi^2/6+8*gam*zeta+27*pi^4/(90*2)
h22<- dev2+dev3+dev4/4-(-gam+3*dev2/2+dev3/2)^2/((1-gam)^2+pi^2/6)
sqrt(h22)

h22_1<-2*pi^2/3+4*(1-gam)^2
sqrt(h22_1)

h<- function(xi){
  if(xi+1/2<1e-4) return(sqrt(2*pi^2/3))
  p<- (1+xi)^2*gamma(2*xi+1)
  q<- gamma(2+xi)*(digamma(1+xi)+1/xi+1)
  s<-1-0.5772157+1/xi
  # r<-q-p/xi
  tmp <- (pi^2/6)+(1-0.5772157)^2-(gamma(2+xi)/xi-q+1-0.5772157)^2/(1+p-2*gamma(xi+2))
  h22<-sqrt(tmp)/abs(xi)
  return(h22)
}

h<-Vectorize(h)
curve(h, from=-0.49,to=10,ylim=c(0,4))
points(0,sqrt(h22))
points(-0.5,sqrt(h22_1))
tmp<-function(xi) {6/(xi+exp(1))}
tmp<-Vectorize(tmp)
curve(tmp,from=0,to=30,add=TRUE,col='red',lty=2)

curve(h, from=-0.5,to=30,xlab=expression(xi), ylab=expression(h[11]^{1/2}))
tmp<-function(xi) {if(xi<0.24) return(NA) else sqrt(pi^2/6+(1-gam)^2)/xi}
tmp<-Vectorize(tmp)
curve(tmp,from=0,to=30,add=TRUE,col='red',lty=2)

library(ggplot2)
dat<-data.frame(xi=c(-0.5,0),y=c(sqrt(h22_1),sqrt(h22)))
# dat1<-data.frame(xi=c(0,0),y=c(-1,sqrt(h11)))
# dat2<-data.frame(xi=c(-1,0),y=c(sqrt(h11),sqrt(h11)))
x<-c(seq(-0.4999,-0.05,by=0.001),seq(0.05,10,by=0.01))
h_approx<-approxfun(x=c(-0.5,x),y=c(sqrt(h22_1),h(x)))
H22_1<-h_approx


ggplot(data.frame(xi=c(-0.5, 10)), aes(xi)) + 
  stat_function(fun=h_approx, col='#deaf04',alpha=0.7, size=1.1) +
  stat_function(fun=tmp, linetype="dashed", col='#0a0001',alpha=0.9) +
  # stat_function(fun=fun1, linetype="dashed", col='#0a0001',alpha=0.9) +
  geom_point(data=dat, size=2, colour = "black", fill = "white", aes(y=y)) +
  xlab(expression(xi))+ylab(expression(sqrt(h[22])))+
  coord_cartesian(xlim=c(-0.5, 10), ylim=c(0.1,4.3)) +
  scale_y_continuous(breaks=c(0,1,sqrt(h22),2,sqrt(h22_1),3,4),labels = c(0,1,'',2,'',3,4)) +
  scale_x_continuous(breaks =c(-0.5,0,2.5,5,7.5,10),labels = c(-0.5,'  0',2.5,5,7.5,10)) +
  theme(panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major.x = element_line(colour = "grey95"),
        panel.grid.major.y = element_line(colour = "grey95"),
        legend.position = 'none', legend.justification = c("right","top"), legend.title=element_blank(), 
        legend.text=element_text(size=10), legend.spacing.x = unit(0.2,'cm'),
        legend.box.background = element_rect(colour = "black",size=0.7),
        legend.text.align=0.5,
        plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"))
ggsave("~/Desktop/Illustration/h22_1.pdf",width = 4, height = 4)
# ggsave("~/Desktop/Illustration/h11_1.eps",width = 4, height = 4, units="in", device="eps")



## h22_2
zeta=1.20205
gam=-digamma(1)
dev2<-gam^2+pi^2/6
dev3<--gam^3-3*gam*pi^2/6-2*zeta
dev4<-gam^4+6*gam^2*pi^2/6+8*gam*zeta+27*pi^4/(90*2)
h22<- pi^2*(1-gam)^2/6+2*(gam-1)*zeta+11*pi^4/360
sqrt(h22)

h22_1<-2*pi^2/3+4*(1+gam)^2
sqrt(h22_1)

h<- function(xi){
  if(xi+1/2<1e-4) return(sqrt(2*pi^2/3))
  p<- (1+xi)^2*gamma(2*xi+1)
  q<- gamma(2+xi)*(digamma(1+xi)+1/xi+1)
  s<-1-0.5772157+1/xi
  # r<-q-p/xi
  tmp <- (pi^2/6)+s^2-q^2/p
  h22<-sqrt(tmp)/abs(xi)
  return(h22)
}

h<-Vectorize(h)
curve(h, from=-0.499,to=10,ylim=c(0,4))
points(0,sqrt(h22))
points(-0.5,sqrt(h22_1))
tmp<-function(xi) {6/(xi+exp(1))}
tmp<-Vectorize(tmp)
curve(tmp,from=0,to=30,add=TRUE,col='red',lty=2)

curve(h, from=-0.5,to=30,xlab=expression(xi), ylab=expression(h[11]^{1/2}))
tmp<-function(xi) {if(xi<0.24) return(NA) else sqrt(pi^2/6+(1-gam)^2)/xi}
tmp<-Vectorize(tmp)
curve(tmp,from=0,to=30,add=TRUE,col='red',lty=2)

library(ggplot2)
dat<-data.frame(xi=c(-0.5,0),y=c(sqrt(h22_1),sqrt(h22)))
x<-c(seq(-0.4999,-0.05,by=0.001),seq(0.05,10,by=0.01))
h_approx<-approxfun(x=c(-0.5,x),y=c(sqrt(h22_1),h(x)))
H22_2<-h_approx


ggplot(data.frame(xi=c(-0.5, 10)), aes(xi)) + 
  stat_function(fun=h_approx, col='#2c20b0',alpha=0.7, size=1.1) +
  stat_function(fun=tmp, linetype="dashed", col='#0a0001',alpha=0.9) +
  # stat_function(fun=fun1, linetype="dashed", col='#0a0001',alpha=0.9) +
  geom_point(data=dat, size=2, colour = "black", fill = "white", aes(y=y)) +
  xlab(expression(xi))+ylab(expression(sqrt(h[22])))+
  coord_cartesian(xlim=c(-0.5, 10), ylim=c(0.1,4.3)) +
  scale_y_continuous(breaks=c(0,1,sqrt(h22),2,sqrt(h22_1),3,4),labels = c(0,1,'',2,'',3,4)) +
  scale_x_continuous(breaks =c(-0.5,0,2.5,5,7.5,10),labels = c(-0.5,'  0',2.5,5,7.5,10)) +
  theme(panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major.x = element_line(colour = "grey95"),
        panel.grid.major.y = element_line(colour = "grey95"),
        legend.position = 'none', legend.justification = c("right","top"), legend.title=element_blank(), 
        legend.text=element_text(size=10), legend.spacing.x = unit(0.2,'cm'),
        legend.box.background = element_rect(colour = "black",size=0.7),
        legend.text.align=0.5,
        plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"))
ggsave("~/Desktop/Illustration/h22_2.pdf",width = 4, height = 4)
# ggsave("~/Desktop/Illustration/h11_1.eps",width = 4, height = 4, units="in", device="eps")







##########################################################################################
## ------------------------ Compare h function to other functions ------------------------
##########################################################################################

## h22 compare
tmp<-function(xi) {if(xi<0) return(NA) else 6/(xi+exp(1))}
tmp<-Vectorize(tmp)
ggplot(data.frame(xi=c(-0.5, 10)), aes(xi)) + 
  stat_function(fun=H22_1, alpha=0.5, size=1.1, aes(color = "line1")) +
  stat_function(fun=H22_2 ,alpha=0.5, size=1.1, aes(color = "line2")) +
  stat_function(fun=tmp, linetype="dashed", col='#0a0001',alpha=0.9) +
  xlab(expression(xi))+ylab(expression(sqrt(h[22])))+
  coord_cartesian(xlim=c(-0.5, 10), ylim=c(0.1,4.5)) +
  scale_colour_manual("Lgend title", values = c("line1"='#deaf04',"line2"='#2c20b0'),
                      labels=c(expression(paste("(",mu,', ',xi,', ', tau,')')),expression(paste("(",tau,', ',xi,', ', mu,')')))) +
  theme(panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major.x = element_line(colour = "grey95"),
        panel.grid.major.y = element_line(colour = "grey95"),
        legend.position =c(0.95,0.95), legend.justification = c("right","top"), legend.title=element_blank(), 
        legend.text=element_text(size=15), legend.spacing.x = unit(0.2,'cm'),
        legend.box.background = element_rect(colour = "black",size=0.7),
        legend.text.align=0.5,
        plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"))
ggsave("~/Desktop/Illustration/h22_comp.pdf",width = 4, height = 4)


## h11 compare
tmp<-function(xi) {if(xi<0) return(NA) else 6/(xi+exp(1))}
tmp<-Vectorize(tmp)
ggplot(data.frame(xi=c(-0.5, 10)), aes(xi)) + 
  stat_function(fun=H11, alpha=0.7, size=1.1, aes(color = "line1")) +
  stat_function(fun=tmp, linetype="dashed", col='#0a0001',alpha=0.9) +
  xlab(expression(xi))+ylab(expression(sqrt(h[11])))+
  coord_cartesian(xlim=c(-0.5, 10), ylim=c(0.1,4.5)) +
  scale_colour_manual("Lgend title", values = c("line1"='#036180'),labels=expression(paste("(",xi,', ',tau,', ', mu,')'))) +
  theme(panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major.x = element_line(colour = "grey95"),
        panel.grid.major.y = element_line(colour = "grey95"),
        legend.position =c(0.95,0.95), legend.justification = c("right","top"), legend.title=element_blank(), 
        legend.text=element_text(size=15), legend.spacing.x = unit(0.2,'cm'),
        legend.box.background = element_rect(colour = "black",size=0.7),
        legend.text.align=0.5,
        plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"))
ggsave("~/Desktop/Illustration/h11_comp.pdf",width = 4, height = 4)


##h33 compare
tmp <- function(xi) if(xi<0.29) return(NA) else sqrt(2^{2/xi-1.5}*xi^{xi-2}*gamma(xi+1))
tmp<-Vectorize(tmp)
ggplot(data.frame(xi=c(-0.499, 3.5)), aes(xi)) + 
  geom_vline(xintercept=3,   colour="grey86") +
  geom_vline(xintercept=-0.52,colour="grey86") +
  stat_function(fun=H33, alpha=0.7, size=1.1, aes(color = "line1")) +
  stat_function(fun=tmp, linetype="dashed", col='#0a0001',alpha=0.9) +
  # stat_function(fun=fun1, linetype="dashed", col='#0a0001',alpha=0.9) +
  xlab(expression(xi))+ylab(expression(sqrt(h[33])))+
  coord_cartesian(xlim=c(-0.499, 3.5), ylim=c(0.5,15)) +
  scale_colour_manual("Lgend title", values = c("line1"='#d113be'),labels=expression(paste("(",mu,', ',tau,', ', xi,')'))) +
  # scale_y_continuous(breaks=c(0,sqrt(h33),3,6,9,12,15),labels = c(0,'',3,6,9,12,15)) +
  scale_x_continuous(breaks =c(-0.52,0,1,2,3),labels = c(-0.5,0,1,2,3)) +
  theme(panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major.x = element_line(colour = "grey95"),
        panel.grid.major.y = element_line(colour = "grey95"),
        legend.position = c(0.95,0.95), legend.justification = c("right","top"), legend.title=element_blank(), 
        legend.text=element_text(size=15), legend.spacing.x = unit(0.2,'cm'),
        legend.box.background = element_rect(colour = "black",size=0.7),
        legend.text.align=0.5,
        plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"))
ggsave("~/Desktop/Illustration/h33_comp.pdf",width = 4, height = 4)


## Jeffreys: sqrt(I(theta))
h<- function(xi){
  p<- (1+xi)^2*gamma(2*xi+1)
  q<- gamma(2+xi)*(digamma(1+xi)+1/xi+1)
  s<-1-0.5772157+1/xi
  # r<-q-p/xi
  tmp <- pi^2/6*(p-gamma(xi+2)^2)-(q-s*gamma(xi+2))^2
  h11<-sqrt(tmp)/xi^2
  return(h11)
}
x=c(-0.04,0.04);y<-c(h(x[1]),h(x[2]))
linear<-approxfun(x=x, y=y)
h<- function(xi){
  if(xi<x[2] & xi>x[1]) return(linear(xi))
  p<- (1+xi)^2*gamma(2*xi+1)
  q<- gamma(2+xi)*(digamma(1+xi)+1/xi+1)
  s<-1-0.5772157+1/xi
  # r<-q-p/xi
  tmp <- pi^2/6*(p-gamma(xi+2)^2)-(q-s*gamma(xi+2))^2
  h11<-sqrt(tmp)/xi^2
  return(h11)
}

Jeffrey<-Vectorize(h)
curve(Jeffrey, from=-0.49, to =0.4)
curve(tmp, add=TRUE)

tmp <- function(xi) if(xi<0.29) return(NA) else sqrt(2^{2/xi-1.5}*xi^{xi-2}*gamma(xi+1))
tmp<-Vectorize(tmp)
ggplot(data.frame(xi=c(-0.499, 3.5)), aes(xi)) + 
  geom_vline(xintercept=3,   colour="grey86") +
  geom_vline(xintercept=-0.52,colour="grey86") +
  stat_function(fun=Jeffrey, alpha=0.7, size=1.1, aes(color = "line1")) +
  stat_function(fun=tmp, linetype="dashed", col='#5c5657',alpha=0.9) +
  # stat_function(fun=function(xi) sqrt(4/(2*xi+1)), linetype="dashed", col='#0a0001',alpha=0.9) +
  xlab(expression(xi))+ylab(expression(tau^2*sqrt(I(theta))))+
  coord_cartesian(xlim=c(-0.499, 3.5), ylim=c(0.5,15)) +
  scale_colour_manual("Lgend title", values = c("line1"='#5c5657'),labels="Jeffreys prior") +
  # scale_y_continuous(breaks=c(0,sqrt(h33),3,6,9,12,15),labels = c(0,'',3,6,9,12,15)) +
  scale_x_continuous(breaks =c(-0.52,0,1,2,3),labels = c(-0.5,0,1,2,3)) +
  theme(panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major.x = element_line(colour = "grey95"),
        panel.grid.major.y = element_line(colour = "grey95"),
        legend.position = c(0.95,0.95), legend.justification = c("right","top"), legend.title=element_blank(), 
        legend.text=element_text(size=15), legend.spacing.x = unit(0.2,'cm'),
        legend.box.background = element_rect(colour = "black",size=0.7),
        legend.text.align=0.5,
        plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"))
ggsave("~/Desktop/Illustration/Jeffreys_comp.pdf",width = 4, height = 4)




##########################################################################################
## ---------------------------------- Analysis plots -------------------------------------
##########################################################################################

m<-1.1; fun<-function(x) gamma(2*x+1+1/m)-(1+m*x)*gamma(2*x+1)*gamma(1+1/m)
curve(fun, from=0.01,to=2, ylim=c(-8,8)); abline(h=0,lty=2,col='red')

m<-1.2; fun<-function(x) gamma(2*x+1+1/m)-(1+m*x)*gamma(2*x+1)*gamma(1+1/m)
curve(fun, from=0.01,to=2, add=TRUE)

m<-1.3; fun<-function(x) gamma(2*x+1+1/m)-(1+m*x)*gamma(2*x+1)*gamma(1+1/m)
curve(fun, from=0.01,to=2, add=TRUE)

m<-1.4; fun<-function(x) gamma(2*x+1+1/m)-(1+m*x)*gamma(2*x+1)*gamma(1+1/m)
curve(fun, from=0.01,to=2, add=TRUE)

m<-1.5; fun<-function(x) gamma(2*x+1+1/m)-(1+m*x)*gamma(2*x+1)*gamma(1+1/m)
curve(fun, from=0.01,to=2, add=TRUE)

m<-1.6; fun<-function(x) gamma(2*x+1+1/m)-(1+m*x)*gamma(2*x+1)*gamma(1+1/m)
curve(fun, from=0.01,to=2, add=TRUE)






m <- 1.1; fun <- function(x) (m*x+1)*gamma(x+1)-gamma(x+1/m+1)/gamma(1/m+1)
m <- 2; fun <- function(x) sqrt(pi)*(2*x+1)*gamma(x+1)/2-(x+1/2)*gamma(x+1/2) # get rid of denominators
m <- 2; fun <- function(x) sqrt(pi)*gamma(x+1)/gamma(x+1/2)-1 # get rid of denominators

curve(fun, from=0.01,to=0.2)
curve(fun, from=0.01,to=2)



fun <- function(m) -log(m)-log(gamma(1+m))+(1-1/m)*gam
curve(fun, from=1,to=2)

first_dev_gamma<-function(x)  digamma(x)*gamma(x)
curve(first_dev_gamma,from=2.5,to=4)

second_dev_gamma<-function(x)  (trigamma(x) + digamma(x)^2)*gamma(x)
curve(second_dev_gamma,from=1,to=4)


fun <- function(xi) {xi0=0.5;-xi0*gamma(xi0/xi+xi0+1)*first_dev_gamma(xi0/xi+1)+
  xi0*gamma(xi0/xi+1)*first_dev_gamma(xi0/xi+xi0+1)+xi^2*gamma(xi0/xi+1)^2*gamma(xi0+1)}
curve(fun,from=0.2, to=0.5)

fun <- function(xi) {xi0=2.4;-xi0*gamma(xi0/xi+xi0+1)*first_dev_gamma(xi0/xi+1)+
  xi0*gamma(xi0/xi+1)*first_dev_gamma(xi0/xi+xi0+1)+xi^2*gamma(xi0/xi+1)^2*gamma(xi0+1)}
curve(fun,from=0.8, to=2.4)

fun <- function(xi) {xi0=0.5;gamma(xi0/xi+1)*(gamma(xi0/xi+2*xi0+1)-gamma(xi0+1)*gamma(xi0/xi+xi0+1)-
  xi*gamma(xi0/xi+1)*gamma(2*xi0+1))}
curve(fun,from=0.2, to=0.5)
abline(h=0)

fun <- function(xi) {xi0=2.4;gamma(xi0/xi+1)*(gamma(xi0/xi+2*xi0+1)-gamma(xi0+1)*gamma(xi0/xi+xi0+1)-
                                                xi*gamma(xi0/xi+1)*gamma(2*xi0+1))}
curve(fun,from=0.8, to=2.4)

fun <- function(xi) {xi0=0.5;(2+1/xi)*gamma(xi0/xi+3*xi0+1)}
curve(fun,from=0.2, to=0.5)

fun <- function(xi) {xi0=2.4;gamma(xi0/xi+2*xi0+1)-gamma(xi0+1)*gamma(xi0/xi+xi0+1)-
  xi*gamma(xi0/xi+1)*gamma(2*xi0+1)}
curve(fun,from=0.8, to=2.4)

fun1 <- function(xi0) {-xi0*gamma(xi0+4)*first_dev_gamma(4)+
  xi0*gamma(4)*first_dev_gamma(xi0+4)+(xi0/3)^2*gamma(4)^2*gamma(xi0+1)}

fun2 <- function(xi0) {gamma(4)*(gamma(2*xi0+4)-gamma(xi0+1)*gamma(xi0+4)-
                                   (xi0/3)*gamma(4)*gamma(2*xi0+1))}
fun<-function(xi0) fun1(xi0)/fun2(xi0)
curve(fun, from=0.2, to=10)


##########################################################################################
## ---------------------------------- Derivative on xi -----------------------------------
##########################################################################################
# calculate the likelihhood
Lik<-function(par,Y){
  mu<-par[1]; xi<-par[2]; tau<-par[3]
  if(tau<=0) return(NA)
  if(xi==0) return(-length(Y)*log(tau)-sum((Y-mu)/tau)-sum(exp(-(Y-mu)/tau)))
  W <- 1+ (xi)*(Y-mu)/(tau)
  if(any(W<0)) return(NA) else return(-length(Y)*log(tau)-(xi+1)*sum(log(W))/xi-sum(W^{-1/xi}))
}


bound <- function(tau,xi,alpha=0.1) tau/xi+min(Y)
bound_max <- function(tau,xi,alpha=0.1) tau/xi+max(Y)
beta_hat_line <- function(beta_hat,tau, xi) beta_hat + tau/xi
beta_bottom_line <- function(tau, xi) 16.4 + tau/xi

# generate the samples
set.seed(123)
xi_try<-0.2
mu <- 20
tau <- 0.5
n <- 1000
Y <- extRemes::revd(n, loc=mu, scale=tau, shape=xi_try)
delta<-0.02
mu+delta-(tau+delta)/(xi_try-delta)
Res <- optim(c(mu+delta,xi_try+delta,tau+delta), Lik, Y=Y, method = "L-BFGS-B", lower = c(mu-delta,xi_try-delta,tau-delta), upper = c(mu+delta,xi_try+delta,tau+delta), control=list(fnscale=-1))
beta_hat <- Res$par[1]-Res$par[3]/Res$par[2]


M<-seq(-0.4,2,length.out = 20)  # Controls the cross section: xi_cross = m*xi_MLE
optim_Par_cross_section <- matrix(NA, nrow=length(M),ncol=3)
optim_Lik_cross_section <- rep(NA,length(M))
for(j in 1:length(M)){
  m<-M[j]
  xi_MLE <- Res$par[2]
  Tau=seq(0,1.7,length.out = 400)
  Mu=seq(16,24,length.out = 400)
  # Mu=seq(0.7*max(Y),1.4*max(Y),length.out = 200)
  Par=expand.grid(tau=Tau,mu=Mu)
  Par<-cbind(mu=Par$mu, xi=m*xi_MLE, tau=Par$tau)
  lik_Par<-rep(NA,nrow(Par))
  for(i in 1:nrow(Par)){
    tmp<-Lik(Par[i,],Y)
    lik_Par[i]<-tmp
  }
  cat("j=",j,'\n')
  optim_Par_cross_section[j,]<-Par[which.max(lik_Par),]
  optim_Lik_cross_section[j]<-max(lik_Par,na.rm=TRUE)
}
# save(optim_Lik_cross_section,file='~/Desktop/temp.RData')

pdf("~/Desktop/Illustration/PL_pos.pdf",width=5, height=5)
plot(optim_Par_cross_section[,2],optim_Lik_cross_section,type='l',col='red',lwd=2,
     xlab=expression(xi),ylab=expression(PL[n](xi)),main=expression(xi[0]==0.2))
abline(v=0.2,lty=2)
grid()
dev.off()



xi_cross<- seq(-0.999999,0.9,length=300) # Controls the cross section: xi_cross = m*xi_MLE
Tau=seq(0,12,length.out = 500)
Mu=seq(16,50,length.out = 750)
max_Lik<-rep(NA,length(xi_cross))
max_mu<-rep(NA,length(xi_cross))
max_tau<-rep(NA,length(xi_cross))
for(j in 1:length(xi_cross)){
  Par=expand.grid(tau=Tau,mu=Mu)
  Par<-cbind(mu=Par$mu, xi=xi_cross[j], tau=Par$tau)
  lik_Par<-rep(NA,nrow(Par))
  for(i in 1:nrow(Par)){
    tmp<-Lik(Par[i,],Y)
    lik_Par[i]<-tmp
  }
  max_Lik[j]<-max(lik_Par,na.rm=TRUE)
  max_mu[j]<-Par[which.max(lik_Par),1]
  max_tau[j]<-Par[which.max(lik_Par),3]
  cat('j=',j,'\n')
}

pdf("~/Desktop/Illustration/PL_neg.pdf",width=4, height=4)
which<-seq(1,100,by=20)
which <- c(which,101:length(xi_cross))
plot(c(-1,xi_cross[which]),c(n*log(n)-n*log(sum(max(Y)-Y)[-which.max(Y)])-n,max_Lik[which]),type='l',
         xlab=expression(xi), ylab=expression(PL[n](xi)))
abline(v=-0.2,lty=2)
dev.off()

plot(c(-1,xi_cross),c(n*log(n)-n*log(sum(max(Y)-Y)[-which.max(Y)])-n,max_Lik),type='l',
    xlab=expression(xi), ylab=expression(PL[n](xi)))
points(-1,n*log(n)-n*log(sum(max(Y)-Y)[-which.max(Y)])-n,pch=20)
plot(xi_cross,max_mu,type='l',xlab=expression(xi),ylab=expression(mu(xi)))
plot(xi_cross,max_tau,type='l',xlab=expression(xi),ylab=expression(tau(xi)))

range(lik_Par,na.rm = TRUE)
lik_Par[which(lik_Par<quantile(lik_Par,0.25,na.rm = TRUE))]<-quantile(lik_Par,0.25,na.rm = TRUE)

xlim=c(range(Tau)[1]+0.05*diff(range(Tau)),range(Tau)[2]-0.05*diff(range(Tau)))
ylim=c(range(Mu)[1]+0.05*diff(range(Mu)),range(Mu)[2]-0.05*diff(range(Mu)))
x<-cbind(tau=c(0,0,1000,1000,0),mu=c(1000,max(Y),bound(1000,xi_MLE),1000,1000))
# polygon(x=c(xlim[1],xlim[1],xlim[2],xlim[2],xlim[1]),y=c(ylim[1],ylim[2],ylim[2],ylim[1],ylim[1]),lty=2,col=alpha('steelblue',0.05))
plot(x,type='n',xlab=expression(tau),ylab=expression(mu),xlim=xlim,
     ylim=ylim,main=paste("ξ_MLE=",sprintf("%.2f",xi_MLE), ", ξ_cross=",sprintf("%.2f",tail(xi_cross,1)),sep=""))
abline(v=0,col="grey",lty=2)
image(x=Tau,y=Mu,
      z=matrix(lik_Par,length(Tau)),col=terrain.colors(40),add=TRUE)



deriv_on_xi <- rep(NA,length(M))
for(j in 1:length(M)){
  mu_c <- optim_Par_cross_section[j,1]
  xi_c <- optim_Par_cross_section[j,2]
  tau_c<- optim_Par_cross_section[j,3]
  beta_c <-  mu_c - tau_c/xi_c
  V <- (Y-beta_c)^{-1/xi_c}
  deriv_on_xi[j] <- sum(V*log(V))-sum(V)-mean(V)*sum(log(V))
}

plot(optim_Par_cross_section[,2], deriv_on_xi, type='l', lwd=2, col='blue', xlab=expression(xi),ylab="Derivative")
abline(v=xi_MLE,lty=2)


tmp <-cbind(optim_Par_cross_section[,2],optim_Lik_cross_section)
tmp<-rbind(cbind(xi_cross,max_Lik),tmp)
tmp <- rbind(c(-1,-32472.73),tmp)
plot(tmp,type='l',col='red',lwd=2,
     xlab=expression(xi),ylab=expression(PL[n](xi)),main=expression(xi[0]==0.2))

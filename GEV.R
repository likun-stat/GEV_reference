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
Res$par[1]-Res$par[3]/Res$par[2]
mu-tau/xi

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
## ------------------------------------ Illustrations ------------------------------------
##########################################################################################

## ---------------- For positive shape parameter -------------------
mu <- 20
tau <- 0.5
xi <- 0.1
n <- 1000
Y <- revd(n, loc=mu, scale=tau, shape=xi)

library(scales)
bound <- function(tau,xi,alpha=0.1) tau/xi+min(Y)

pdf("~/Desktop/Illustration/illus2.pdf",width=5, height=5)
x<-cbind(tau=c(1000,0,0,1000,1000),mu=c(-1000,-1000,min(Y),bound(1000,xi),-1000))
plot(x,type='n',xlab=expression(tau),ylab='',xlim=c(-0.2,1.5),ylim=c(0.3*min(Y),2.2*min(Y)),yaxt='n',main=expression(xi==hat(xi)[n]))
title(ylab=expression(mu),line=0, cex.lab=1)
abline(v=0,col="grey",lty=2)
text(-0.11,min(Y)-0.05*min(Y),"min(Y)",pos=3)
text(-0.05,mu-tau/xi-0.1*min(Y),expression(hat(beta)[n]),pos=3)

alpha=0.1
lines(x=c(0,-1000),y=c(min(Y),bound(-1000,xi)),col="grey",lty=2)
polygon(x,col=alpha('steelblue',alpha))

x<-cbind(tau=c(1000,0,0,1000,1000),mu=c(-1000,-1000,mu-tau/xi,bound(1000,xi)-min(Y)+mu-tau/xi,-1000))
lines(x=c(0,-1000),y=c(min(Y),bound(-1000,xi)-min(Y)+mu-tau/xi),col="grey",lty=2)
polygon(x,col=alpha('steelblue',alpha+0.3))

points(tau,mu,pch=20)
r=0.08
Theta<-seq(0,2*pi,length.out=500)
circle<-cbind(tau+r*cos(Theta),mu+r*sin(Theta)*23)
points(circle,type='l',lwd=2,col='red')
dev.off()




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
## ------------------------------ Likelihood contour plot---------------------------------
##########################################################################################

# calculate the likelihhood
Lik<-function(par,Y){
  mu<-par[1]; xi<-par[2]; tau<-par[3]
  if(tau<=0) return(NA)
  W <- 1+ (xi)*(Y-mu)/(tau)
  if(any(W<0)) return(NA) else return(-length(Y)*log(tau)-(xi+1)*sum(log(W))/xi-sum(W^{-1/xi}))
}


bound <- function(tau,xi,alpha=0.1) tau/xi+min(Y)

xi_try<-0.2
mu <- 20
tau <- 0.5
n <- 1000
Y <- revd(n, loc=mu, scale=tau, shape=xi_try)
delta<-0.02
mu+delta-(tau+delta)/(xi_try-delta)
Res <- optim(c(mu+delta,xi_try+delta,tau+delta), Lik, Y=Y, method = "L-BFGS-B", lower = c(mu-delta,xi_try-delta,tau-delta), upper = c(mu+delta,xi_try+delta,tau+delta), control=list(fnscale=-1))



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



Tau=seq(0.3,0.7,length.out = 200)
Mu=seq(19.5,20.5,length.out = 200)
# Mu=seq(0.7*max(Y),1.4*max(Y),length.out = 200)
Par=expand.grid(tau=Tau,mu=Mu)
Par<-cbind(mu=Par$mu, xi=xi_try, tau=Par$tau)
lik_Par<-rep(NA,nrow(Par))
hes_Par<-rep(NA,nrow(Par))
for(i in 1:nrow(Par)){
  tmp<-Lik(Par[i,],Y)
  lik_Par[i]<-tmp
  tmp<-hessian(Par[i,],Y,Res)
  hes_Par[i]<-tmp
}
range(lik_Par,na.rm = TRUE)
range(hes_Par,na.rm = TRUE)
lik_Par[which(lik_Par<quantile(lik_Par,0.25,na.rm = TRUE))]<-quantile(lik_Par,0.25,na.rm = TRUE)

xlim=c(range(Tau)[1]+0.05*diff(range(Tau)),range(Tau)[2]-0.05*diff(range(Tau)))
ylim=c(range(Mu)[1]+0.05*diff(range(Mu)),range(Mu)[2]-0.05*diff(range(Mu)))
x<-cbind(tau=c(0,0,1000,1000,0),mu=c(1000,max(Y),bound(1000,xi_try),1000,1000))
# polygon(x=c(xlim[1],xlim[1],xlim[2],xlim[2],xlim[1]),y=c(ylim[1],ylim[2],ylim[2],ylim[1],ylim[1]),lty=2,col=alpha('steelblue',0.05))
plot(x,type='n',xlab=expression(tau),ylab=expression(mu),xlim=xlim,
     ylim=ylim,main=paste("ξ=",sprintf("%.1f",xi_try)))
abline(v=0,col="grey",lty=2)
text(-1.7,max(Y)-0.1*max(Y),"max(Y)",pos=3)

alpha=0.1
x<-cbind(tau=c(0,0,1000,1000,0),mu=c(1000,max(Y),bound(1000,xi_try),1000,1000))
lines(x=c(0,-1000),y=c(max(Y),bound(-1000,xi_try)),col="grey",lty=2)

image(x=Tau,y=Mu,
      z=matrix(lik_Par,length(Tau)),col=terrain.colors(40),add=TRUE)

polygon(x,col=alpha('steelblue',0.05))
points(tau,mu,pch=20)
points(Par[which.max(lik_Par),3],Par[which.max(lik_Par),1],pch=20,col='red')





# Speculation
Tau_1d<-seq(0.35,100,length.out = 200)
Mu_1d<-22.5+Tau_1d/xi
lik_1d<-rep(NA, length(Tau_1d))
dist_1d<-rep(NA, length(Tau_1d))
for(i in 1:length(Tau_1d)){
  # lik_1d[i]<-Lik(c(Mu_1d[i],xi,Tau_1d[i]),Y)-Lik(c(mu,xi,tau),Y)
  lik_1d[i]<--log(tau/Tau_1d[i])/xi+1-(tau/Tau_1d[i])^{-1/xi}#-(tau/Tau_1d[i])^{-1/xi}*{xi*(-22.5)/tau}^{-1/xi}
  dist_1d[i]<-sqrt((Tau_1d[i]-tau)^2+(Mu_1d[i]-mu)^2)
}
plot(Tau_1d,lik_1d,type='l',xlab=expression(tau),ylab='Lik',ylim=c(-25000,0))
points(Tau_1d,-dist_1d,type='l')


plot(Tau_1d[-(1:3)],lik_1d[-(1:3)],type='l',xlab=expression(tau),ylab='Lik',ylim=c(-100,0))
points(Tau_1d,-dist_1d^0.01,type='l')



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
tmp <- function(xi) if(xi<0.28) return(NA) else sqrt(2^{2/xi-1.5}*xi^{xi-2}*gamma(xi+1))
tmp<-Vectorize(tmp)
ggplot(data.frame(xi=c(-0.49, 4)), aes(xi)) + 
  stat_function(fun=H33, alpha=0.7, size=1.1, aes(color = "line1")) +
  stat_function(fun=tmp, linetype="dashed", col='#0a0001',alpha=0.9) +
  # stat_function(fun=fun1, linetype="dashed", col='#0a0001',alpha=0.9) +
  xlab(expression(xi))+ylab(expression(sqrt(h[33])))+
  coord_cartesian(xlim=c(-0.49, 4), ylim=c(0.5,15)) +
  scale_colour_manual("Lgend title", values = c("line1"='#d113be'),labels=expression(paste("(",mu,', ',tau,', ', xi,')'))) +
  # scale_y_continuous(breaks=c(0,sqrt(h33),3,6,9,12,15),labels = c(0,'',3,6,9,12,15)) +
  scale_x_continuous(breaks =c(-0.52,0,1,2,3,4),labels = c(-0.5,0,1,2,3,4)) +
  theme(panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major.x = element_line(colour = "grey95"),
        panel.grid.major.y = element_line(colour = "grey95"),
        legend.position = c(0.95,0.95), legend.justification = c("right","top"), legend.title=element_blank(), 
        legend.text=element_text(size=15), legend.spacing.x = unit(0.2,'cm'),
        legend.box.background = element_rect(colour = "black",size=0.7),
        legend.text.align=0.5,
        plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"))
ggsave("~/Desktop/Illustration/h33_comp.pdf",width = 4, height = 4)






##########################################################################################
## ---------------------------------- Assorted plots -------------------------------------
##########################################################################################

f1 <- function(tau) {xi0=2; tau0=10; -log(tau0/tau)/xi0-(tau0/tau)^{-1/xi0}+1}
curve(f1, from=0, to=20, ylab=expression(L[n](theta)-L[n](hat(theta)[n])),xlab=expression(tau))
ggplot(data.frame(xi=c(0, 20)), aes(xi)) + 
  stat_function(fun=f1, col='#036180',alpha=0.5, size=1.1) +
  geom_vline(xintercept=10, linetype='dashed') +
  xlab(expression(tau))+ylab(expression(L[n](theta)-L[n](hat(theta)[n])))+
  coord_cartesian(xlim=c(-0.51, 20), ylim=c(0,-8)) +
  scale_y_continuous(breaks =c(0,-8),labels = c('0',expression(-infinity))) +
  scale_x_continuous(breaks =c(-0.4,10,20),labels = c('0',expression(hat(tau)[n]),expression(infinity))) +
  theme(panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA),
        # panel.grid.major.x = element_line(colour = "grey95"),
        # panel.grid.major.y = element_line(colour = "grey95"),
        legend.position = 'none', legend.justification = c("right","top"), legend.title=element_blank(), 
        legend.text=element_text(size=10), legend.spacing.x = unit(0.2,'cm'),
        axis.text=element_text(size=11,color=alpha('black',0.6)),
        legend.box.background = element_rect(colour = "black",size=0.7),
        legend.text.align=0.5,
        plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"))
ggsave("~/Desktop/Illustration/caseI1.pdf",width = 4, height = 4)

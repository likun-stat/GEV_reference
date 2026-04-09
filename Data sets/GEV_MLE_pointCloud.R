##########################################################################################
## ----------------------------- GEV maximum likelihood ----------------------------------
##########################################################################################
## MAIN IDEA:
## 1. PL_n(xi) stands for profile likelihood, which is the maximum of the GEV joint likelihood 
##     over the (mu, tau) surface while fixing the shape xi.
## 2. We can prove that the GEV joint likelihood is a convex function of (mu, tau) while 
##     fixing xi. Therefore, it is better to first calculate the profile likelihood, and then
##     optimize PL_n(xi) over xi. 



## -------------------------------
##       Getting ready
## -------------------------------
# calculate the GEV joint likelihhood
Lik<-function(par,Y){
  mu<-par[1]; xi<-par[2]; tau<-par[3]
  if(tau<=0) return(NA)
  if(xi==0) return(-length(Y)*log(tau)-sum((Y-mu)/tau)-sum(exp(-(Y-mu)/tau)))
  W <- 1+ (xi)*(Y-mu)/(tau)
  if(any(W<0)) return(NA) else return(-length(Y)*log(tau)-(xi+1)*sum(log(W))/xi-sum(W^{-1/xi}))
}


bound_min <- function(tau,xi,Y) tau/xi+min(Y)
bound_max <- function(tau,xi,Y) tau/xi+max(Y)
## Value initial value for tau when fixing mu and xi
tau_initial_value <- function(mu=0, xi, Y){
  if(xi>0) tau_init <- max(0.1, xi*(mu-min(Y))+0.01)
  if(xi<0) tau_init <- max(0.1, xi*(mu-max(Y))+0.01)
  if(xi==0) tau_init <- 0.1
  return(tau_init)
}


## -------------------------------
##    Simulate Y: Positive shape
## -------------------------------
set.seed(12)
xi_true<-0.2
mu <- 0
tau <- 1
n <- 1000  #the number of seasonal maxima
MLEs <- matrix(NA, nrow=1000, ncol=3)
for(iter in 1:nrow(MLEs)){
  Y <- extRemes::revd(n, loc=mu, scale=tau, shape=xi_true)
  res_tmp <- optim(par=c(0, 0.3, tau_initial_value(0, 0.3, Y)), Lik, Y=Y, control=list(fnscale=-1))
  MLEs[iter,] <- res_tmp$par
}

plot_data <- data.frame(mu = MLEs[,1], xi = MLEs[,2], tau = MLEs[,3])
plot1 <- ggplot(plot_data) + geom_point(aes(x=mu, y=xi), alpha=0.3) +
  geom_point(aes(x=0, y=0.2), colour="red", size=7, shape="+") +
  labs(x= expression(mu), y = expression(xi), title = paste("sample size = ", n)) +
  scale_x_continuous(expand = c(0, 0), limits=c(-1,1)) +
  scale_y_continuous(expand = c(0, 0), limits=c(-0.5, 0.9))  +
  theme(plot.title = element_text(face="bold", hjust=0.5))
plot1


plot2 <- ggplot(plot_data) + geom_point(aes(x=mu, y=tau), alpha=0.3) +
  geom_point(aes(x=0, y=1), colour="red", size=7, shape="+") +
  labs(x= expression(mu), y = expression(tau), title = paste("sample size = ", n)) +
  scale_x_continuous(expand = c(0, 0), limits=c(-1,1)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 2))  +
  theme(plot.title = element_text(face="bold", hjust=0.5))
plot2






zeta <- 1.20205
gam <- -digamma(1)

dev2 <- gam^2+pi^2/6
dev3 <- -gam^3-3*gam*pi^2/6-2*zeta
dev4 <- gam^4+6*gam^2*pi^2/6+8*gam*zeta+27*pi^4/(90*2)


##-------------------------------------------------------
##-------------------- Beta Prior -----------------------
##-------------------------------------------------------

beta_prior<-function(xi){
  if(xi>=0.5) return(-Inf)
  if(xi<=-0.5) return(-Inf)
  log((0.5+xi)^{5}*(0.5-xi)^{8}/beta(6,9))
  }
# curve(beta_prior, from=-0.5, to=0.5)



##-------------------------------------------------------
##------------------- (xi,tau,mu) -----------------------
##-------------------------------------------------------
H11<- function(xi){
  if(xi==-1/2) return(log(sqrt(2*pi^2/3)))
  if(xi==0) return(log(sqrt(11*pi^4/360-6*zeta^2/pi^2)))
  
  p<- (1+xi)^2*gamma(2*xi+1)
  q<- gamma(2+xi)*(digamma(1+xi)+1/xi+1)
  s<-1-0.5772157+1/xi
  # r<-q-p/xi
  tmp <- (pi^2/6)-(q-s*gamma(xi+2))^2/(p-gamma(xi+2)^2)
  h11<-sqrt(tmp)/abs(xi)
  return(log(h11))
}
x <- c(seq(-0.4999,-0.03,by=0.001),0,seq(0.025,10,by=0.001))
H11 <- Vectorize(H11)
h_approx<-approxfun(x=c(-0.5,x),y=c(H11(-0.5),H11(x)))
h11<-function(xi) {
  if(xi<=-0.5) return(-Inf)
  else if(xi <10) return(h_approx(xi))
  else return(H11(xi))}

# curve(h11, from=-0.3, to=0.3)
# curve(h11, from=-0.5, to=10)



##-------------------------------------------------------
##------------------- (mu,xi,tau) -----------------------
##-------------------------------------------------------
H22_1 <- function(xi){
  if(xi==-1/2) return(log(sqrt(2*pi^2/3+4*(1-gam)^2)))
  if(xi==0) return(log(sqrt(dev2+dev3+dev4/4-(-gam+3*dev2/2+dev3/2)^2/((1-gam)^2+pi^2/6))))
  p<- (1+xi)^2*gamma(2*xi+1)
  q<- gamma(2+xi)*(digamma(1+xi)+1/xi+1)
  s<-1-0.5772157+1/xi
  # r<-q-p/xi
  tmp <- (pi^2/6)+(1-0.5772157)^2-(gamma(2+xi)/xi-q+1-0.5772157)^2/(1+p-2*gamma(xi+2))
  h22<-sqrt(tmp)/abs(xi)
  return(log(h22))
}
x <- c(seq(-0.4999,-0.03,by=0.001),0,seq(0.025,10,by=0.001))
H22_1 <- Vectorize(H22_1)
h_approx<-approxfun(x=c(-0.5,x),y=c(H22_1(-0.5),H22_1(x)))
h22_1<-function(xi){
  if(xi<=-0.5) return(-Inf)
  else if(xi <10) return(h_approx(xi))
  else return(H22_1(xi))}

# curve(h22_1, from=-0.3, to=0.3)
# curve(h22_1, from=-0.5, to=10,ylim=c(-2,2))



##-------------------------------------------------------
##------------------- (tau,xi,mu) -----------------------
##-------------------------------------------------------
H22_2<- function(xi){
  if(xi+1/2<1e-4) return(log(sqrt(2*pi^2/3+4*(1+gam)^2)))
  if(xi==0) return(log(sqrt(pi^2*(1-gam)^2/6+2*(gam-1)*zeta+11*pi^4/360)))
  p<- (1+xi)^2*gamma(2*xi+1)
  q<- gamma(2+xi)*(digamma(1+xi)+1/xi+1)
  s<-1-0.5772157+1/xi
  # r<-q-p/xi
  tmp <- (pi^2/6)+s^2-q^2/p
  h22<-sqrt(tmp)/abs(xi)
  return(log(h22))
}
x <- c(seq(-0.4999,-0.03,by=0.001),0,seq(0.025,10,by=0.001))
H22_2 <- Vectorize(H22_2)
h_approx<-approxfun(x=c(-0.5,x),y=c(H22_2(-0.5),H22_2(x)))
h22_2<-function(xi){
  if(xi<=-0.5) return(-Inf)
  else if(xi <10) return(h_approx(xi))
  else return(H22_1(xi))}

# curve(h22_2, from=-0.3, to=0.3)
# curve(h22_2, from=-0.5, to=10, col='red', add=TRUE)




##-------------------------------------------------------
##----------------------- MDI ---------------------------
##-------------------------------------------------------
MDI<-function(xi){
  if(xi<=-0.5) return(-Inf)
  -gam*(1+xi)
}



##-------------------------------------------------------
##------------------- Jefferys split --------------------
##-------------------------------------------------------
jeffreys_split <- function(xi){
  p<- (1+xi)^2*gamma(2*xi+1)
  q<- gamma(2+xi)*(digamma(1+xi)+1/xi+1)
  s<-1-0.5772157+1/xi
  # r<-q-p/xi
  tmp <- p*(1+p-2*gamma(xi+2))*(pi^2/6+s^2-2*q/xi+p/xi^2)
  jef<-sqrt(tmp)/xi^2
  return(log(jef))
}



##-------------------------------------------------------
##------------------- Prior for tau ---------------------
##-------------------------------------------------------
tau_prior <- function(tau){
  if(tau<=0) return(-Inf) else return(-log(tau))
}


##-------------------------------------------------------
##------------------- Prior for mu ---------------------
##-------------------------------------------------------
mu_prior <- function(mu){
  return(0)
}


##----------------------------------------------------------
##------------------- Joint likelihood ---------------------
##----------------------------------------------------------
# For optim()
Lik<-function(par,Y){
  mu<-par[1]; xi<-par[2]; tau<-par[3]
  if(tau<=0) return(NA)
  W <- 1+ (xi)*(Y-mu)/(tau)
  if(any(W<0)) return(NA) else return(-length(Y)*log(tau)-(xi+1)*sum(log(W))/xi-sum(W^{-1/xi}))
}


Lik_beta <- function(par,Y){
  xi<-par[2]; tau<-par[3]
  return(Lik(par,Y)+beta_prior(xi = xi)-log(tau))
}


Lik_h11 <- function(par,Y){
  xi<-par[2]; tau<-par[3]
  return(Lik(par,Y)+h11(xi)-log(tau))
}

Lik_h22_1 <- function(par,Y){
  xi<-par[2]; tau<-par[3]
  return(Lik(par,Y)+h22_1(xi)-log(tau))
}

Lik_h22_2 <- function(par,Y){
  xi<-par[2]; tau<-par[3]
  return(Lik(par,Y)+h22_2(xi)-log(tau))
}

Lik_MDI <- function(par,Y){
  xi<-par[2]; tau<-par[3]
  return(Lik(par,Y)+MDI(xi = xi)-log(tau))
}


# For Metropolis
Lik2<-function(Y,mu,xi,tau){
  if(tau<=0) return(-Inf)
  W <- 1+ (xi)*(Y-mu)/(tau)
  if(xi==0) {
    tmp <- (Y-mu)/tau
    return(-length(Y)*log(tau)-sum(tmp)-sum(exp(-tmp)))}
  if(any(W<0)) return(-Inf) else 
    return(-length(Y)*log(tau)-(xi+1)*sum(log(W))/xi-sum(W^{-1/xi}))
}


##---------------------------------------------------------------------
##------------------- Generate GEV random samples ---------------------
##---------------------------------------------------------------------
revd <- function (n, loc = 0, scale = 1, shape = 0, threshold = 0, type = c("GEV", "GP")) 
{
  type <- match.arg(type)
  type <- tolower(type)
  if (type == "gev") 
    z <- rexp(n)
  else if (type == "gp") {
    z <- runif(n)
    loc <- threshold
  }
  else stop("revd: invalid type argument.")
  out <- numeric(n) + NA
  if (min(scale) < 0) 
    stop("revd: scale parameter(s) must be positively valued.")
  if (length(loc) == 1) 
    loc <- rep(loc, n)
  if (length(scale) == 1) 
    sc <- rep(scale, n)
  else sc <- scale
  if (length(shape) == 1) 
    sh <- rep(shape, n)
  else sh <- shape
  if (length(loc) != n || length(sc) != n || length(sh) != 
      n) 
    stop("revd: parameters must have length equal to 1 or n.")
  id <- sh == 0
  if (any(id)) {
    if (type == "gev") 
      out[id] <- loc[id] - sc[id] * log(z[id])
    else if (type == "gp") 
      out[id] <- loc[id] + rexp(n, rate = 1/sc[id])
  }
  if (any(!id)) 
    out[!id] <- loc[!id] + sc[!id] * (z[!id]^(-sh[!id]) - 
                                        1)/sh[!id]
  return(out)
}




##-----------------------------------------------------------------
##------------------- Priors in the same plot ---------------------
##-----------------------------------------------------------------
library(ggplot2)
prior_plot<-ggplot(data.frame(xi=c(-0.5, 0.5)), aes(xi)) + 
  stat_function(fun=h11, aes(color ='line1'),alpha=0.7, size=1.1) +
  stat_function(fun=h22_1, aes(color ='line2'),alpha=0.7, size=1.1)+
  stat_function(fun=h22_2, aes(color ='line3'),alpha=0.7, size=1.1)+
  stat_function(fun=MDI,  aes(color ='line4'),alpha=0.9) +
  xlab(expression(xi))+ylab(expression(pi(xi)))+
  coord_cartesian(xlim=c(-0.5, 0.5)) +
  # scale_y_continuous(breaks=c(0,1,sqrt(h22),2,sqrt(h22_1),3,4),labels = c(0,1,'',2,'',3,4)) +
  # scale_x_continuous(breaks =c(-0.5,0,2.5,5,7.5,10),labels = c(-0.5,'  0',2.5,5,7.5,10)) +
  scale_colour_manual("Lgend title", values = c("line1"='#036180', 'line2'='#deaf04','line3'='#2c20b0','line4'='red'),
                      labels=c(expression(paste("(",xi,', ',tau,', ', mu,')')),
                               expression(paste("(",mu,', ',xi,', ', tau,')')),
                               expression(paste("(",tau,', ',xi,', ', mu,')')),
                               'MDI')) +
  theme(panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major.x = element_line(colour = "grey95"),
        panel.grid.major.y = element_line(colour = "grey95"),
        legend.position = c(0.95,0.95), legend.justification = c("right","top"), legend.title=element_blank(), 
        legend.text=element_text(size=10), legend.spacing.x = unit(0.2,'cm'),
        legend.box.background = element_rect(colour = "black",size=0.7),
        legend.text.align=0.5,
        plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"))

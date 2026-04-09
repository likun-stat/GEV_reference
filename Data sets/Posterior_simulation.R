setwd("~/Desktop/Paper Submissions/Statistica Sinica & EJS & SJS & CJS (Pos_Nor)/Simulation_output_figures/")
source("../../../Past research projects/Research/GEV_reference/generic_samplers.R")
source("../../../Past research projects/Research/GEV_reference/ReferencePrior_utils.R")
source("../../../Past research projects/Research/GEV_reference/ReferencePrior_sampler.R")


set.seed(123)
N=1000
true_Par <- c(10, 0.2, 2)
Current_data <- extRemes::revd(N, loc = true_Par[1], scale = true_Par[3], shape = true_Par[2])
id <- gsub(".","_", true_Par[2], fixed=TRUE)

Tau <- seq(0.1,10,length.out = 35)
Mu <- seq(-20,20,length.out = 35)
Xi <- seq(-0.5,0.5,length.out = 35)
Par=expand.grid(tau=Tau,mu=Mu,xi=Xi)
Par<-cbind(mu=Par$mu, xi=Par$xi, tau=Par$tau)
lik_Par <- rep(NA,nrow(Par))

xi_prior_type="h11" #xi,mu,tau
# xi_prior_type="h22_2" #tau,xi,mu
# xi_prior_type="MDI"
# xi_prior_type="Beta"

for(i in 1:nrow(Par)){
  tmp <- Lik(Par[i,],Current_data)
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


## Initial setup
n.updates <- 10000
thin <- 10
echo.interval <- 1000

## Run MCMC
initial.values <- list(mu=mu_init, xi=xi_init, tau=tau_init)
true_params <- list(mu = true_Par[1], xi = true_Par[2], tau = true_Par[3])
res <- GEV.sampler.02(Y=Current_data, iter.num = id,
                      xi_prior_type=xi_prior_type,
                      initial.values=initial.values,
                      n.updates=n.updates, thin=thin,
                      experiment.name="station",
                      echo.interval=echo.interval,
                      true.params=true_params, sd.ratio=NULL, lower.prob.lim=0.5)

choose_pos_num <- 1000
mu_vec <- tail(res$mu.trace, choose_pos_num)
tau_vec <- tail(res$tau.trace, choose_pos_num)
xi_vec <- tail(res$xi.trace, choose_pos_num)


# setwd("~/Desktop/Statistica Sinica (Pos_Nor)/Simulation_output_figures_neg_zero")
## Visualize
pal <- hcl.colors(10, "Spectral")[10:1]


if(N<500) h=0.3 else h=0.2
f1 <- MASS::kde2d(mu_vec, tau_vec, h=h, n = 200, lims = c(9,11.5,1.5,2.5))
f2 <- MASS::kde2d(mu_vec, xi_vec, h=h, n = 200, lims = c(9,11.5,-0.25,0.4))
f3 <- MASS::kde2d(tau_vec, xi_vec,h=h, n = 200, lims = c(1.5,2.5,-0.25,0.4))
upper <- max(f1$z, f2$z, f3$z)

plot_dat1 <- expand.grid(x=f1$x, y=f1$y)
plot_dat1$den <- as.vector(f1$z)
plot1 <- ggplot(plot_dat1) + geom_raster(aes(x=x, y=y, fill=den)) +
  geom_contour(col="black", aes(x=x, y=y, z=den)) +
  labs(x= expression(mu), y = expression(tau)) +
  scale_x_continuous(expand = c(0, 0), limits=c(9,11.5)) +
  scale_y_continuous(expand = c(0, 0), limits=c(1.5,2.5)) +
  scale_fill_gradientn(colours = pal, name = "Density", na.value = NA, limits = c(0,36.9)) 

plot_dat2 <- expand.grid(x=f2$x, y=f2$y)
plot_dat2$den <- as.vector(f2$z)
plot2 <- ggplot(plot_dat2) + geom_raster(aes(x=x, y=y, fill=den)) +
  geom_contour(col="black", aes(x=x, y=y, z=den)) +
  labs(x= expression(mu), y = expression(xi)) +
  scale_x_continuous(expand = c(0, 0), limits=c(9,11.5)) +
  scale_y_continuous(expand = c(0, 0), limits=c(-0.25,0.4)) +
  scale_fill_gradientn(colours = pal, name = "Density", na.value = NA, limits = c(0,36.9)) 

plot_dat3 <- expand.grid(x=f3$x, y=f3$y)
plot_dat3$den <- as.vector(f3$z)
plot3 <- ggplot(plot_dat3) + geom_raster(aes(x=x, y=y, fill=den)) +
  geom_contour(col="black", aes(x=x, y=y, z=den)) +
  labs(x= expression(tau), y = expression(xi)) +
  scale_x_continuous(expand = c(0, 0), limits=c(1.5,2.5)) +
  scale_y_continuous(expand = c(0, 0), limits=c(-0.25,0.4)) +
  scale_fill_gradientn(colours = pal, name = paste0("Density\n(n = ",N, ")"), na.value = NA, limits = c(0,36.9))
  # scale_fill_gradientn(colours = pal, name = paste0("Density\n(prior = ",xi_prior_type, ")"), na.value = NA, limits = c(0,upper))


library(cowplot)
prow <- plot_grid(plot1+ theme(legend.position="none"), plot2+ theme(legend.position="none"), plot3+ theme(legend.position="none"),
          align = "h", nrow = 1)
legend <- get_legend(
  # create some space to the left of the legend
  plot3 + theme(legend.box.margin = margin(0, 0, 0, 0))
)


plot_grid(prow, legend, rel_widths = c(3, .4))
ggsave(file = paste0("./posterior_density_xi_", id, "_N_", N, ".pdf"), width = 8, height =2.5)
# ggsave(file = paste0("./posterior_density_xi_", id, "_N_", N, "_neg.pdf"), width = 8, height =2.5)
# ggsave(file = paste0("./posterior_density_xi_", id, "_prior_", xi_prior_type, ".pdf"), width = 8, height =2.5)

## Normality test
library(energy)
pos_data <- data.frame(tail(mu_vec,200), tail(xi_vec,200),tail(tau_vec,200))
mvnorm.etest(pos_data, R=1000)
mvnorm.etest(pos_data[,1], R=1000)
mvnorm.etest(pos_data[,2], R=1000)
mvnorm.etest(pos_data[,3], R=1000)

pos_data <- data.frame(mu_vec, xi_vec,tau_vec)
mvnorm.etest(pos_data, R=1000)
mvnorm.etest(pos_data[,1], R=1000)
mvnorm.etest(pos_data[,2], R=1000)
mvnorm.etest(pos_data[,3], R=1000)




save(Current_data, mu_vec, tau_vec, xi_vec, file = paste0("./posterior_sample_xi_", id, "_N_", N, ".RData"))


## calculate the MLE
Lik<-function(par,Y){
  mu<-par[1]; xi<-par[2]; tau<-par[3]
  if(tau<=0) return(-1e10)
  if(xi==0) return(-length(Y)*log(tau)-sum((Y-mu)/tau)-sum(exp(-(Y-mu)/tau)))
  W <- 1+ (xi)*(Y-mu)/(tau)
  if(any(W<0)) return(-1e10) else return(-length(Y)*log(tau)-(xi+1)*sum(log(W))/xi-sum(W^{-1/xi}))
}

tau_mat <- array(NA, dim=c(1000,4))
mu_mat <- array(NA, dim=c(1000,4))
xi_mat <- array(NA, dim=c(1000,4))
tau_mat_mle <- array(NA, dim=c(1000,4))
mu_mat_mle <- array(NA, dim=c(1000,4))
xi_mat_mle <- array(NA, dim=c(1000,4))
counter = 1
for(N in c(50, 100, 500, 1000)){
  load(paste0("./posterior_sample_xi_", id, "_N_", N, ".RData"))
  delta <-0.1
  Res <- optim(c(true_Par[1],true_Par[2],true_Par[3]), Lik, Y=Current_data, method = "L-BFGS-B", 
               lower = c(true_Par[1]-delta,true_Par[2]-delta,true_Par[3]-delta), 
               upper = c(true_Par[1]+delta,true_Par[2]+delta,true_Par[3]+delta), control=list(fnscale=-1))
  cov <- normal_Cov_GEV(mu = true_Par[1], xi = true_Par[2], tau = true_Par[3], n=N)
  tau_mat[,counter] <- (tau_vec-true_Par[3])/sqrt(cov[2,2])
  mu_mat[,counter] <- (mu_vec-true_Par[1])/sqrt(cov[1,1])
  xi_mat[,counter] <- (xi_vec-true_Par[2])/sqrt(cov[3,3])
  
  tau_mat_mle[,counter] <- (tau_vec-Res$par[3])/sqrt(cov[2,2])
  mu_mat_mle[,counter] <- (mu_vec-Res$par[1])/sqrt(cov[1,1])
  xi_mat_mle[,counter] <- (xi_vec-Res$par[2])/sqrt(cov[3,3])
  
  counter <- counter + 1
}


## ------------------------ qqplots ------------------------
pdf(paste0("./QQ_mu_", id, ".pdf"), width = 10, height =2.5)
par(mfrow=c(1,4))
qqnorm(mu_mat[,1], main = paste0("n = ", 50), ylim=c(-4,4), xlim=c(-4,4))
qqline(mu_mat[,1], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(mu_mat[,2], main = paste0("n = ", 100), ylim=c(-4,4), xlim=c(-4,4))
qqline(mu_mat[,2], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(mu_mat[,3], main = paste0("n = ", 500), ylim=c(-4,4), xlim=c(-4,4))
qqline(mu_mat[,3], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(mu_mat[,4], main = paste0("n = ", 1000), ylim=c(-4,4), xlim=c(-4,4))
qqline(mu_mat[,4], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))
par(mfrow=c(1,1))
dev.off()


pdf(paste0("./QQ_tau_", id, ".pdf"), width = 10, height =2.5)
par(mfrow=c(1,4))
qqnorm(tau_mat[,1], main = paste0("n = ", 50), ylim=c(-4,4), xlim=c(-4,4))
qqline(tau_mat[,1], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(tau_mat[,2], main = paste0("n = ", 100), ylim=c(-4,4), xlim=c(-4,4))
qqline(tau_mat[,2], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(tau_mat[,3], main = paste0("n = ", 500), ylim=c(-4,4), xlim=c(-4,4))
qqline(tau_mat[,3], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(tau_mat[,4], main = paste0("n = ", 1000), ylim=c(-4,4), xlim=c(-4,4))
qqline(tau_mat[,4], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))
par(mfrow=c(1,1))
dev.off()

pdf(paste0("./QQ_xi_", id, ".pdf"), width = 10, height =2.5)
par(mfrow=c(1,4))
qqnorm(xi_mat[,1], main = paste0("n = ", 50), ylim=c(-4,4), xlim=c(-4,4))
qqline(xi_mat[,1], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(xi_mat[,2], main = paste0("n = ", 100), ylim=c(-4,4), xlim=c(-4,4))
qqline(xi_mat[,2], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(xi_mat[,3], main = paste0("n = ", 500), ylim=c(-4,4), xlim=c(-4,4))
qqline(xi_mat[,3], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(xi_mat[,4], main = paste0("n = ", 1000), ylim=c(-4,4), xlim=c(-4,4))
qqline(xi_mat[,4], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))
par(mfrow=c(1,1))
dev.off()


## ------------------------ qqplots_mle ------------------------
pdf(paste0("./QQ_mu_", id, "_mle.pdf"), width = 10, height =2.5)
par(mfrow=c(1,4))
qqnorm(mu_mat_mle[,1], main = paste0("n = ", 50), ylim=c(-4,4), xlim=c(-4,4))
qqline(mu_mat_mle[,1], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(mu_mat_mle[,2], main = paste0("n = ", 100), ylim=c(-4,4), xlim=c(-4,4))
qqline(mu_mat_mle[,2], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(mu_mat_mle[,3], main = paste0("n = ", 500), ylim=c(-4,4), xlim=c(-4,4))
qqline(mu_mat_mle[,3], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(mu_mat_mle[,4], main = paste0("n = ", 1000), ylim=c(-4,4), xlim=c(-4,4))
qqline(mu_mat_mle[,4], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))
par(mfrow=c(1,1))
dev.off()


pdf(paste0("./QQ_tau_", id, "_mle.pdf"), width = 10, height =2.5)
par(mfrow=c(1,4))
qqnorm(tau_mat_mle[,1], main = paste0("n = ", 50), ylim=c(-4,4), xlim=c(-4,4))
qqline(tau_mat_mle[,1], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(tau_mat_mle[,2], main = paste0("n = ", 100), ylim=c(-4,4), xlim=c(-4,4))
qqline(tau_mat_mle[,2], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(tau_mat_mle[,3], main = paste0("n = ", 500), ylim=c(-4,4), xlim=c(-4,4))
qqline(tau_mat_mle[,3], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(tau_mat_mle[,4], main = paste0("n = ", 1000), ylim=c(-4,4), xlim=c(-4,4))
qqline(tau_mat_mle[,4], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))
par(mfrow=c(1,1))
dev.off()

pdf(paste0("./QQ_xi_", id, "_mle.pdf"), width = 10, height =2.5)
par(mfrow=c(1,4))
qqnorm(xi_mat_mle[,1], main = paste0("n = ", 50), ylim=c(-4,4), xlim=c(-4,4))
qqline(xi_mat_mle[,1], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(xi_mat_mle[,2], main = paste0("n = ", 100), ylim=c(-4,4), xlim=c(-4,4))
qqline(xi_mat_mle[,2], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(xi_mat_mle[,3], main = paste0("n = ", 500), ylim=c(-4,4), xlim=c(-4,4))
qqline(xi_mat_mle[,3], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(xi_mat_mle[,4], main = paste0("n = ", 1000), ylim=c(-4,4), xlim=c(-4,4))
qqline(xi_mat_mle[,4], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))
par(mfrow=c(1,1))
dev.off()










## xi = 0.2, N = 1000
c(0.124, 0.412, 0.157, 0.394) 
c(0.581, 0.453, 0.55,  0.323)

## xi = 0.2, N = 500
c(0.032, 0.801, 0.058, 0.322)
c(0.604, 0.873, 0.892, 0.418)

## xi = 0.2, N = 100
c(2.2e-16, 0.082, 0.005, 2.2e-16)
c(2.2e-16, 0.282, 0.132, 0.004)

## xi = 0.2, N = 50
c(2.2e-16, 0.568, 2.2e-16, 2.2e-16)
c(2.2e-16, 0.663, 0.002, 0.311)





## ------ For neg xi --------
pdf(paste0("./QQ_neg.pdf"), width = 7.5, height =2.5)
par(mfrow=c(1,3))
qqnorm(tau_mat[,3], main = expression(paste(tau, " (n = 500, ", xi[0]==-0.2,")")), ylim=c(-4,4), xlim=c(-4,4))
qqline(tau_mat[,3], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))


qqnorm(mu_mat[,3], main = expression(paste(mu, " (n = 500, ", xi[0]==-0.2,")")), ylim=c(-4,4), xlim=c(-4,4))
qqline(mu_mat[,3], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(xi_mat[,3], main = expression(paste(xi, " (n = 500, ", xi[0]==-0.2,")")), ylim=c(-4,4), xlim=c(-4,4))
qqline(xi_mat[,3], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))
par(mfrow=c(1,1))
dev.off()

pdf(paste0("./QQ_neg_mle.pdf"), width = 7.5, height =2.5)
par(mfrow=c(1,3))
qqnorm(tau_mat_mle[,3], main = expression(paste(tau, " (n = 500, ", xi[0]==-0.2,")")), ylim=c(-4,4), xlim=c(-4,4))
qqline(tau_mat_mle[,3], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))


qqnorm(mu_mat_mle[,3], main = expression(paste(mu, " (n = 500, ", xi[0]==-0.2,")")), ylim=c(-4,4), xlim=c(-4,4))
qqline(mu_mat_mle[,3], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(xi_mat_mle[,3], main = expression(paste(xi, " (n = 500, ", xi[0]==-0.2,")")), ylim=c(-4,4), xlim=c(-4,4))
qqline(xi_mat_mle[,3], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))
par(mfrow=c(1,1))
dev.off()


## ------ For zero xi --------
pdf(paste0("./QQ_zero.pdf"), width = 7.5, height =2.5)
par(mfrow=c(1,3))
qqnorm(tau_mat[,3], main = expression(paste(tau, " (n = 500, ", xi[0]==0,")")), ylim=c(-4,4), xlim=c(-4,4))
qqline(tau_mat[,3], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))


qqnorm(mu_mat[,3], main = expression(paste(mu, " (n = 500, ", xi[0]==0,")")), ylim=c(-4,4), xlim=c(-4,4))
qqline(mu_mat[,3], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(xi_mat[,3], main = expression(paste(xi, " (n = 500, ", xi[0]==0,")")), ylim=c(-4,4), xlim=c(-4,4))
qqline(xi_mat[,3], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))
par(mfrow=c(1,1))
dev.off()

pdf(paste0("./QQ_zero_mle.pdf"), width = 7.5, height =2.5)
par(mfrow=c(1,3))
qqnorm(tau_mat_mle[,3], main = expression(paste(tau, " (n = 500, ", xi[0]==0,")")), ylim=c(-4,4), xlim=c(-4,4))
qqline(tau_mat_mle[,3], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))


qqnorm(mu_mat_mle[,3], main = expression(paste(mu, " (n = 500, ", xi[0]==0,")")), ylim=c(-4,4), xlim=c(-4,4))
qqline(mu_mat_mle[,3], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(xi_mat_mle[,3], main = expression(paste(xi, " (n = 500, ", xi[0]==0,")")), ylim=c(-4,4), xlim=c(-4,4))
qqline(xi_mat_mle[,3], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))
par(mfrow=c(1,1))
dev.off()





## ------ For one xi --------
pdf(paste0("./QQ_one.pdf"), width = 7.5, height =2.5)
par(mfrow=c(1,3))
qqnorm(tau_mat[,1], main = expression(paste(tau, " (n = 500, ", xi[0]==1,")")), ylim=c(-4,4), xlim=c(-4,4))
qqline(tau_mat[,1], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))


qqnorm(mu_mat[,1], main = expression(paste(mu, " (n = 500, ", xi[0]==1,")")), ylim=c(-4,4), xlim=c(-4,4))
qqline(mu_mat[,1], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(xi_mat[,1], main = expression(paste(xi, " (n = 500, ", xi[0]==1,")")), ylim=c(-4,4), xlim=c(-4,4))
qqline(xi_mat[,1], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))
par(mfrow=c(1,1))
dev.off()

pdf(paste0("./QQ_one_mle.pdf"), width = 7.5, height =2.5)
par(mfrow=c(1,3))
qqnorm(tau_mat_mle[,1], main = expression(paste(tau, " (n = 500, ", xi[0]==1,")")), ylim=c(-4,4), xlim=c(-4,4))
qqline(tau_mat_mle[,1], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))


qqnorm(mu_mat_mle[,1], main = expression(paste(mu, " (n = 500, ", xi[0]==1,")")), ylim=c(-4,4), xlim=c(-4,4))
qqline(mu_mat_mle[,1], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))

qqnorm(xi_mat_mle[,1], main = expression(paste(xi, " (n = 500, ", xi[0]==1,")")), ylim=c(-4,4), xlim=c(-4,4))
qqline(xi_mat_mle[,1], col='red', lwd=2)
abline(a = 0, b = 1, lwd=2, lty=3, col='blue')
legend("bottomright", lty=c(1,3), col=c('red', 'blue'), legend=c("QQ line", "1-1 line"))
par(mfrow=c(1,1))
dev.off()


source("~/Desktop/Research/GEV_reference/ReferencePrior_utils.R")

experiment.name <- "Beta-prior"
n.experiments <- 100
xi <- 0.15
which <- 5001:10000

L50 <- qevd(0.98, loc=0, shape=xi, scale=1)
L100 <- qevd(0.99, loc=0, shape=xi, scale=1)
returnLevels_50 <- rep(NA, length(which))
returnLevels_100 <- rep(NA, length(which))
QuantileScore50 <- rep(NA, length(which))
QuantileScore100 <- rep(NA, length(which))

Coverage05 <- matrix(NA, nrow=n.experiments, ncol=3)
IntervalScore05 <- matrix(NA, nrow=n.experiments, ncol=3)
Coverage01 <- matrix(NA, nrow=n.experiments, ncol=3)
IntervalScore01 <- matrix(NA, nrow=n.experiments, ncol=3)
Coverage1 <- matrix(NA, nrow=n.experiments, ncol=3)
IntervalScore1 <- matrix(NA, nrow=n.experiments, ncol=3)
PosteriorMean <- matrix(NA, nrow=n.experiments, ncol=3)
QuantileScoreAverage50 <- rep(NA, n.experiments)
QuantileScoreAverage100 <- rep(NA, n.experiments)

for (iter in 1:n.experiments){
  file=sprintf("%s_progress_%d.RData", experiment.name, iter)
  
  load(file=file)
  
  ## 1. for xi
  alpha <- 0.01
  Intvl <- quantile(out.obj$xi.trace[which], probs = c(alpha/2, 1-alpha/2), names=FALSE)
  Coverage01[iter,1] <- xi<=Intvl[2] & xi>=Intvl[1]
  IntervalScore01[iter,1] <- (Intvl[2]-Intvl[1]) + 
                   2*(Intvl[1]-xi)*as.numeric(xi<Intvl[1])/alpha +
                   2*(xi-Intvl[2])*as.numeric(xi>Intvl[2])/alpha
  
  alpha <- 0.05
  Intvl <- quantile(out.obj$xi.trace[which], probs = c(alpha/2, 1-alpha/2), names=FALSE)
  Coverage05[iter,1] <- xi<=Intvl[2] & xi>=Intvl[1]
  IntervalScore05[iter,1] <- (Intvl[2]-Intvl[1]) + 
    2*(Intvl[1]-xi)*as.numeric(xi<Intvl[1])/alpha +
    2*(xi-Intvl[2])*as.numeric(xi>Intvl[2])/alpha
  
  alpha <- 0.1
  Intvl <- quantile(out.obj$xi.trace[which], probs = c(alpha/2, 1-alpha/2), names=FALSE)
  Coverage1[iter,1] <- xi<=Intvl[2] & xi>=Intvl[1]
  IntervalScore1[iter,1] <- (Intvl[2]-Intvl[1]) + 
    2*(Intvl[1]-xi)*as.numeric(xi<Intvl[1])/alpha +
    2*(xi-Intvl[2])*as.numeric(xi>Intvl[2])/alpha
  
  PosteriorMean[iter,1] <- mean(out.obj$xi.trace[which])
  
  
  
  ## For quantile scores
  for (i in 1:length(which)){
    mu_tmp <- out.obj$mu.trace[which[i]]
    xi_tmp <- out.obj$xi.trace[which[i]]
    tau_tmp <- out.obj$tau.trace[which[i]]
    returnLevels_50[i]  <- qevd(0.98, loc=mu_tmp, shape=xi_tmp, scale=tau_tmp)
    returnLevels_100[i] <- qevd(0.99, loc=mu_tmp, shape=xi_tmp, scale=tau_tmp)
    QuantileScore50[i] <- (L50-returnLevels_50[i])*(as.numeric(L50<=returnLevels_50[i])-0.98)
    QuantileScore100[i] <- (L100-returnLevels_100[i])*(as.numeric(L100<=returnLevels_100[i])-0.99)
  }
  
  
  ## 2. for 50 year
  alpha <- 0.01
  Intvl <- quantile(returnLevels_50, probs = c(alpha/2, 1-alpha/2), names=FALSE)
  Coverage01[iter,2] <- L50<=Intvl[2] & L50>=Intvl[1]
  IntervalScore01[iter,2] <- (Intvl[2]-Intvl[1]) + 
    2*(Intvl[1]-L50)*as.numeric(L50<Intvl[1])/alpha +
    2*(L50-Intvl[2])*as.numeric(L50>Intvl[2])/alpha
  
  alpha <- 0.05
  Intvl <- quantile(returnLevels_50, probs = c(alpha/2, 1-alpha/2), names=FALSE)
  Coverage05[iter,2] <- L50<=Intvl[2] & L50>=Intvl[1]
  IntervalScore05[iter,2] <- (Intvl[2]-Intvl[1]) + 
    2*(Intvl[1]-L50)*as.numeric(L50<Intvl[1])/alpha +
    2*(L50-Intvl[2])*as.numeric(L50>Intvl[2])/alpha
  
  alpha <- 0.1
  Intvl <- quantile(returnLevels_50, probs = c(alpha/2, 1-alpha/2), names=FALSE)
  Coverage1[iter,2] <- L50<=Intvl[2] & L50>=Intvl[1]
  IntervalScore1[iter,2] <- (Intvl[2]-Intvl[1]) + 
    2*(Intvl[1]-L50)*as.numeric(L50<Intvl[1])/alpha +
    2*(L50-Intvl[2])*as.numeric(L50>Intvl[2])/alpha
  
  PosteriorMean[iter,2] <- mean(returnLevels_50)
  QuantileScoreAverage50[iter] <- mean(QuantileScore50)
  
  ## 3. for 100 year
  alpha <- 0.01
  Intvl <- quantile(returnLevels_100, probs = c(alpha/2, 1-alpha/2), names=FALSE)
  Coverage01[iter,3] <- L100<=Intvl[2] & L100>=Intvl[1]
  IntervalScore01[iter,3] <- (Intvl[2]-Intvl[1]) + 
    2*(Intvl[1]-L100)*as.numeric(L100<Intvl[1])/alpha +
    2*(L100-Intvl[2])*as.numeric(L100>Intvl[2])/alpha
  
  alpha <- 0.05
  Intvl <- quantile(returnLevels_100, probs = c(alpha/2, 1-alpha/2), names=FALSE)
  Coverage05[iter,3] <- L100<=Intvl[2] & L100>=Intvl[1]
  IntervalScore05[iter,3] <- (Intvl[2]-Intvl[1]) + 
    2*(Intvl[1]-L100)*as.numeric(L100<Intvl[1])/alpha +
    2*(L100-Intvl[2])*as.numeric(L100>Intvl[2])/alpha
  
  alpha <- 0.1
  Intvl <- quantile(returnLevels_100, probs = c(alpha/2, 1-alpha/2), names=FALSE)
  Coverage1[iter,3] <- L100<=Intvl[2] & L100>=Intvl[1]
  IntervalScore1[iter,3] <- (Intvl[2]-Intvl[1]) + 
    2*(Intvl[1]-L100)*as.numeric(L100<Intvl[1])/alpha +
    2*(L100-Intvl[2])*as.numeric(L100>Intvl[2])/alpha
  
  PosteriorMean[iter,3] <- mean(returnLevels_100)
  QuantileScoreAverage100[iter] <- mean(QuantileScore100)
  
  rm(out.obj);rm(state)
}


FILE <- sprintf("%s_ProperScores.RData", experiment.name)
save(Coverage05, IntervalScore05, Coverage01, IntervalScore01, Coverage1, IntervalScore1, 
     PosteriorMean, QuantileScoreAverage50, QuantileScoreAverage100, file = FILE)

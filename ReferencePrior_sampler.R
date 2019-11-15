###################################################################################
## Main sampler

## Y ...................................................  observations on GEV scale
## mu .......................................................... location parameter
## xi ............................................................. shape parameter
## tau ............................................................ scale parameter
## xi_prior_type ................... One of: "Beta", "h11", "h22_1", "h22_2", "MDI"
## initial.values .................. a list: mu, xi, tau
## n.updates .................................................... number of updates
## thin ............................................. number of runs in each update
## experiment.name
## echo.interval ......................... echo process every echo.interval updates
## sigma.m
## prop.Sigma
## true.params ..................... a list:mu, xi, tau


GEV.sampler.02 <- function(Y, iter.num = 1,   
                              xi_prior_type='Beta',
                              initial.values,
                              n.updates, thin=10,
                              experiment.name="Huser-wadsworth",
                              echo.interval=50,
                              sigma.m=NULL, prop.Sigma=NULL, 
                              true.params=NULL, sd.ratio=NULL, lower.prob.lim=0.5) {
  
  #library(doParallel)
  #library(foreach)
  # save.bit <- TRUE
  
  # Constants to control how many Metropolis steps get executed at each
  # iteration of the main loop
  n.metr.updates.mu <- 4
  n.metr.updates.xi <- 4
  n.metr.updates.tau <- 4
  
  
  # Constants to control adaptation of the Metropolis sampler
  c.0 <- 10
  c.1 <- 0.8
  k <- 3  # the iteration offset
  metr.opt.1d <- 0.41

  
  # A small number
  eps <- 1e-06
  
  # Bookkeeping
  n <- length(Y)
  
  # Load initial values
  mu <- initial.values$mu
  xi <- initial.values$xi
  tau <- initial.values$tau
  
  

  # Initialize trace objects
 
  mu.trace <- rep(NA, n.updates)
  xi.trace <- rep(NA, n.updates)
  tau.trace <- rep(NA, n.updates)
  mu.trace[1] <- mu
  xi.trace[1] <- xi
  tau.trace[1] <- tau
  
 
  
  # For tuning Metropolis updates of theta
  if (is.null(sigma.m$mu)) sigma.m$mu <- 2.4^2
  if (is.null(sigma.m$xi)) sigma.m$xi <- 2.4^2
  if (is.null(sigma.m$tau)) sigma.m$tau <- 2.4^2
 
  
  r.hat.mu <- NA
  r.hat.xi <- NA
  r.hat.tau <- NA
  
  # Choose prior type for xi
  if(xi_prior_type=="Beta")  xi_prior <- beta_prior
  if(xi_prior_type=="h11")  xi_prior <- h11
  if(xi_prior_type=="h22_1")  xi_prior <- h22_1
  if(xi_prior_type=="h22_2")  xi_prior <- h22_2
  if(xi_prior_type=="MDI")  xi_prior <- MDI
  
  
  
  for (i in 2:n.updates) {
    
    ################################################################
    ## Update Metropolis adaptation parameters
    ################################################################
    gamma1 <- c.0 / (i + k)^(c.1)
    gamma2 <- 1 / (i + k)^(c.1)

    for (j in 1:thin) {
      
        
      ################################################################
      ## Update mu
      ################################################################
      metr.out.mu <- static.metr(starting.theta = mu,
                                  likelihood.fn = Lik2, prior.fn = mu_prior,
                                  n.updates = n.metr.updates.mu, prop.Sigma = 1, sigma.m=sigma.m$mu, verbose=FALSE,
                                  Y=Y, xi=xi, tau=tau)
      r.hat.mu <- metr.out.mu$acc.prob
      mu <- metr.out.mu$trace[n.metr.updates.mu]
      sigma.m$mu <- exp(log(sigma.m$mu) + gamma1*(r.hat.mu - metr.opt.1d))
      
      
      ################################################################
      ## Update xi
      ################################################################
      metr.out.xi <- static.metr(starting.theta = xi,
                                 likelihood.fn = Lik2, prior.fn = xi_prior,
                                 n.updates = n.metr.updates.xi, prop.Sigma = 1, sigma.m=sigma.m$xi, verbose=FALSE,
                                 Y=Y, mu=mu, tau=tau)
      r.hat.xi <- metr.out.xi$acc.prob
      xi <- metr.out.xi$trace[n.metr.updates.xi]
      sigma.m$xi<- exp(log(sigma.m$xi) + gamma1*(r.hat.xi - metr.opt.1d))
      
      
      ################################################################
      ## Update tau
      ################################################################
      metr.out.tau <- static.metr(starting.theta = tau,
                                 likelihood.fn = Lik2, prior.fn = tau_prior,
                                 n.updates = n.metr.updates.tau, prop.Sigma = 1, sigma.m=sigma.m$tau, verbose=FALSE,
                                 Y=Y, mu=mu, xi=xi)
      r.hat.tau <- metr.out.tau$acc.prob
      tau <- metr.out.tau$trace[n.metr.updates.tau]
      sigma.m$tau<- exp(log(sigma.m$tau) + gamma1*(r.hat.tau - metr.opt.1d))
    }
    
    
    # ------------------------- Fill in trace objects ---------------------------
    mu.trace[i] <- mu
    xi.trace[i] <- xi
    tau.trace[i] <- tau
    
      
    
    # ----------------------------- Echo progress --------------------------------
    if ((i %% echo.interval) == 0) {
      cat("Done with", i, "updates,\n")
      pdf(file=sprintf("%s_progress_%d.pdf", experiment.name, iter.num))
      par(mfrow=c(3,2))
      plot(mu.trace, type="l", ylab=expression(mu))
      if (!is.null(true.params)) abline(h=true.params$mu, lty=2, col=2, lwd=3)
      plot(xi.trace, type="l", ylab=expression(xi))
      if (!is.null(true.params)) abline(h=true.params$xi, lty=2, col=2, lwd=3)
      plot(tau.trace, type="l", ylab=expression(tau))
      if (!is.null(true.params)) abline(h=true.params$tau, lty=2, col=2, lwd=3)
      dev.off()
      
      state <- list(Y=Y, mu=mu, xi=xi, tau=tau,
                    i=i, sigma.m=sigma.m, prop.Sigma=prop.Sigma)
      out.obj <- list(Y=Y, i=i,
                      mu.trace=mu.trace,
                      xi.trace=xi.trace,
                      tau.trace=tau.trace)
      save(state, out.obj, file=sprintf("%s_progress_%d.RData", experiment.name, iter.num))
      # save.bit <- !save.bit
    }
  }
  
  return(out.obj)
  
}

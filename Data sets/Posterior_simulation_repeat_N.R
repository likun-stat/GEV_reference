source("./generic_samplers.R")
source("./ReferencePrior_utils.R")
source("./ReferencePrior_sampler.R")

library(extRemes)
library(energy)

set.seed(123)

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------
true_Par <- c(10, 0.2, 2)   # mu, xi, tau
xi_prior_type <- "h11"

Tau <- seq(0.1, 10, length.out = 35)
Mu  <- seq(-20, 20, length.out = 35)
Xi  <- seq(-0.5, 0.5, length.out = 35)

Par <- expand.grid(tau = Tau, mu = Mu, xi = Xi)
Par <- cbind(mu = Par$mu, xi = Par$xi, tau = Par$tau)

n.updates <- 10000
thin <- 10
echo.interval <- 1000
choose_pos_num <- 500

n.rep <- 100
N.vec <- c(50, 100, 500, 1000)

# ------------------------------------------------------------
# Helper: posterior mode initialization on grid
# ------------------------------------------------------------
get_init_vals <- function(Current_data, Par, xi_prior_type) {
  lik_Par <- rep(NA_real_, nrow(Par))
  
  for (i in seq_len(nrow(Par))) {
    tmp <- Lik(Par[i, ], Current_data)
    xi_tmp  <- Par[i, 2]
    tau_tmp <- Par[i, 3]
    
    if (xi_prior_type == "Beta")  lik_Par[i] <- tmp + beta_prior(xi = xi_tmp) - log(tau_tmp)
    if (xi_prior_type == "h11")   lik_Par[i] <- tmp + h11(xi_tmp) - log(tau_tmp)
    if (xi_prior_type == "h22_1") lik_Par[i] <- -tmp + h22_1(xi_tmp) - log(tau_tmp)
    if (xi_prior_type == "h22_2") lik_Par[i] <- -tmp + h22_2(xi_tmp) - log(tau_tmp)
    if (xi_prior_type == "MDI")   lik_Par[i] <- tmp + MDI(xi = xi_tmp) - log(tau_tmp)
  }
  
  idx <- which.max(lik_Par)
  
  list(
    mu  = Par[idx, 1],
    xi  = Par[idx, 2],
    tau = Par[idx, 3]
  )
}

# ------------------------------------------------------------
# One replicate
# Returns 4 p-values:
#   p_mv  = multivariate normality test
#   p_mu  = univariate test for mu
#   p_xi  = univariate test for xi
#   p_tau = univariate test for tau
# ------------------------------------------------------------
run_one_rep <- function(N, rep_id,
                        true_Par,
                        Par,
                        xi_prior_type,
                        n.updates,
                        thin,
                        echo.interval,
                        choose_pos_num) {
  
  Current_data <- extRemes::revd(
    N,
    loc   = true_Par[1],
    scale = true_Par[3],
    shape = true_Par[2]
  )
  
  init_vals <- get_init_vals(Current_data, Par, xi_prior_type)
  
  iter_id <- paste0("N", N, "_rep", rep_id)
  
  res <- GEV.sampler.02(
    Y = Current_data,
    iter.num = iter_id,
    xi_prior_type = xi_prior_type,
    initial.values = init_vals,
    n.updates = n.updates,
    thin = thin,
    experiment.name = "station",
    echo.interval = echo.interval,
    true.params = list(mu = true_Par[1], xi = true_Par[2], tau = true_Par[3]),
    sd.ratio = NULL,
    lower.prob.lim = 0.5
  )
  
  mu_vec  <- tail(res$mu.trace[seq(1, length(res$tau.trace), by = 10)],  
                  choose_pos_num)
  xi_vec  <- tail(res$xi.trace[seq(1, length(res$tau.trace), by = 10)],
                  choose_pos_num)
  tau_vec <- tail(res$tau.trace[seq(1, length(res$tau.trace), by = 10)],
                  choose_pos_num)
  
  pos_data <- data.frame(mu = mu_vec, xi = xi_vec, tau = tau_vec)
  pos_data <- apply(pos_data, 2, function(x) (x-mean(x))/sd(x))
  # energy::mvnorm.etest returns an htest-like object with $p.value
  p_mv  <- mvnorm.test(pos_data,      R = 1000)$p.value
  p_mu  <- normal.test(pos_data[, 1], R = 1000)$p.value
  p_xi  <- normal.test(pos_data[, 2], R = 1000)$p.value
  p_tau <- normal.test(pos_data[, 3], R = 1000)$p.value
  
  c(p_mv = p_mv, p_mu = p_mu, p_xi = p_xi, p_tau = p_tau)
}

# ------------------------------------------------------------
# Repeat for one N
# Saves a 1000 x 4 data frame
# ------------------------------------------------------------
run_experiment_for_N <- function(N, n.rep = 1000) {
  out <- matrix(NA_real_, nrow = n.rep, ncol = 4)
  colnames(out) <- c("p_mv", "p_mu", "p_xi", "p_tau")
  
  for (r in seq_len(n.rep)) {
    cat("Running N =", N, ", replicate", r, "of", n.rep, "\n")
    
    tmp <- tryCatch(
      run_one_rep(
        N = N,
        rep_id = r,
        true_Par = true_Par,
        Par = Par,
        xi_prior_type = xi_prior_type,
        n.updates = n.updates,
        thin = thin,
        echo.interval = echo.interval,
        choose_pos_num = choose_pos_num
      ),
      error = function(e) {
        message("Error at N = ", N, ", replicate = ", r, ": ", e$message)
        rep(NA_real_, 4)
      }
    )
    
    out[r, ] <- tmp
  }
  
  out_df <- as.data.frame(out)
  out_df$replicate <- seq_len(n.rep)
  out_df <- out_df[, c("replicate", "p_mv", "p_mu", "p_xi", "p_tau")]
  
  return(out_df)
}

# ------------------------------------------------------------
# Main loop over N = 50, 100, 500, 1000
# ------------------------------------------------------------
all_results <- vector("list", length(N.vec))
names(all_results) <- paste0("N", N.vec)

for (j in seq_along(N.vec)) {
  N <- N.vec[j]
  cat("\n==============================\n")
  cat("Starting experiments for N =", N, "\n")
  cat("==============================\n")
  
  res_N <- run_experiment_for_N(N, n.rep = n.rep)
  all_results[[j]] <- res_N
  
  # Save each 1000 x 4 (+ replicate column) data frame
  write.csv(res_N, file = paste0("pvalues_N", N, "_pos.csv"), row.names = FALSE)
  saveRDS(res_N, file = paste0("pvalues_N", N, "_pos.rds"))
}

# Save all results together as a named list
saveRDS(all_results, file = "all_pvalues_by_N_pos.rds")

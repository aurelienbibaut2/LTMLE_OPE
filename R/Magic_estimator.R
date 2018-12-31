library(MASS)
library(quadprog)
library(boot)
source('Estimators.R')

# One helper function
distance_to_interval <- function(interval, x){
  lb <- interval[1]; ub <- interval[2]
  if(x >= lb && x <= ub){
    return(0)
  }else{
    return(min(abs(x-lb), abs(x-ub)))
  }
}

# Fast bootstrap of g^(horizon)(D) based on the observation that WDR estimator is the mean of 
# trajectory-specific statistics (called xi below), so it's just a matter of resampling these statistics rather 
# than resampling the entire dataset and recomputing the WDR estimator each time
bootstrap_WDR <- function(D, Q_hat, V_hat, gamma, n_bootstrap=1000, alpha=0.1){
  n <- dim(D)[1]
  horizon <- dim(D)[2]
  D_star <- matrix(NA, nrow=n, ncol=horizon)
  # Compute mean rho_t to be used as denominator in the stabilized weights
  w_t <- apply(D[, , 'rho_t'], 2, mean)
  
  t <- horizon
  D_star[, t] <- gamma^t * D[, t, 'rho_t'] / w_t[t] * (D[, t, 'r'] + 
                                                         - apply(D[, t, ], 1, function(x) Q_hat[t, x['s'], x['a']]))
  for(t in (horizon-1):1){
    D_star[, t] <- gamma^t * D[, t, 'rho_t'] / w_t[t] * (D[, t, 'r'] + gamma * V_hat[t+1, D[, t+1, 's']]
                                                         - apply(D[, t, ], 1, function(x) Q_hat[t, x['s'], x['a']]))
  }
  xi <- V_hat[1, D[, 1, 's']] + apply(D_star, 1, sum)
  bootstrapped_estimates <- boot(data = xi, statistic = function(data, indices) mean(data[indices]), R=n_bootstrap)$t
  quantile(bootstrapped_estimates, probs = c(alpha / 2, 1 - alpha / 2))
}

# The MAGIC estimator from Thomas and Brunskill 2016
MAGIC_estimator <- function(D, Q_hat, V_hat, gamma, horizon, n_bootstrap=1000, force_PD=T){
  # J: set of indices j
  if(horizon >= 10){
    J <- (floor(seq(1, horizon, length.out=10))) + 1
  }else{
    J <- (1:horizon) + 1
  }
  # Get g^(j)s
  WDR_results <- WDR_estimator_TB(D, Q_hat=Q_hat, V_hat=V_hat, gamma=gamma, j=horizon, compute_covariance=T)
  g_js <- WDR_results$g_js
  
  # Get bias by bootstrapping g^(horizon)
  bootstrap_CI <- bootstrap_WDR(D, Q_hat=Q_hat, V_hat=V_hat, gamma=gamma, n_bootstrap=n_bootstrap, alpha=0.1)
  b_n <- sapply(g_js, Vectorize(function(g_j) distance_to_interval(bootstrap_CI, g_j)) )
  
  # Solving x^\top D x under the constraint that A^\top x >= b0. 
  # First row of A is actually an equality constraint. This is specified by setting meq=1 in solve.QP
  n <- dim(D)[1]
  Dmat <- WDR_results$Omega_n[J, J]/n + b_n[J] %*% t(b_n[J])
  if(force_PD)
    Dmat <- Matrix::nearPD(Dmat, eig.tol=1e-10)$mat
  Amat <- t(rbind(rep(1, length(J)), diag(length(J))))
  dvec <- rep(0, length(J))
  b0 <- c(1, rep(0, length(J)))
  x_star <- solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=b0, meq=1)$solution
  
  # Compute the MAGIC estimate as the weighted sum of the g^(j)'s, that is x_star^\top b_n[2:horizon]
  estimate <- t(x_star) %*% g_js[J]
  list(estimate=estimate, x_star=x_star, g_js=g_js, b_n=b_n, Omega_n=WDR_results$Omega_n, bootstrap_CI=bootstrap_CI, J=J)
}
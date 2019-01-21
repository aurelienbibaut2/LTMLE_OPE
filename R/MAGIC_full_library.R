library(boot)
source('Estimators.R')
source('utils.R')
source('partial_LTMLE.R')
source('Magic_estimator.R')

# MAGIC-WDR selector + partial-softened-LTMLE -----------------------------
MAGIC_full_library_estimator <- function(D, Q_hat, V_hat, evaluation_action_matrix, gamma, horizon, n_bootstrap=1000, force_PD=T){
  n <- dim(D)[1]
  
  # J: set of indices j
  if(horizon >= 10){
    J <- (floor(seq(1, horizon, length.out=10))) + 1
  }else{
    J <- (1:horizon) + 1
  }
  # Get g^(j)s
  WDR_results <- WDR_estimator_TB(D, Q_hat=Q_hat, V_hat=V_hat, gamma=gamma, j=horizon, compute_covariance=T)
  g_js_WDR <- WDR_results$g_js[J]
  centered_xi_WDR <- WDR_results$centered_xi[, J]
  
  # Compute the corresponding sequence of partial softened LTMLE
  n_alphas <- 10
  alphas_bis <- c(rep(0, n_ids/2), seq(0, 1, length.out = n_alphas/2))
  lambdas_bis <- c(rev(seq(0, 5e-5, length.out = n_ids/2)), rep(0, n_alphas/2))
  
  js <- (ceiling(seq(1, horizon, length.out=n_alphas)))
  alphas <- seq(0, 1, length.out = n_alphas)
  lambdas <- seq(0, 1e-4, length.out = n_alphas)
  LTMLE_results <- sapply(1:n_alphas, function(j) partial_LTMLE_estimator(D, Q_hat=Q_hat, V_hat=V_hat, gamma=gamma, 
                                                                       evaluation_action_matrix = evaluation_action_matrix, 
                                                                       alpha=alphas_bis[j], j=js[j],
                                                                       lambda=lambdas_bis[j]))
  xi_LTMLE <- sapply(1:n_alphas, function(j) evaluate_stabilized_EIC(D, epsilons = LTMLE_results[, j]$epsilons, 
                                                                      Q_hat, V_hat, evaluation_action_matrix, gamma, j=js[j]))
  xi_LTMLE_bar <- apply(xi_LTMLE, 2, mean)
  centered_xi_LTMLE <- xi_LTMLE - rep(1, n) %*% t(xi_LTMLE_bar)
  g_js_LTMLE <- unlist(LTMLE_results[1, ])
  
  # Compute global Omega_n
  centered_xi <- cbind(
    #centered_xi_WDR, 
    centered_xi_LTMLE)
  Omega_n <- t(centered_xi) %*% centered_xi / n
  
  
  # Get bias by bootstrapping g^(horizon)
  g_js <- c(
    #g_js_WDR, 
    g_js_LTMLE)
  
  bootstrap_CI <- bootstrap_WDR(D, Q_hat=Q_hat, V_hat=V_hat, gamma=gamma, n_bootstrap=n_bootstrap, alpha=1)
  b_n <- sapply(g_js, Vectorize(function(g_j) distance_to_interval(bootstrap_CI, g_j)) )
  
  # Solving x^\top D x under the constraint that A^\top x >= b0. 
  # First row of A is actually an equality constraint. This is specified by setting meq=1 in solve.QP
  
  Dmat <- Omega_n/n + b_n %*% t(b_n)
  if(force_PD)
    Dmat <- Matrix::nearPD(Dmat, eig.tol=1e-10)$mat
  Amat <- t(rbind(rep(1, length(b_n)), diag(length(b_n))))
  dvec <- rep(0, length(b_n))
  b0 <- c(1, rep(0, length(b_n)))
  x_star <- solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=b0, meq=1)$solution
  
  # Compute the MAGIC estimate as the weighted sum of the g^(j)'s, that is x_star^\top b_n[2:horizon]
  # estimate <- t(x_star) %*% g_js[J]
  estimate <- t(x_star) %*% g_js
  list(estimate=estimate, x_star=x_star, g_js=g_js, g_js_LTMLE=g_js_LTMLE, 
       b_n=b_n, Omega_n=Omega_n, bootstrap_CI=bootstrap_CI, J=J)
}


# MAGIC full library debugging experiments -------------------------------------
# source('MDP_modelWin.R')
# horizon <- 5; gamma <- 1; n_states <- 3; n_actions <- 2
# V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
#                                   transition_based_rewards,
#                                   evaluation_action_matrix, horizon, gamma = gamma)
# V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0
# Q_hat <- array(dim=dim(Q0)); V_hat <- array(dim=dim(V0))
# Q_hat <- Q0; V_hat <- V0
# b <- 0 * rnorm(1)
# Delta_t <- 0
# 
# D <- generate_discrete_MDP_dataset(100, 1, state_transition_matrix,
#                                    behavior_action_matrix,
#                                    transition_based_rewards,
#                                    horizon)
# 
# 
# # debug(MAGIC_LTMLE_estimator)
# print(MAGIC_full_library_estimator(D, Q_hat, V_hat,
#                         evaluation_action_matrix, 1, horizon=horizon
#                         ))
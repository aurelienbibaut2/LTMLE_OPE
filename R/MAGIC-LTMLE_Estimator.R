library(boot)
source('Estimators.R')
source('utils.R')
source('partial_LTMLE.R')


# MAGIC-WDR selector + partial-softened-LTMLE -----------------------------
MAGIC_LTMLE_estimator <- function(D, Q_hat, V_hat, evaluation_action_matrix, gamma, horizon, n_bootstrap=1000, force_PD=T){
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
  
  # Compute the corresponding sequence of partial softened LTMLE
  n_alphas <- 10
  js <- (ceiling(seq(1, horizon, length.out=n_alphas))) 
  alphas <- seq(0, 1, length.out = n_alphas)
  
  g_js_LTMLE <- sapply(1:n_alphas, function(j) partial_LTMLE_estimator(D, Q_hat=Q_hat, V_hat=V_hat, gamma=gamma, 
                                                                 evaluation_action_matrix = evaluation_action_matrix, alpha=alphas[j], j=js[j])$estimate)
  
  # Compute the MAGIC estimate as the weighted sum of the g^(j)'s, that is x_star^\top b_n[2:horizon]
  # estimate <- t(x_star) %*% g_js[J]
  estimate <- t(x_star) %*% g_js_LTMLE
  list(estimate=estimate, x_star=x_star, g_js=g_js, g_js_LTMLE=g_js_LTMLE, 
       b_n=b_n, Omega_n=WDR_results$Omega_n, bootstrap_CI=bootstrap_CI, J=J)
}



# Estimator ---------------------------------------------------------------
MAGIC_bootstrap_LTMLE <- function(D, Q_hat, V_hat, gamma, evaluation_action_matrix, force_PD=T){
  n <- dim(D)[1]
  horizon <- dim(D)[2]
  
  # # Get g^(j)s
  n_alphas <- 10
  js <- (ceiling(seq(1, horizon, length.out=n_alphas))) 
  alphas <- seq(0, 1, length.out = n_alphas)
  R1 <- horizon * 2
  R2 <- horizon * 2
  g_js <- sapply(1:n_alphas, function(j) partial_LTMLE_estimator(D, Q_hat=Q_hat, V_hat=V_hat, gamma=gamma, 
                                                         evaluation_action_matrix = evaluation_action_matrix, alpha=alphas[j], j=js[j])$estimate)
  
  # # Get bias by bootstrapping g^(horizon)
  bootstrap_CI <- quantile(boot(data=D, statistic=function(data, indices) partial_LTMLE_estimator(D[indices, , ], Q_hat=Q_hat, V_hat=V_hat, gamma=gamma, 
                                                                                          evaluation_action_matrix = evaluation_action_matrix, 
                                                                                          alpha=1, j=horizon)$estimate,
                                R = R1)$t, probs = c(0.1 / 2, 1 - 0.1 / 2))
  b_n <- sapply(g_js, Vectorize(function(g_j) distance_to_interval(bootstrap_CI, g_j)) )
  
  # Get covariance matrix by bootstrapping all g_js
  X <- boot(data=D, 
            statistic=function(data, indices) sapply(1:n_alphas, 
                                                     function(j) partial_LTMLE_estimator(D[indices, , ], 
                                                                                     Q_hat=Q_hat, V_hat=V_hat, gamma=gamma, 
                                                                                     evaluation_action_matrix = evaluation_action_matrix, 
                                                                                     alpha=alphas[j], j=js[j])$estimate),
            R = R2)$t
  Omega_n <- t((X - rep(1, horizon * 2) %*% t(apply(X, 2, mean)))) %*% (X - rep(1, horizon * 2) %*% t(apply(X, 2, mean))) / R2
  
  # Define and solve QP
  Dmat <- Omega_n + b_n %*% t(b_n)
  if(force_PD)
    Dmat <- Matrix::nearPD(Dmat, eig.tol=1e-10)$mat
  Amat <- t(rbind(rep(1, n_alphas), diag(n_alphas)) )
  dvec <- rep(0, n_alphas)
  b0 <- c(1, rep(0, n_alphas))
  x_star <- solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=b0, meq=1)$solution
  
  # Compute the MAGIC estimate as the weighted sum of the g^(j)'s, that is x_star^\top b_n[2:horizon]
  estimate <- t(x_star) %*% g_js
  # Output
  list(estimate=estimate, x_star=x_star, g_js=g_js, b_n=b_n, Omega_n=Omega_n, bootstrap_CI=bootstrap_CI)
}

# # Debugging simulations ---------------------------------------------------
# Parameters of the environment
# n_states <- 3; n_actions <- 2
# # Set DGP parameters
# horizon <- 10; n <- 1e3; gamma <- 1
# V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
#                                   transition_based_rewards,
#                                   evaluation_action_matrix, horizon, gamma)
# # # # Compute truth
# V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0
# b <- 1e-5*rnorm(1)
# # b <- 0
# Q_hat <- Q0 +  b
# V_hat <-  V0 +  b
# 
# D <- generate_discrete_MDP_dataset(n, 1, state_transition_matrix,
#                                    behavior_action_matrix,
#                                    transition_based_rewards,
#                                    horizon)
# # res <- MAGIC_bootstrap_LTMLE(D, Q_hat, V_hat, gamma, evaluation_action_matrix, force_PD=T)
# res <- MAGIC_LTMLE_estimator(D, Q_hat, V_hat, evaluation_action_matrix, gamma, horizon, n_bootstrap=1000, force_PD=T)
# 
# 
# 

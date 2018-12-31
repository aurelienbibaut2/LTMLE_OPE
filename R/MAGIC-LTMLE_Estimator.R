library(boot)
source('Magic_estimator.R')
# Two helper functions
# Two helper functions
# Expit and logit
expit <- function(x) { 1 / (1 + exp(-x)) } 
logit_scalar <- Vectorize(function(x) { 
  if(x > 0 & x < 1){
    return(log(x / (1 - x)))
  }else if(x <= 0){
    return(-Inf)
  }else{
    return(Inf)
  }
})
logit <- function(x){
  if(is.null(dim(x))){
    return(logit_scalar(x))
  }else{
    return(array(logit_scalar(x), dim=dim(x)))
  }
}

# Soften weights
soften <- function(a, alpha){
  a^alpha / sum(a^alpha) 
}

bound <- function(x){
  tol <- 1e-4
  x * as.numeric(x >= tol & x <= 1-tol) + (1-tol) * as.numeric(x > 1-tol) + tol * as.numeric(x < tol)
}

# LTMLE
partial_LTMLE <-  function(D, Q_hat, V_hat, evaluation_action_matrix, gamma, alpha, j){
  # Get dataset dimensions
  n <- dim(D)[1]
  horizon <- dim(D)[2]
  D_star <- matrix(0, nrow=n, ncol=j)
  Delta_t <- 0
  t <- horizon
  while(t > j){
    Delta_t <- 1 + gamma * Delta_t
    t <- t - 1
  }
  
  epsilons <- c()
  if(t == horizon){
    V_evaluated <- rep(0, n) #V_{t+1}(S_{t+1})
  }else{
    V_evaluated <- V_hat[t+1, D[, t+1, 's']]
  }
  for(t in j:1){
    Delta_t <- 1 + gamma * Delta_t
    R <- D[, t, 'r'] # R_t
    U_tilde <- bound((R + gamma * V_evaluated + Delta_t) / (2 * Delta_t)) # U_tilde = R_tilde_t + gamma*V_tilde_{t+1}(S_t) in the notations of the write-up
    Q_t_evaluated <- apply(D[ , t, ], 1, function(x) Q_hat[t, x['s'], x['a']]) # Q_t(A_t, S_t)
    Q_tilde_t_evaluated <- (Q_t_evaluated + Delta_t) / (2 * Delta_t)
    epsilon <- glm(U_tilde ~ offset(logit(Q_tilde_t_evaluated)) + 1, 
                   family=quasibinomial, weights = soften(D[,t, 'rho_t'], alpha) )$coefficients[1]
    epsilons <- c(epsilons, epsilon)
    # Evaluate Q_tilde(s_t, a) for a_t = 1, a_t = 2
    Q_tilde_t_star <- expit(logit( (Q_hat[t, ,] + Delta_t) / (2 * Delta_t) ) + epsilon) # Q_tilde_t^*
    
    # Then set V_tilde(s_t) = \sum_{a} Q_tilde(s_t, a) pi_a(a|s_t)
    V_tilde <- apply(Q_tilde_t_star * evaluation_action_matrix, 1, sum)
    
    # Compute V = 2 * Delta_t * (V_tilde - 1)
    V <- 2 * Delta_t * (V_tilde - 1/2)
    Q_t_star <- 2 * Delta_t * (Q_tilde_t_star - 1/2)
    
    # Compute efficient influence curve at Q_t_star
    # Note that at that point we haven't updated V_evaluated yet, so it is still V^*(S_{t+1})
    # Also note that the case t=horizon is already taken care of since at t=horizon V_evaluated=0
    D_star[, t] <- gamma^t * D[, t, 'rho_t'] * (D[, t, 'r'] + gamma * V_evaluated
                                                - apply(D[, t, ], 1, function(x) Q_t_star[x['s'], x['a']]))
    # Evaluate V
    V_evaluated <-  V[D[, t, 's']]
    
    
    # V_tilde is gonna be V_{t+1} in the next iteration
  }
  # The average of the last V is the LTML estimator of the value. We return the V vector pre-averaging
  list(estimate=mean(V_evaluated), EIC=apply(D_star, 1, sum))
}

# MAGIC LTMLE. Hacky version that uses the x_star vector corresponding to WDR and uses it to weight the equivalent 
# sequence of truncated LTMLEs
MAGIC_LTMLE_estimator_hacky <- function(D, Q_hat, V_hat, evaluation_action_matrix, gamma, n_bootstrap){
  horizon <- dim(D)[2]
  g_js_LTMLE <- c()
  for(t in 1:horizon){
    g_js_LTMLE <- cbind(g_js_LTMLE, partial_LTMLE(D, Q_hat, V_hat, evaluation_action_matrix, gamma=gamma, alpha=1, j=t)$estimate)
  }
  MAGIC_result <- MAGIC_estimator(D, Q_hat, V_hat, gamma, horizon, n_bootstrap=n_bootstrap)
  t(MAGIC_result$x_star) %*% as.vector(g_js_LTMLE[MAGIC_result$J-1])
}

# MAGIC LTMLE. Proper version
MAGIC_LTMLE_estimator <- function(D, Q_hat, V_hat, evaluation_action_matrix, gamma, n_bootstrap){
  horizon <- dim(D)[2]
  g_js_LTMLE <- c()
  xi <- c()
  # Get estimates and efficient influence curve (EIC) for the sequence of truncated LTMLEs
  for(t in 1:horizon){
    partial_LTMLE_result <- partial_LTMLE(D, Q_hat, V_hat, evaluation_action_matrix, gamma=gamma, alpha=1, j=t)
    g_js_LTMLE <- cbind(g_js_LTMLE, partial_LTMLE_result$estimate)
    xi <- cbind(xi, partial_LTMLE_result$EIC - mean(partial_LTMLE_result$EIC))
  }
  # Compute covariance matrix
  Omega_n <- t(xi) %*% xi / n
  
  # Get bias by bootstrapping g^(horizon)
  bootstrapped_estimates <- boot(data = D, 
                                 statistic = function(data, indices) LTMLE_estimator(data[indices, ,],
                                                                                     Q_hat, V_hat, evaluation_action_matrix, gamma=gamma, alpha=1), 
                                 R=20)$t
  bootstrap_CI <- quantile(bootstrapped_estimates, probs = c(0.8 / 2, 1 - 0.8 / 2))
  b_n <- sapply(g_js_LTMLE, Vectorize(function(g_j) distance_to_interval(bootstrap_CI, g_j)) )
  
  # Solving x^\top D x under the constraint that A^\top x >= b0. 
  # First row of A is actually an equality constraint. This is specified by setting meq=1 in solve.QP
  Dmat <- Omega_n + b_n %*% t(b_n)
  Amat <- t(rbind(rep(1, horizon), diag(horizon)))
  dvec <- rep(0, horizon)
  b0 <- c(1, rep(0, horizon))
  x_star <- solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=b0, meq=1)$solution
  
  # Compute the MAGIC estimate as the weighted sum of the g^(j)'s, that is x_star^\top b_n[2:horizon]
  estimate <- t(x_star) %*% as.vector(g_js_LTMLE)
  list(estimate=estimate, x_star=x_star, g_js=g_js_LTMLE, b_n=b_n, Omega_n=Omega_n, bootstrap_CI=bootstrap_CI)
}

# Debugging experiments ---------------------------------------------------
# horizon <- 20; n <- 100; gamma <- 0.9
# V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
#                                   transition_based_rewards,
#                                   evaluation_action_matrix, horizon, gamma = gamma)
# V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0
# 
# b <- 1e-1*rnorm(1)
# Q_hat <- Q0 +  b
# V_hat <-  V0 +  b
# 
# D <- generate_discrete_MDP_dataset(n, 1, state_transition_matrix,
#                                    behavior_action_matrix,
#                                    transition_based_rewards,
#                                    horizon)
# 
# MAGIC_LTMLE_estimator_hacky(D, Q_hat, V_hat, evaluation_action_matrix, gamma = gamma, n_bootstrap = 100)
# MAGIC_LTMLE_estimator(D, Q_hat, V_hat, evaluation_action_matrix, gamma = gamma, n_bootstrap = 20)$estimate
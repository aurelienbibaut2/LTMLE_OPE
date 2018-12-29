source('Estimators.R')

# Two helper functions
# Expit and logit
expit <- function(x) { 1 / (1 + exp(-x)) } 
logit <- Vectorize(function(x) { 
  if(x > 0 & x < 1){
    return(log(x / (1 - x)))
  }else if(x <= 0){
    return(-Inf)
  }else{
    return(Inf)
  }
})

# Evaluate EIC, potentially on a split that is not the one one which epsilons have been fitted
evaluate_EIC <- function(D, epsilons, Q_hat, V_hat, evaluation_action_matrix, gamma){
  n <- dim(D)[1]; horizon <- dim(D)[2]
  D_star <- matrix(0, nrow=n, ncol=horizon)
  V_t_evaluated <- rep(0, n)
  Delta_t <- 0
  for(t in horizon:1){
    Delta_t <- 1 + gamma * Delta_t
    # Compute the fluctuated estimator \hat{Q}^*_t from \hat{Q}_t and epsilon_t
    Q_tilde_t_star <- expit( array(logit( (Q_hat[t, ,] + Delta_t) / (2 * Delta_t) ), dim=dim(Q_hat[t, ,]) ) 
                             + epsilons[t])
    V_tilde_t_star <- as.vector(apply(Q_tilde_t_star * evaluation_action_matrix, 1, sum))
    Q_t_star <- 2 * Delta_t * (Q_tilde_t_star - 1/2)
    V_t_star <- 2 * Delta_t * (V_tilde_t_star - 1/2)
    # Compute component t of the EIC
    # Note that at this point we haven't updated V_t_evaluated yet, so it is still V_{t+1}(S_{t+1})
    D_star[, t] <- gamma^t * D[, t, 'rho_t'] * (D[, t, 'r'] + gamma * V_t_evaluated
                                                - apply(D[, t, ], 1, function(x) Q_t_star[x['s'], x['a']]))
    # Evaluate V
    V_t_evaluated <-  V_t_star[D[, t, 's']]
  }
  apply(D_star, 1, sum)
}

# Collaborative TMLE that uses a sequence of decreasingly softened weights
C_LTMLE_softening <- function(D, Q_hat, V_hat, evaluation_action_matrix, gamma, V=3){
  # Split the dataset. Compute sequence of epsilons for each softening coeff. Pick the one that minimizes a cross-validated risk
  # D_split1 <- D[1:floor(0.5 * n), ,]; D_split2 <- D[(floor(0.5 * n) + 1):n, ,]
  n <- dim(D)[1]
  softening_coeffs <- 10^seq(from=-1, to=0, length.out = 10)
  CV_doubly_roubust_risks <- rep(0, length(softening_coeffs))
  epsilons_log <- c()
  for(i in 1:length(softening_coeffs)){
    for(v in 1:V){
      test_ids <- (floor( (v-1) / V * n) + 1) : floor( v / V * n)
      epsilons <- LTMLE_estimator(D[setdiff(1:n, test_ids), ,], Q_hat, V_hat, evaluation_action_matrix, gamma, alpha=softening_coeffs[i])$epsilons
      CV_doubly_roubust_risks[i] <- (CV_doubly_roubust_risks[i] +
                                       var(evaluate_EIC(D[test_ids, ,], epsilons, Q_hat, V_hat, evaluation_action_matrix, gamma))
                                     )
    }
  }
  # plot(softening_coeffs, CV_doubly_roubust_risks)
  # print(CV_doubly_roubust_risks)
  # Output the estimate computed on the full dataset that uses the softening coeff that minimizes the CV risk
  list(estimate=LTMLE_estimator(D, Q_hat, V_hat, evaluation_action_matrix, gamma, 
                  alpha=softening_coeffs[which.min(CV_doubly_roubust_risks)]
                  )$estimate,
       softening_coeff=softening_coeffs[which.min(CV_doubly_roubust_risks)])
}

# Debugging experiments ---------------------------------------------------
# Set DGP parameters
# horizon <- 20; n <- 1e2; gamma <- 0.9
# V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
#                                   transition_based_rewards,
#                                   evaluation_action_matrix, horizon, gamma)
# V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0
# b <- 0.1e-1*rnorm(1)
# Q_hat <- Q0 +  b
# V_hat <-  V0 +  b
# D <- generate_discrete_MDP_dataset(n, 1, state_transition_matrix,
#                                    behavior_action_matrix,
#                                    transition_based_rewards,
#                                    horizon)
# # debug(C_LTMLE_softening)
# print(C_LTMLE_softening(D, Q_hat, V_hat, evaluation_action_matrix, gamma, V=3))
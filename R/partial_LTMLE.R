targeting_risk <- function(U_tilde, Q_tilde_t_evaluated, omega, lambda, epsilon){
  (-mean(omega * ( U_tilde * log(expit(logit(Q_tilde_t_evaluated) + epsilon)) + 
                     ( 1- U_tilde) * log(1 - expit(logit(Q_tilde_t_evaluated) + epsilon) ) )
         )
    + lambda  * abs(epsilon)
  )
}
  
partial_LTMLE_estimator <-  function(D, Q_hat, V_hat, evaluation_action_matrix, gamma=1, alpha, j=horizon, lambda=0){
  # Get dataset dimensions
  n <- dim(D)[1]
  horizon <- dim(D)[2]
  
  epsilons <- rep(0, horizon)
  Delta_t <- 0
  V_evaluated <- rep(0, n) #V_{t+1}(S_{t+1})
  
  for(t in horizon:1){
    R <- D[, t, 'r'] # R_t
    
    #Address issue with GridWorld
    #max(D[, , 'r']) > 1
    Delta_t <- max(D[, , 'r']) + gamma * Delta_t

    U_tilde <- (R + gamma * V_evaluated + Delta_t) / (2 * Delta_t) # U_tilde = R_tilde_t + gamma*V_tilde_{t+1}(S_t) in the notations of the write-up
    Q_t_evaluated <- apply(D[ , t, ], 1, function(x) Q_hat[t, x['s'], x['a']]) # Q_t(A_t, S_t)
    Q_tilde_t_evaluated <- (Q_t_evaluated + Delta_t) / (2 * Delta_t)
    Q_tilde_t_evaluated <- bound(Q_tilde_t_evaluated)
    # Targeting step
    if(t <= j){
      # epsilons[t] <- glm(U_tilde ~ offset(logit(Q_tilde_t_evaluated)) + 1, 
                         # family=quasibinomial, weights = soften(D[,t, 'rho_t'], alpha) )$coefficients[1]
      epsilons[t] <- optimize(function(x) targeting_risk(U_tilde, Q_tilde_t_evaluated, 
                                                         omega=soften(D[,t, 'rho_t'], alpha), 
                                                         lambda, x), lower=-6, upper=6)$minimum
    }
    # Evaluate Q_tilde(s_t, a) for a_t = 1, a_t = 2
    Q_tilde_t_star <- expit( logit( (Q_hat[t, , ] + Delta_t) / (2 * Delta_t)  )
                             + epsilons[t]) # Q_tilde_t^*
    # Then set V_tilde(s_t) = \sum_{a} Q_tilde(s_t, a) pi_a(a|s_t)
    V_tilde <- apply(Q_tilde_t_star * evaluation_action_matrix, 1, sum)
    # Compute V = 2 * Delta_t * (V_tilde - 1)
    V <- 2 * Delta_t * (V_tilde - 1/2)
    # Evaluate V
    V_evaluated <-  V[D[, t, 's']]
    # V_tilde is gonna be V_{t+1} in the next iteration
  }
  # The average of the last V is the LTML estimator of the value
  V_hat_LTMLE <- mean(V_evaluated)
  list(estimate=V_hat_LTMLE, epsilons=epsilons)
}

# Partial LTMLE debugging experiments -------------------------------------
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
# # debug(partial_LTMLE_estimator)
# print(partial_LTMLE_estimator(D, Q_hat, V_hat,
#                         evaluation_action_matrix, 1,
#                         alpha=1, j=horizon, lambda=1e-4))

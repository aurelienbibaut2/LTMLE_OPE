source('MDP_modelWin.R')
source('utils.R')
source('evaluate_EIC.R')

# Estimator definition ----------------------------------------------------
one_step_LTMLE <- function(D, Q_hat, V_hat, evaluation_action_matrix, gamma=1, verbose=F, max_it=20, step_size=1e-1){
  # Get data dimensions
  n <- dim(D)[1]; horizon <- dim(D)[2]

  epsilon <- 0; epsilons <- epsilon
  score_eqs <- mean(evaluate_EIC(D, rep(epsilon, horizon), Q_hat, V_hat, evaluation_action_matrix, gamma))
  
  for(it in 1:max_it){
    # Take a gradient step
    score_eq <- mean(evaluate_EIC(D, rep(epsilon, horizon), Q_hat, V_hat, evaluation_action_matrix, gamma))
    epsilon <- epsilon + step_size * score_eq
    
    score_eqs <- c(score_eqs, score_eq) 
    epsilons <- c(epsilons, epsilon)
  }
  
  Delta_t <- 1
  for(t in horizon:1) Delta_t <- 1 + gamma * Delta_t
  # Compute the estimate
  Q_tilde_t_star <- expit( logit( (Q_hat[1, ,] + Delta_t) / (2 * Delta_t) )
                           + epsilon )
  V_tilde_t_star <- as.vector(apply(Q_tilde_t_star * evaluation_action_matrix, 1, sum))
  V_t_star <- 2 * Delta_t * (V_tilde_t_star - 1/2)
  V_t_star_evaluated <- V_t_star[D[, 1, 's']]
  estimate <- mean(V_t_star_evaluated)
  
  # Communicate a bit
  if(verbose){
    par(mfrow=c(2,1))
    plot(epsilons, ylim = c(-max(abs(epsilons)), max(abs(epsilons))) )
    abline(h=0)
    plot(score_eqs, ylim = c(-max(abs(score_eqs)), max(abs(score_eqs))) )
    print(score_eqs)
    abline(h=0)
  }
  
  list(estimate=estimate, epsilon=epsilon, score_eq=score_eq)
}
# 
# # Debugging simulations ---------------------------------------------------
# # Parameters of the environment
# n_states <- 3; n_actions <- 2
# # Set DGP parameters
# horizon <- 7; n <- 1e3; gamma <- 1
# V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
#                                   transition_based_rewards,
#                                   evaluation_action_matrix, horizon, gamma)
# # Compute truth
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
# print(one_step_LTMLE(D, Q_hat, V_hat, evaluation_action_matrix, gamma=1, verbose=T, max_it=20, step_size=1e-1))
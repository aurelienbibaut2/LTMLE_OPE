source('Estimators.R')
library(boot)
source('MDP_modelWin.R')
source('C_TMLE.R')

# ModelWin parameters
horizon <- 5; gamma <- 1; n_states <- 3; n_actions <- 2
V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
                                  transition_based_rewards,
                                  evaluation_action_matrix, horizon, gamma = gamma)
V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0


Q_hat <- array(dim=dim(Q0)); V_hat <- array(dim=dim(V0))
b <- 5e-2 * rnorm(1)
Delta_t <- 0
for(t in horizon:1){
  Delta_t <- 1 + gamma * Delta_t
  Q_hat[t, ,] <- 2 * Delta_t * (expit(logit( (Q0[t, ,] + Delta_t) / (2*Delta_t) ) + b) - 1/2)
  V_hat[t, ] <- 2 * Delta_t * (expit(logit( (V0[t, ] + Delta_t) / (2 * Delta_t) ) + b) - 1/2)
}

n <- 1e2; alpha <- 0.1

estimates <- c()
for(i in 1:100){
  D <- generate_discrete_MDP_dataset(n, 1, state_transition_matrix,
                                     behavior_action_matrix,
                                     transition_based_rewards,
                                     horizon)
  estimates <- c(estimates, LTMLE_estimator(D, Q_hat, V_hat, evaluation_action_matrix, gamma = 1, alpha=alpha)$estimate)
}

cat('bias: ', mean(estimates - V0[1,1]), '\n')
cat('var: ', var(estimates), '\n')
boostrapped_estimates <- boot(data=D, statistic=function(data, indices) LTMLE_estimator(data[indices, , ], Q_hat, V_hat, evaluation_action_matrix, gamma = 1, alpha=alpha)$estimate, R = 100)$t
cat('bootstrap bias: ', mean(boostrapped_estimates) - estimates[length(estimates)], '\n')
cat('bootstrap variance:', var(boostrapped_estimates), '\n')
epsilons <- LTMLE_estimator(D, Q_hat, V_hat, evaluation_action_matrix, gamma = 1, alpha=alpha)$epsilons
cat('EIC variance:', var(evaluate_EIC(D, epsilons, Q_hat, V_hat, evaluation_action_matrix, gamma, loss_softening_coeff=1)) / n)

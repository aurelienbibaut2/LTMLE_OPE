library(MASS)
library(quadprog)
library(boot)
source('MDP_modelWin.R')
source('Estimators.R')
source('Q_learning_discrete_state_space.R')
source('Magic_estimator.R')

set.seed(1)
# Parameters of the environment
n_states <- 3; n_actions <- 2

# Set DGP parameters
horizon <- 20; n <- 1e3; gamma <- 0.8
V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
                                  transition_based_rewards,
                                  evaluation_action_matrix, horizon, gamma)
# Compute truth
V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0
# Generate data
D <- generate_discrete_MDP_dataset(n, 1, state_transition_matrix,
                                   behavior_action_matrix,
                                   transition_based_rewards,
                                   horizon)
# D_split1 <- D[1:(floor(n/2)), ,]
# D_split2 <- D[(1+floor(n/2)):n, , ]
D_split1 <- D; D_split2 <- D
# Fit Q_function
cat('Start Bellman iterations\n')
transitions_dataset <- make_transitions_dataset(D_split1)
Q_learning_results <- bellman_iterations(transitions_dataset, evaluation_action_matrix, 
                                         gamma=gamma, n_iterations=30, V0=V0[1, ])
cat('Done with Bellman iterations\n')
plot(Q_learning_results$l2_errors)
# Replicate accross time points the Bellman iterations based Q-function
# That would be correct under infinite horizon.
# Here's it's only approximately correct. Better for early time points than late time points.
Q_hat <- array(0, dim=c(horizon, n_states, n_actions)); V_hat <- array(0, dim=c(horizon, n_states))
for(t in 1:horizon){
  Q_hat[t, ,] <-  Q_learning_results$Q_hat
  V_hat[t, ] <-  Q_learning_results$V_hat
}

MAGIC_estimates <- c(); WDR_estimates <- c(); LTMLE_estimates <- c()
for(i in 1:100){
  D <- generate_discrete_MDP_dataset(n, 1, state_transition_matrix,
                                     behavior_action_matrix,
                                     transition_based_rewards,
                                     horizon)
  MAGIC_estimates <- c(MAGIC_estimates, 
                       MAGIC_estimator(D, Q_hat, V_hat, gamma = gamma, horizon = horizon, n_bootstrap = 100)$estimate)
  WDR_estimates <- c(WDR_estimates, 
                     WDR_estimator_TB_old(D, Q_hat, V_hat, gamma = gamma))
  LTMLE_estimates <- c(LTMLE_estimates,
                       LTMLE_estimator(D, Q_hat, V_hat, evaluation_action_matrix, gamma))
}

cat('WDR MSE:', mean( (WDR_estimates - V0[1,1])^2), '\n')
cat('MAGIC MSE:', mean( (MAGIC_estimates - V0[1,1])^2), '\n')
cat('LTMLE MSE:', mean( (LTMLE_estimates - V0[1,1])^2), '\n')
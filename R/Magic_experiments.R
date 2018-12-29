library(MASS)
library(quadprog)
library(boot)
source('MDP_modelWin.R')
source('Estimators.R')
source('Q_learning_discrete_state_space.R')
source('Magic_estimator.R')
source('MAGIC-LTMLE_Estimator.R')
source('C_TMLE.R')
# set.seed(1)
# Parameters of the environment
n_states <- 3; n_actions <- 2

# Set DGP parameters
horizon <- 20; n <- 1e2; gamma <- 0.9
V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
                                  transition_based_rewards,
                                  evaluation_action_matrix, horizon, gamma)
# Compute truth
V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0

# Q-learning of the initial estimator -------------------------------------
# Generate data
# D <- generate_discrete_MDP_dataset(n, 1, state_transition_matrix,
#                                    behavior_action_matrix,
#                                    transition_based_rewards,
#                                    horizon)
# D_split1 <- D[1:(floor(n/2)), ,]
# D_split2 <- D[(1+floor(n/2)):n, , ]
# D_split1 <- D; D_split2 <- D
# Fit Q_function
# cat('Start Bellman iterations\n')
# transitions_dataset <- make_transitions_dataset(D_split1)
# Q_learning_results <- bellman_iterations(transitions_dataset, evaluation_action_matrix, 
#                                          gamma=gamma, n_iterations=30, V0=V0[1, ])
# cat('Done with Bellman iterations\n')
# plot(Q_learning_results$l2_errors)
# # Replicate accross time points the Bellman iterations based Q-function
# # That would be correct under infinite horizon.
# # Here's it's only approximately correct. Better for early time points than late time points.
# Q_hat <- array(0, dim=c(horizon, n_states, n_actions)); V_hat <- array(0, dim=c(horizon, n_states))
# for(t in 1:horizon){
#   Q_hat[t, ,] <-  Q_learning_results$Q_hat
#   V_hat[t, ] <-  Q_learning_results$V_hat
# }
b <- 5e-1*rnorm(1)
# b <- 0
Q_hat <- Q0 +  b
V_hat <-  V0 +  b

MAGIC_estimates <- c();
WDR_estimates <- c(); LTMLE_estimates <- c(); MAGIC_LTMLE_estimates <- c(); MAGIC_LTMLE_hacky_estimates <- c()
C_TMLE_estimates <- c()
for(i in 1:50){
  D <- generate_discrete_MDP_dataset(n, 1, state_transition_matrix,
                                     behavior_action_matrix,
                                     transition_based_rewards,
                                     horizon)
  MAGIC_estimates <- c(MAGIC_estimates, 
                       MAGIC_estimator(D, Q_hat, V_hat, gamma = gamma, horizon = horizon, n_bootstrap = 100)$estimate)
  WDR_estimates <- c(WDR_estimates, 
                     WDR_estimator_TB_old(D, Q_hat, V_hat, gamma = gamma, j = horizon))
  LTMLE_estimates <- c(LTMLE_estimates,
                       LTMLE_estimator(D, Q_hat, V_hat, evaluation_action_matrix, gamma, alpha=0.1)$estimate)
  C_TMLE_estimates <- c(C_TMLE_estimates,
                        C_LTMLE_softening(D, Q_hat, V_hat, evaluation_action_matrix, gamma, plot_risks=F)$estimate)
  # MAGIC_LTMLE_hacky_estimates <- c(MAGIC_LTMLE_hacky_estimates, 
                                   # MAGIC_LTMLE_estimator_hacky(D, Q_hat, V_hat, evaluation_action_matrix, gamma = gamma, n_bootstrap = 100))
  # MAGIC_LTMLE_estimates <- c(MAGIC_LTMLE_estimates, 
                             # MAGIC_LTMLE_estimator(D, Q_hat, V_hat, evaluation_action_matrix, gamma = gamma, n_bootstrap = 20)$estimate)
}

cat('WDR MSE:', mean( (WDR_estimates - V0[1,1])^2), '\n')
cat('MAGIC MSE:', mean( (MAGIC_estimates - V0[1,1])^2), '\n')
cat('LTMLE MSE:', mean( (LTMLE_estimates - V0[1,1])^2), '\n')
cat('C-LTMLE MSE:', mean( (C_TMLE_estimates - V0[1,1])^2), '\n')

# cat('MAGIC LTMLE hacky MSE:', mean( (MAGIC_LTMLE_hacky_estimates - V0[1,1])^2), '\n')
# cat('MAGIC LTMLE MSE:', mean( (MAGIC_LTMLE_estimates - V0[1,1])^2), '\n')
cat('\n')
cat('WDR bias', mean(WDR_estimates - V0[1,1]), '\n')
cat('MAGIC bias',  mean(MAGIC_estimates - V0[1,1]), '\n')
cat('LTMLE bias:', mean(LTMLE_estimates - V0[1,1]), '\n')
cat('C-LTMLE bias:', mean(C_TMLE_estimates - V0[1,1]), '\n')
# cat('MAGIC LTMLE bias:', mean(MAGIC_LTMLE_estimates - V0[1,1]), '\n')
cat('\n')
cat('WDR variance', var(WDR_estimates), '\n')
cat('MAGIC variance',  var(MAGIC_estimates), '\n')
cat('LTMLE variance:', var(LTMLE_estimates), '\n')
cat('C-LTMLE variance:', var(C_TMLE_estimates), '\n')
# cat('MAGIC LTMLE variance:', var(MAGIC_LTMLE_estimates), '\n')
library(MASS)
source('MDP_modelWin.R')
source('Estimators.R')

# Debugging experiments ---------------------------------------------------
# set.seed(1)
horizon <- 100; n <- 100; gamma <- 0.8
V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
                                  transition_based_rewards,
                                  evaluation_action_matrix, horizon, gamma = gamma)
V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0
Q_hat <- Q0; V_hat <- V0
D <- generate_discrete_MDP_dataset(n, 1, state_transition_matrix,
                                   behavior_action_matrix,
                                   transition_based_rewards,
                                   horizon)

# Save D, V0 and Q0
D_flattened <- array(D, dim=c(dim(D)[1] * dim(D)[2], dim(D)[3]) )
Q_hat_flattened <- array(Q_hat, dim=c(dim(Q_hat)[1] * dim(Q_hat)[2], dim(Q_hat)[3]) )
write.matrix(D_flattened, 'data.csv', sep=',')
write.matrix(Q_hat_flattened, 'Q_hat.csv', sep=',')
write.matrix(V_hat, 'V_hat.csv', sep=',')

cat('True V0: ', V0[1, 1], '\n')
# cat('IS: ', IS_estimator(D), '\n')
# cat('WIS: ', WIS_estimator(D), '\n')
# cat('stepIS: ', stepIS_estimator(D), '\n')
# cat('stepWIS: ', stepWIS_estimator(D, horizon), '\n')
# cat('DR: ', DR_estimator_TB(D, Q_hat=Q0, V_hat=V0), '\n')
# cat('WDR old: ', WDR_estimator_TB_old(D, Q_hat=Q0, V_hat=V0, gamma=gamma, j=horizon), '\n')
# cat('WDR: ', WDR_estimator_TB(D, Q_hat=Q0, V_hat=V0, gamma=gamma, j=horizon)$g_js[horizon+1], '\n')
# cat('LTMLE: ', LTMLE_estimator(D, Q_hat=Q0, V_hat=V0, evaluation_action_matrix=evaluation_action_matrix, gamma=gamma), '\n')

# cat('g^(5):', WDR_estimator_TB(D, Q_hat=Q0, V_hat=V0, gamma=gamma, j=5)$g_js[6], '\n')
# cat('WDR(j=5)', WDR_estimator_TB_old(D, Q_hat=Q0, V_hat=V0, gamma=gamma, j=5), '\n')
# 
# cat('g^(4):', WDR_estimator_TB(D, Q_hat=Q0, V_hat=V0, gamma=gamma, j=4)$g_js[5], '\n')
# cat('WDR(j=4)', WDR_estimator_TB_old(D, Q_hat=Q0, V_hat=V0, gamma=gamma, j=4), '\n')
# 
# cat('g^(3):', WDR_estimator_TB(D, Q_hat=Q0, V_hat=V0, gamma=gamma, j=3)$g_js[4], '\n')
# cat('WDR(j=3)', WDR_estimator_TB_old(D, Q_hat=Q0, V_hat=V0, gamma=gamma, j=3), '\n')
# 
# cat('g^(2):', WDR_estimator_TB(D, Q_hat=Q0, V_hat=V0, gamma=gamma, j=2)$g_js[3], '\n')
# cat('WDR(j=2)', WDR_estimator_TB_old(D, Q_hat=Q0, V_hat=V0, gamma=gamma, j=2), '\n')
# 
# cat('g^(1):', WDR_estimator_TB(D, Q_hat=Q0, V_hat=V0, gamma=gamma, j=1)$g_js[2], '\n')
# cat('WDR(j=1)', WDR_estimator_TB_old(D, Q_hat=Q0, V_hat=V0, gamma=gamma, j=1), '\n')
# 
# cat('g^(0):', WDR_estimator_TB(D, Q_hat=Q0, V_hat=V0, gamma=gamm, j=0)$g_js[1], '\n')
# cat('WDR(j=0)', WDR_estimator_TB_old(D, Q_hat=Q0, V_hat=V0, gamma=gamma, j=0), '\n')

# debug(WDR_estimator_TB)


# Q learning of an initial estimator --------------------------------------
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

# Run Monte Carlo simulation ----------------------------------------------
# Some small MSE experiment
# DR_estimates <- c(); 
WDR_estimates <- c(); LTMLE_estimates <- c()
n <- 100; gamma <- 1
for(i in 1:500){
  D <- generate_discrete_MDP_dataset(n, 1, state_transition_matrix,
                                     behavior_action_matrix,
                                     transition_based_rewards,
                                     horizon)
  # DR_estimates <- c(DR_estimates,
                    # DR_estimator_TB(D, Q_hat=Q_hat, V_hat=Q_hat))
  WDR_estimates <- c(WDR_estimates,
                     WDR_estimator_TB(D, Q_hat=Q_hat, V_hat=V_hat, gamma=gamma, j=horizon)$g_js[horizon+1])
  LTMLE_estimates <- c(LTMLE_estimates,
                       LTMLE_estimator(D, Q_hat=Q_hat, V_hat=Q_hat, evaluation_action_matrix=evaluation_action_matrix, gamma=gamma))
}

# cat('DR MSE:', mean( (DR_estimates - V0[1,1])^2 ), '\n')
cat('WDR MSE:', mean((WDR_estimates - V0[1,1])^2), '\n')
cat('LTMLE MSE:', mean((LTMLE_estimates - V0[1,1])^2), '\n')

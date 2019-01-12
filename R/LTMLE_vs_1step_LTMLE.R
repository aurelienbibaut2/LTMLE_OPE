source('Estimators.R')
source('single_epsilon_LTMLE.R')

n_states <- 3; n_actions <- 2
# Set DGP parameters
horizon <- 10; n <- 1e3; gamma <- 1
V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
                                  transition_based_rewards,
                                  evaluation_action_matrix, horizon, gamma)
# Compute truth
V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0
b <- 1e-1*rnorm(1)
# b <- 0
Q_hat <- Q0 +  b
V_hat <-  V0 +  b

D <- generate_discrete_MDP_dataset(n, 1, state_transition_matrix,
                                   behavior_action_matrix,
                                   transition_based_rewards,
                                   horizon)


one_step_result <- one_step_LTMLE(D, Q_hat, V_hat, evaluation_action_matrix, gamma=1, verbose=T, max_it=20, step_size=5e-2)
LTMLE_result <- LTMLE_estimator(D, Q_hat, V_hat, evaluation_action_matrix, gamma=1, alpha=1)

cat('1step EIC eq:', mean(evaluate_EIC(D, rep(one_step_result$epsilon, horizon), Q_hat, V_hat, evaluation_action_matrix, gamma)), '\n')
cat('LTMLE EIC eq:', mean(evaluate_EIC(D, LTMLE_result$epsilons, Q_hat, V_hat, evaluation_action_matrix, gamma)), '\n')

D_large <- generate_discrete_MDP_dataset(1e4, 1, state_transition_matrix,
                                         behavior_action_matrix,
                                         transition_based_rewards,
                                         horizon)
cat('1step var(EIC):', var(evaluate_EIC(D_large, rep(one_step_result$epsilon, horizon), Q_hat, V_hat, evaluation_action_matrix, gamma)), '\n')
cat('LTMLE var(EIC):', var(evaluate_EIC(D_large, LTMLE_result$epsilons, Q_hat, V_hat, evaluation_action_matrix, gamma)), '\n')

cat('LTMLE epsilon[1]:', LTMLE_result$epsilons[1], '\n')
cat('1 step epsilon:', one_step_result$epsilon, '\n')

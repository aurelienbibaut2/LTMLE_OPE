source('MDP_modelWin.R')
source('Estimators.R')
source('Q_learning_discrete_state_space.R')

set.seed(1)
horizon <- 50; n <- 1e3; gamma = 0.9
V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
                                  transition_based_rewards,
                                  evaluation_action_matrix, horizon, gamma=gamma)
V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0

# Perform Bellman iterations
ns <- c(100)
n_its <- c(); final_l2_errors <- c()

for(n in ns){
  D <- generate_discrete_MDP_dataset(n, 1, state_transition_matrix,
                                     behavior_action_matrix,
                                     transition_based_rewards,
                                     horizon)
  transitions_dataset <- make_transitions_dataset(D)
  t0 <- Sys.time()
  bellman_results <- bellman_iterations(transitions_dataset, evaluation_action_matrix, gamma=gamma, max_it=200, V0=V0[1, ],
                                        relative_tol = 0.0001, verbose=T, start_new_plot=(n == ns[1]))
  t1 <- Sys.time()
  print(t1 - t0)
  n_its <- c(n_its, bellman_results$n_it)
  final_l2_errors <- c(final_l2_errors, bellman_results$l2_errors[bellman_results$n_it])
}

# plot(ns, n_its)
# plot(final_l2_errors)
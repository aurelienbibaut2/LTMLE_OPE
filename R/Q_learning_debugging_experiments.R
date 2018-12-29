source('MDP_modelWin.R')
source('Estimators.R')
source('Q_learning_discrete_state_space.R')

set.seed(1)
horizon <- 100; n <- 10; gamma = 0.9
V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
                                  transition_based_rewards,
                                  evaluation_action_matrix, horizon, gamma=gamma)
V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0

Q_hat <- Q0; V_hat <- V0
D <- generate_discrete_MDP_dataset(n, 1, state_transition_matrix,
                                   behavior_action_matrix,
                                   transition_based_rewards,
                                   horizon)

# Get transitions dataset
transitions_dataset <- make_transitions_dataset(D)

# Perform Bellman iterations
bellman_results <- bellman_iterations(transitions_dataset, evaluation_action_matrix, gamma=gamma, n_iterations=40, V0=V0[1, ])
plot(1:length(bellman_results$l2_errors), bellman_results$l2_errors)


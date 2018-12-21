source('MDP_modelWin.R')
source('Estimators.R')

# Debugging experiments ---------------------------------------------------
horizon <- 15; n <- 1e3
V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
                                  transition_based_rewards,
                                  evaluation_action_matrix, horizon)
V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0
Q_hat <- Q0; V_hat <- V0
D <- generate_discrete_MDP_dataset(n, 1, state_transition_matrix,
                                   behavior_action_matrix,
                                   transition_based_rewards,
                                   horizon)
cat('True V0: ', V0[1, 1], '\n')
cat('IS: ', IS_estimator(D), '\n')
cat('WIS: ', WIS_estimator(D), '\n')
cat('stepIS: ', stepIS_estimator(D), '\n')
cat('stepWIS: ', stepWIS_estimator(D), '\n')
cat('DR: ', DR_estimator(D, Q_hat=Q0, V_hat=V0), '\n')
cat('LTMLE: ', LTMLE_estimator(D, Q_hat=Q0, V_hat=V0), '\n')
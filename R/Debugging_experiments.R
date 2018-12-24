source('MDP_modelWin.R')
source('Estimators.R')

# Debugging experiments ---------------------------------------------------
set.seed(1)
horizon <- 5; n <- 100
V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
                                  transition_based_rewards,
                                  evaluation_action_matrix, horizon)
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
cat('IS: ', IS_estimator(D), '\n')
cat('WIS: ', WIS_estimator(D), '\n')
cat('stepIS: ', stepIS_estimator(D), '\n')
cat('stepWIS: ', stepWIS_estimator(D), '\n')
cat('DR: ', DR_estimator_TB(D, Q_hat=Q0, V_hat=V0), '\n')
cat('LTMLE: ', LTMLE_estimator(D, Q_hat=Q0, V_hat=V0), '\n')


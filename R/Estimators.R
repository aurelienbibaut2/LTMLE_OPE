D <- generate_discrete_MDP_dataset(7, 1, state_transition_matrix,
                                          behavior_action_matrix,
                                          transition_based_rewards,
                                          horizon)

IS_estimator <- function(D){
  horizon <- dim(D)[2]
  mean(D[, horizon, 'rho_t'] * apply(D[, , 'r'], 1, sum))
}

WIS_estimator <- function(D){
  horizon <- dim(D)[2]
  w_H <- mean(D[, horizon, 'rho_t'])
  mean(D[, horizon, 'rho_t'] / w_H * apply(D[, , 'r'], 1, sum))
}

stepIS_estimator <- function(D){
  horizon <- dim(D)[2]
  mean(apply(D[, , 'r'] * D[, , 'rho_t'], 1, sum))
}

WIS_estimator <- function(D){
  n <- dim(D)[1]
  horizon <- dim(D)[2]
  w_t <- apply(D[, , 'rho_t'], 2, mean)
  mean(apply(D[, , 'r'] * D[, , 'rho_t'] / (rep(1, n) %*% t(w_t)), 1, sum))
}


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

stepWIS_estimator <- function(D){
  n <- dim(D)[1]
  horizon <- dim(D)[2]
  w_t <- apply(D[, , 'rho_t'], 2, mean)
  mean(apply(D[, , 'r'] * D[, , 'rho_t'] / (rep(1, n) %*% t(w_t)), 1, sum))
}

DR_estimator <- function(D, Q_hat, V_hat){
  n <- dim(D)[1]
  horizon <- dim(D)[2]
  V_DR <- matrix(NA, nrow=n, ncol=horizon)
  V_DR[, t] <- (V_hat[t, D[, t, 's']]
                +  D[, t, 'pi_a'] / D[, t, 'pi_b'] * (D[, t, 'r'] - apply(D[, t, ], 1, function(x) Q_hat[t, x['s'], x['a']]))
                )
  for(t in (horizon-1):1){
    V_DR[, t] <- (V_hat[t, D[, t, 's']] 
                  +  D[, t, 'pi_a'] / D[, t, 'pi_b'] * (D[, t, 'r'] + V_DR[, t+1] - apply(D[, t, ], 1, function(x) Q_hat[t, x['s'], x['a']]))
                  )
  }
  mean(V_DR[, 1])
}

WDR_estimator <- function(D, Q_hat, V_hat){
  n <- dim(D)[1]
  horizon <- dim(D)[2]
  V_DR <- matrix(NA, nrow=n, ncol=horizon)
  V_DR[, t] <- (V_hat[t, D[, t, 's']]
                +  D[, t, 'pi_a'] / D[, t, 'pi_b'] * (D[, t, 'r'] - apply(D[, t, ], 1, function(x) Q_hat[t, x['s'], x['a']]))
  )
  for(t in (horizon-1):1){
    V_DR[, t] <- (V_hat[t, D[, t, 's']] 
                  +  D[, t, 'pi_a'] / D[, t, 'pi_b'] * (D[, t, 'r'] + V_DR[, t+1] - apply(D[, t, ], 1, function(x) Q_hat[t, x['s'], x['a']]))
    )
  }
  mean(V_DR[, 1])
}

# Debugging experiments ---------------------------------------------------
source('MDP_modelWin.R')
horizon <- 5; n <- 1e4
V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix, 
                                  transition_based_rewards, 
                                  evaluation_action_matrix, horizon)
V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0

D <- generate_discrete_MDP_dataset(n, 1, state_transition_matrix,
                                   behavior_action_matrix,
                                   transition_based_rewards,
                                   horizon)

cat('True V0: ', V0[1, 1], '\n')
cat('IS: ', IS_estimator(D), '\n')
cat('WIS: ', WIS_estimator(D), '\n')
cat('stepIS: ', stepIS_estimator(D), '\n')
cat('stepWIS: ', stepWIS_estimator(D), '\n')
cat('DR: ', DR_estimator(D, Q_hat=Q0, V_hat=V0))
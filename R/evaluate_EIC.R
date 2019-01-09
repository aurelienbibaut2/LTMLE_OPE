# Evaluate EIC, potentially on a split that is not the one one which epsilons have been fitted
evaluate_EIC <- function(D, epsilons, Q_hat, V_hat, evaluation_action_matrix, gamma){
  n <- dim(D)[1]; horizon <- dim(D)[2]
  D_star <- matrix(0, nrow=n, ncol=horizon)
  V_t_evaluated <- rep(0, n)
  Delta_t <- 0
  for(t in horizon:1){
    Delta_t <- 1 + gamma * Delta_t
    # Compute the fluctuated estimator \hat{Q}^*_t from \hat{Q}_t and epsilon_t
    Q_tilde_t_star <- expit( logit( (Q_hat[t, ,] + Delta_t) / (2 * Delta_t) )
                             + epsilons[t])
    V_tilde_t_star <- as.vector(apply(Q_tilde_t_star * evaluation_action_matrix, 1, sum))
    Q_t_star <- 2 * Delta_t * (Q_tilde_t_star - 1/2)
    V_t_star <- 2 * Delta_t * (V_tilde_t_star - 1/2)
    # Compute component t of the EIC
    # Note that at this point we haven't updated V_t_evaluated yet, so it is still V_{t+1}(S_{t+1})
    D_star[, t] <- gamma^t * D[, t, 'rho_t'] * (D[, t, 'r'] + gamma * V_t_evaluated
                                                - apply(D[, t, ], 1, function(x) Q_t_star[x['s'], x['a']]))
    # Evaluate V
    V_t_evaluated <-  V_t_star[D[, t, 's']]
  }
  apply(D_star, 1, sum)
}

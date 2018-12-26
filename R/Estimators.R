AM_estimator <- function(D,V_hat,d_hat=NA,t){ # Thomas and Brunskill's Formula for Approximate Model pg.4
  #Use data to get the estimate of d and V_hat
  if(is.na(d_hat)){
    d_hat<-data.frame(s=unique(as.vector(D[,,'s'])))  
    inter<-table(D[, 1, 's'])/length(D[, 1, 's'])
    d_hat[(d_hat$s %in% attr(inter, which = "dimnames")[[1]]),2]<-table(D[, 1, 's'])/length(D[, 1, 's'])
    d_hat[is.na(d_hat)]<-0
  }
  sum(V_hat[t,]*d_hat$V2)
}

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

stepWIS_estimator <- function(D,t){
  n <- dim(D)[1]
  w_t <- apply(D[, , 'rho_t'], 2, mean)
  mean(apply(D[, , 'r'][,1:t] * D[, , 'rho_t'][,1:t] / (rep(1, n) %*% t(w_t[1:t])), 1, sum))
}

DR_estimator_JL <- function(D, Q_hat, V_hat){ # Jiang and Li's DR estimator, based on a recursive formula
  n <- dim(D)[1]
  horizon <- dim(D)[2]
  V_DR <- matrix(NA, nrow=n, ncol=horizon)
  t <- horizon
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

DR_estimator_TB <- function(D, Q_hat, V_hat){ # Thomas and Brunskill's DR estimator. Turns out that it is exactly the same as Jiang and Li's
  n <- dim(D)[1]
  horizon <- dim(D)[2]
  D_star <- matrix(NA, nrow=n, ncol=horizon)
  t <- horizon
  D_star[, t] <- D[, t, 'rho_t'] * (D[, t, 'r'] + 
                                    - apply(D[, t, ], 1, function(x) Q_hat[t, x['s'], x['a']]))
  for(t in (horizon-1):1){
    D_star[, t] <- D[, t, 'rho_t'] * (D[, t, 'r'] + V_hat[t, D[, t+1, 's']]
                                                      - apply(D[, t, ], 1, function(x) Q_hat[t, x['s'], x['a']]))
  }
  mean(V_hat[t, D[, 1, 's']] + apply(D_star, 1, sum))
}

WDR_estimator_TB <- function(D, Q_hat, V_hat){ # Thomas and Brunskill's Weighted DR estimator
  n <- dim(D)[1]
  horizon <- dim(D)[2]
  D_star <- matrix(NA, nrow=n, ncol=horizon)
  t <- horizon
  
  w_t <- apply(D[, , 'rho_t'], 2, mean)

  D_star[, t] <- D[, t, 'rho_t'] / w_t[t] * (D[, t, 'r'] + 
                                      - apply(D[, t, ], 1, function(x) Q_hat[t, x['s'], x['a']]))
  for(t in (horizon-1):1){
    D_star[, t] <- D[, t, 'rho_t'] / w_t[t] * (D[, t, 'r'] + V_hat[t, D[, t+1, 's']]
                                      - apply(D[, t, ], 1, function(x) Q_hat[t, x['s'], x['a']]))
  }
  mean(V_hat[t, D[, 1, 's']] + apply(D_star, 1, sum))
}

# Two helper functions
# Expit and logit
expit <- function(x) { 1 / (1 + exp(-x)) } 
logit <- function(x) { log(x / (1 - x)) }

# LTMLE
LTMLE_estimator <-  function(D, Q_hat, V_hat){
  # Get dataset dimensions
  n <- dim(D)[1]
  horizon <- dim(D)[2]
  
  epsilons <- c()
  V_evaluated <- rep(0, n) #V_{t+1}(S_{t+1})
  for(t in horizon:1){
    Delta_t <- horizon + 1 - t
    R <- D[, t, 'r'] # R_t
    U_tilde <- (R + V_evaluated + Delta_t) / (2 * Delta_t) # U_tilde = R_tilde_t + V_tilde_{t+1}(S_t) in the notations of the write-up
    Q_t_evaluated <- apply(D[ , t, ], 1, function(x) Q_hat[t, x['s'], x['a']]) # Q_t(A_t, S_t)
    Q_tilde_t_evaluated <- (Q_t_evaluated + Delta_t) / (2 * Delta_t)
    epsilon <- glm(U_tilde ~ offset(logit(Q_tilde_t_evaluated)) + 1, 
                   family=quasibinomial, weights = D[,t, 'rho_t'])$coefficients[1]
    epsilons <- c(epsilons, epsilon)
    # Evaluate Q_tilde(s_t, a) for a_t = 1, a_t = 2
    Q_tilde_t_star <- expit(logit((Q_hat[t, ,] + Delta_t) / (2 * Delta_t)) + epsilon) # Q_tilde_t^*
    # Then set V_tilde(s_t) = \sum_{a} Q_tilde(s_t, a) pi_a(a|s_t)
    V_tilde <- apply(Q_tilde_t_star * evaluation_action_matrix, 1, sum)
    # Compute V = 2 * Delta_t * (V_tilde - 1)
    V <- 2 * Delta_t * (V_tilde - 1/2)
    # Evaluate V
    V_evaluated <-  V[D[, t, 's']]
    # V_tilde is gonna be V_{t+1} in the next iteration
  }
  # The average of the last V is the LTML estimator of the value
  V_hat_LTMLE <- mean(V_evaluated)
}



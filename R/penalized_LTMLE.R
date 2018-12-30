# Two helper functions
# Expit and logit
expit <- function(x) { 1 / (1 + exp(-x)) } 
logit_scalar <- Vectorize(function(x) { 
  if(x > 0 & x < 1){
    return(log(x / (1 - x)))
  }else if(x <= 0){
    return(-Inf)
  }else{
    return(Inf)
  }
})
logit <- function(x){
  if(is.null(dim(x))){
    return(logit_scalar(x))
  }else{
    return(array(logit_scalar(x), dim=dim(x)))
  }
}

# Soften weights
soften <- function(a, alpha){
  a^alpha / sum(a^alpha) 
}


targeting_loss <- function(epsilon, Delta_t, t, gamma,  D_t, U_tilde, Q_tilde_t_evaluated, V_tplus1_evaluated, penalty){
  Q_tilde_t_epsilon <- expit(logit(Q_tilde_t_evaluated) + epsilon)
  Q_t_epsilon <- 2 * Delta_t * (Q_tilde_t_epsilon  - 1/2)
  D_star_t <- gamma^t * D_t[, 'rho_t'] * (D_t[, 'r'] + gamma * V_tplus1_evaluated
                                          - Q_t_epsilon)
  ll <- mean(D_t[, 'rho_t'] * (U_tilde * log(Q_tilde_t_epsilon) + (1-U_tilde) * log(1 - Q_tilde_t_epsilon) ) )
  
  #-ll + penalty * var(D_star_t)
  -ll + penalty * abs(epsilon)
}



# Evaluate EIC, potentially on a split that is not the one one which epsilons have been fitted
evaluate_EIC <- function(D, epsilons, Q_hat, V_hat, evaluation_action_matrix, gamma, loss_softening_coeff=1){
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
    D_star[, t] <- gamma^t * soften(D[, t, 'rho_t'], loss_softening_coeff) * (D[, t, 'r'] + gamma * V_t_evaluated
                                                                              - apply(D[, t, ], 1, function(x) Q_t_star[x['s'], x['a']]))
    # Evaluate V
    V_t_evaluated <-  V_t_star[D[, t, 's']]
  }
  apply(D_star, 1, sum)
}

# Targeting loss bis
targeting_loss_bis <- function(epsilon, Delta_t, t, gamma, D, 
                               U_tilde, Q_tilde_t_evaluated, V_tplus1_evaluated, 
                               Q_hat, V_hat,
                               evaluation_action_matrix, penalty){
  D_t <- D[, t, ]
  Q_tilde_t_epsilon <- expit(logit(Q_tilde_t_evaluated) + epsilon)
  Q_t_epsilon <- 2 * Delta_t * (Q_tilde_t_epsilon  - 1/2)
  D_star_t <- gamma^t * D_t[, 'rho_t'] * (D_t[, 'r'] + gamma * V_tplus1_evaluated
                                          - Q_t_epsilon)
  ll <- mean(D_t[, 'rho_t'] * (U_tilde * log(Q_tilde_t_epsilon) + (1-U_tilde) * log(1 - Q_tilde_t_epsilon) ) )
  
  epsilons <- c(rep(epsilon, t), rep(0, horizon-t))
  -ll + penalty * var(evaluate_EIC(D, epsilons, Q_hat, V_hat, evaluation_action_matrix, gamma, loss_softening_coeff=1))
  #-ll + penalty * abs(epsilon)
}


# LTMLE
penalized_LTMLE_estimator <-  function(D, Q_hat, V_hat, evaluation_action_matrix, gamma=1, alpha=1, penalty=1){
  # Get dataset dimensions
  n <- dim(D)[1]
  horizon <- dim(D)[2]
  
  epsilons <- c()
  Delta_t <- 0
  V_evaluated <- rep(0, n) #V_{t+1}(S_{t+1})
  for(t in horizon:1){
    Delta_t <- 1 + gamma * Delta_t
    R <- D[, t, 'r'] # R_t
    U_tilde <- (R + gamma * V_evaluated + Delta_t) / (2 * Delta_t) # U_tilde = R_tilde_t + gamma*V_tilde_{t+1}(S_t) in the notations of the write-up
    Q_t_evaluated <- apply(D[ , t, ], 1, function(x) Q_hat[t, x['s'], x['a']]) # Q_t(A_t, S_t)
    Q_tilde_t_evaluated <- (Q_t_evaluated + Delta_t) / (2 * Delta_t)
    # Targeting step
    #epsilon <- glm(U_tilde ~ offset(logit(Q_tilde_t_evaluated)) + 1, 
    #               family=quasibinomial, weights = soften(D[,t, 'rho_t'], alpha) )$coefficients[1]
    loss <- function(epsilon) targeting_loss(epsilon, Delta_t, t, gamma,  
                           D[, t, ], U_tilde, Q_tilde_t_evaluated, V_evaluated, 
                           penalty=penalty)
    loss_bis <- function(epsilon) targeting_loss_bis(epsilon, Delta_t, t, gamma,  
                                                     D, U_tilde, Q_tilde_t_evaluated, V_evaluated, 
                                                     Q_hat, V_hat,
                                                     evaluation_action_matrix, penalty)
    epsilon <- optimize(f = Vectorize(loss_bis), lower=-1, upper=1)$minimum
    epsilons <- c(epsilons, epsilon)
    # Evaluate Q_tilde(s_t, a) for a_t = 1, a_t = 2
    Q_tilde_t_star <- expit( logit( (Q_hat[t, , ] + Delta_t) / (2 * Delta_t)  )
                             + epsilon) # Q_tilde_t^*
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
  list(estimate=V_hat_LTMLE, epsilons=epsilons)
}
# debug(penalized_LTMLE_estimator)
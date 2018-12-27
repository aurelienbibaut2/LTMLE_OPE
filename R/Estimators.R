AM_estimator <- function(D, V_hat, d_hat=NA, t){ # Thomas and Brunskill's Formula for Approximate Model pg.4
  #Use data to get the estimate of d and V_hat
  if(is.na(d_hat)){
    d_hat <- data.frame(s=unique(as.vector(D[,,'s'])))  
    inter <- table(D[, 1, 's']) / length(D[, 1, 's'])
    d_hat[(d_hat$s %in% attr(inter, which = "dimnames")[[1]]),2] <- table(D[, 1, 's'])/length(D[, 1, 's'])
    d_hat[is.na(d_hat)] <- 0
  }
  sum(V_hat[1,]*d_hat$V2)
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

stepWIS_estimator <- function(D, t=NULL){
  if(is.null(t)) t <- dim(D)[2] # By default, compute the WIS estimator for the value of the entire trajectory
  n <- dim(D)[1]
  w_t <- apply(D[, , 'rho_t'], 2, mean)
  mean(apply(D[, 1:t, 'r'] * D[, 1:t, 'rho_t'] / (rep(1, n) %*% t(w_t[1:t])), 1, sum))
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
    D_star[, t] <- D[, t, 'rho_t'] * (D[, t, 'r'] + V_hat[t+1, D[, t+1, 's']]
                                                      - apply(D[, t, ], 1, function(x) Q_hat[t, x['s'], x['a']]))
  }
  mean(V_hat[t, D[, 1, 's']] + apply(D_star, 1, sum))
}

# Thomas and Brunskill's Weighted DR estimator
# Needs to be provided with a Q_hat and a V_hat that have dimension
# horizon x n_states x n_actions and horizon x n_states, respectively
# With the function implemented below, g^(j)(D) is just WDR_estimator(D, Q_hat, V_hat, gamma, j)
WDR_estimator_TB <- function(D, Q_hat, V_hat, gamma=1, j=NULL, gjD=FALSE){ 
  n <- dim(D)[1]
  horizon <- dim(D)[2]
  if(is.null(j)) j <- horizon
  D_star <- matrix(NA, nrow=n, ncol=j)
  # Compute mean rho_t to be used as denominator in the stabilized weights
  w_t <- apply(D[, , 'rho_t'], 2, mean)
  
  # D_star below is computed slightly differently whether t=horizon or t <= horizon-1
  # The reason is that V_hat[horizon+1, s] is zero in reality for every s, but the
  # row V_hat[horizon+1, ] does not exist in the V_hat passed as argument 
  # (there is no need to have a row for V_hat[horizon+1, ] as we know it's zero.
  if(j == horizon){
    t <- j
    D_star[, t] <- gamma^t * D[, t, 'rho_t'] / w_t[t] * (D[, t, 'r'] + 
                                                           - apply(D[, t, ], 1, function(x) Q_hat[t, x['s'], x['a']]))
  }
  
  if(j > 0){
    for(t in (min(horizon-1, j)):1){
      D_star[, t] <- gamma^t * D[, t, 'rho_t'] / w_t[t] * (D[, t, 'r'] + gamma * V_hat[t+1, D[, t+1, 's']]
                                                           - apply(D[, t, ], 1, function(x) Q_hat[t, x['s'], x['a']]))
    }
  }
  if(gjD){
    return(V_hat[1, D[, 1, 's']] + apply(D_star, 1, sum))
  }else{
    mean(V_hat[1, D[, 1, 's']] + apply(D_star, 1, sum))
  }
}

# Two helper functions
# Expit and logit
expit <- function(x) { 1 / (1 + exp(-x)) } 
logit <- function(x) { log(x / (1 - x)) }

# LTMLE
LTMLE_estimator <-  function(D, Q_hat, V_hat, evaluation_action_matrix){
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

##MAGIC

library(quadprog)
library(Matrix)

bootstrap_bias <- function(D, Q0, V0, number_bootstrap = 1e3, alpha = 0.1) {
  n <- dim(D)[1]
  wdr_bootstrap <- c()
  for (i in 1:number_bootstrap) {
    index_bootstrap <- sample(x = 1:n, size = n, replace = TRUE)
    D_bootstrap <- D[index_bootstrap, , ]
    wdr_once <- WDR_estimator_TB(D = D_bootstrap, Q_hat = Q0, V_hat = V0)
    wdr_bootstrap <- c(wdr_bootstrap, wdr_once)
  }
  wdr_quantiles <- quantile(wdr_bootstrap, probs = c(alpha / 2, 1 - alpha / 2))
  return(wdr_quantiles)
}

MAGIC_estimator<-function(D,Q_hat,V_hat,gamma=1,k=1e4,alpha=0.1,N=100){
  
  horizon <- dim(D)[2]
  n <- dim(D)[1]
  
  #g^(j)(D) for until horizon
  #0 is AM-based, horizon IS-based
  gjD<-foreach::foreach(i=0:horizon, .combine = rbind, .inorder=T) %dopar% 
  {WDR_estimator_TB(D=D,Q_hat=Q_hat,V_hat=V_hat,gamma=gamma,j=i,gjD=TRUE)}
  attr(gjD, "dimnames") <- NULL
  gjD_est<-data.frame(apply(gjD,1,mean))
  
  #Calculate sample covariance:
  Sigma<-cov(t(gjD)) #See variance increase as more IS
  
  #Estimate bias:
  CI<-bootstrap_bias(D=D,Q0=Q_hat,V0=V_hat,number_bootstrap=k,alpha=alpha)
  bn<-apply(gjD_est,1,function(x) min(abs(x-CI[1]),abs(x-CI[2])) )
  
  res<-quadprog::solve.QP(Dmat=Matrix::nearPD(Sigma + bn %*% t(bn), eig.tol=1e-10)$mat, dvec = rep(0,nrow(gjD_est)), 
                          Amat = cbind(rep(1,nrow(gjD_est)), diag(nrow(gjD_est))), bvec = c(1,rep(0,nrow(gjD_est))))
  
  as.numeric(t(res$solution) %*% gjD_est[,1])
}

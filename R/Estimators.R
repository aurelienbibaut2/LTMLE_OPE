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

  D_star[, t] <- D[, , 'rho_t'] / (rep(1, n) %*% t(w_t)) * (D[, t, 'r'] + 
                                      - apply(D[, t, ], 1, function(x) Q_hat[t, x['s'], x['a']]))
  for(t in (horizon-1):1){
    D_star[, t] <- D[, , 'rho_t'] / (rep(1, n) %*% t(w_t)) * (D[, t, 'r'] + V_hat[t, D[, t+1, 's']]
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

# Debugging experiments ---------------------------------------------------
# source('MDP_modelWin.R')
# horizon <- 15; n <- 1e3
# V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
#                                   transition_based_rewards,
#                                   evaluation_action_matrix, horizon)
# V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0
# Q_hat <- Q0; V_hat <- V0
# D <- generate_discrete_MDP_dataset(n, 1, state_transition_matrix,
#                                    behavior_action_matrix,
#                                    transition_based_rewards,
#                                    horizon)
# cat('True V0: ', V0[1, 1], '\n')
# cat('IS: ', IS_estimator(D), '\n')
# cat('WIS: ', WIS_estimator(D), '\n')
# cat('stepIS: ', stepIS_estimator(D), '\n')
# cat('stepWIS: ', stepWIS_estimator(D), '\n')
# cat('DR: ', DR_estimator(D, Q_hat=Q0, V_hat=V0), '\n')
# cat('LTMLE: ', LTMLE_estimator(D, Q_hat=Q0, V_hat=V0), '\n')

# Simulations -------------------------------------------------------------
horizon <- 5
V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
                                  transition_based_rewards,
                                  evaluation_action_matrix, horizon)
V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0
# Specify jobs ------------------------------------------------------------
library(foreach); library(doParallel)
nb_repeats <- (detectCores() - 1) * 10
ns <- c(50, 100, 200, 500, 1000
        #, 5000, 10000
        )
jobs <- expand.grid(n = ns, repeat.id = 1:nb_repeats)


cat(detectCores(), 'cores detected\n')
cl <- makeCluster(getOption("cl.cores", detectCores()-1), outfile = '')
registerDoParallel(cl)



results <- foreach(i=1:nrow(jobs), .combine = rbind,
                             #.packages = c('speedglm'),
                             #.export = c('TMLE_EY1_speedglm', 'CV_AIPTW_EY1', 'expit', 'logit', 'g_to_g_delta'),
                             .verbose = T, .inorder = T) %dopar% {
                               D <- generate_discrete_MDP_dataset(jobs[i, ]$n, 1, state_transition_matrix,
                                                                  behavior_action_matrix,
                                                                  transition_based_rewards,
                                                                  horizon)
                               rbind(c(n=jobs[i, ]$n, estimator='IS', estimate=IS_estimator(D)),
                                     c(n=jobs[i, ]$n, estimator='WIS', estimate=WIS_estimator(D)),
                                     c(n=jobs[i, ]$n, estimator='stepIS', estimate=stepIS_estimator(D)),
                                     c(n=jobs[i, ]$n, estimator='stepWIS', estimate=stepWIS_estimator(D)),
                                     c(n=jobs[i, ]$n, estimator='DR',  estimate=DR_estimator_JL(D, Q_hat=Q0, V_hat=V0)),
                                     # c(n=jobs[i, ]$n, estimator='DR_TB',  estimate=DR_estimator_TB(D, Q_hat=Q0, V_hat=V0)),
                                     c(n=jobs[i, ]$n, estimator='WDR',  estimate=WDR_estimator_TB(D, Q_hat=Q0, V_hat=V0)),
                                     c(n=jobs[i, ]$n, estimator='LTMLE', estimate=LTMLE_estimator(D, Q_hat=Q0, V_hat=V0))
                                     )
                             }
results_df <- transform(as.data.frame(results), 
                        n=as.numeric(as.character(n)),
                        estimate=as.numeric(as.character(estimate)))
results_df$squared_error <- (results_df$estimate - V0[1,1])^2

MSE_table <- aggregate(results_df$squared_error, 
                       list(estimator=results_df$estimator, n=results_df$n), 
                       mean)
colnames(MSE_table)[3] <- 'MSE'
library(ggplot2)
MSE_plot <- ggplot(data=MSE_table, aes(x=log10(n), y=log10(MSE), color=estimator)) + 
  geom_line() + geom_point() + ggtitle(paste('ModelWin, horizon=', horizon, ', number of draws per point=', nb_repeats))
print(MSE_plot)
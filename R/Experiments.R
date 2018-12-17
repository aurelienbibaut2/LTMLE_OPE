library(tensorA)
library(here)
source(here("R/MDP_modelWin.R"))

# Experiments -------------------------------------------------------------
horizon <- 2; M <- 1e4
V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix, 
                                  transition_based_rewards, 
                                  evaluation_action_matrix, horizon)
V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0

## Step-wise Importance Sampling. Following paper 2015 by Nan Jiang and Nihong Li
V_IS_summands <- c()
for(m in 1:M){
  H <- generate_discrete_MDP_trajectory(1, state_transition_matrix,
                                        behavior_action_matrix,
                                        transition_based_rewards,
                                        horizon)
  pi_b <- apply(H, 1, function(x) behavior_action_matrix[x[1], x[2]]) # propensity score 
  pi_a <- apply(H, 1, function(x) evaluation_action_matrix[x[1], x[2]]) # probabilities of actions taken in H under the evaluation policy
  
  V_IS <- rep(NA, horizon)
  t <- horizon
  V_IS[t] <- pi_a[t] / pi_b[t] * (H[t, 'r'])
  
  for(t in (horizon-1):1){
    V_IS[t] <- V_IS[t+1] + pi_a[t] / pi_b[t] * (H[t, 'r'])
  }
  V_IS_summands <- c(V_IS_summands, V_IS[1])
}

## Step-wise Weighted Importance Sampling. Following paper 2015 by Nan Jiang and Nihong Li
# Should be made more efficient
V_WIS_summands <- c()
H <- list()
for(m in 1:M){
  H[[m]]<- generate_discrete_MDP_trajectory(1, state_transition_matrix,
                                            behavior_action_matrix,
                                            transition_based_rewards,
                                            horizon)
}

pi_b <- do.call(rbind, lapply(H, function(h) {apply(h, 1, function(x) behavior_action_matrix[x[1], x[2]])}))
pi_a <- do.call(rbind, lapply(H, function(h) {apply(h, 1, function(x) evaluation_action_matrix[x[1], x[2]])}))

for(m in 1:M){
  t <- horizon
  w <- rep(NA, horizon)
  w[t]<-sum((pi_a[1:M,t] / pi_b[1:M,t]) / M)
  
  V_WIS <- rep(NA, horizon)
  
  V_WIS[t] <- (pi_a[m,t] / pi_b[m,t])/w[t] * (H[[m]][t, 'r'])
  
  for(t in (horizon-1):1){
    w[t]<-sum((pi_a[1:M,t] / pi_b[1:M,t]) / M)
    V_WIS[t] <- V_WIS[t+1] + (pi_a[m,t] / pi_b[m,t])/w[t] * (H[[m]][t, 'r'])
  }
  V_WIS_summands <- c(V_WIS_summands, V_WIS[1])
}

# DR estimator. Following 2015 paper by Nan Jiang and Nihong Li
V_DR_summands <- c()
for(m in 1:M){
  H <- generate_discrete_MDP_trajectory(1, state_transition_matrix,
                                        behavior_action_matrix,
                                        transition_based_rewards,
                                        horizon)
  pi_b <- apply(H, 1, function(x) behavior_action_matrix[x[1], x[2]]) # propensity score 
  pi_a <- apply(H, 1, function(x) evaluation_action_matrix[x[1], x[2]]) # probabilities of actions taken in H under the evaluation policy
  
  V_DR <- rep(NA, horizon)
  t <- horizon
  V_DR[t] <- V0[t, H[t, 's']] +  pi_a[t] / pi_b[t] * (H[t, 'r'] - Q0[t, H[t, 's'], H[t, 'a']])
  for(t in (horizon-1):1){
    V_DR[t] <- V0[t, H[t, 's']] +  pi_a[t] / pi_b[t] * (H[t, 'r'] + V_DR[t+1] - Q0[t, H[t, 's'], H[t, 'a']])
  }
  V_DR_summands <- c(V_DR_summands, V_DR[1])
}

# Weighted DR estimator. Following 2016 paper by Philip Thomas and Emma Brunskill
V_WDR_summands <- c()
H <- list()

for(m in 1:M){
  H[[m]]<- generate_discrete_MDP_trajectory(1, state_transition_matrix,
                                            behavior_action_matrix,
                                            transition_based_rewards,
                                            horizon)
}

pi_b <- do.call(rbind, lapply(H, function(h) {apply(h, 1, function(x) behavior_action_matrix[x[1], x[2]])}))
pi_a <- do.call(rbind, lapply(H, function(h) {apply(h, 1, function(x) evaluation_action_matrix[x[1], x[2]])}))

for(m in 1:M){
  
  t <- horizon
  w <- c(1,rep(NA, horizon))
  w[t+1]<-(pi_a[m,t]/pi_b[m,t]) / sum((pi_a[1:M,t] / pi_b[1:M,t])^2)
  
  for(i in (horizon-1):1){
    w[i+1]<-(pi_a[m,i]/pi_b[m,i]) / sum((pi_a[1:M,i] / pi_b[1:M,i])^2)
  }
  
  V_WDR <- rep(NA, horizon)
  V_WDR[t] <- w[t+1]*H[[m]][t, 'r'] - (w[t+1]*Q0[t, H[[m]][t, 's'], H[[m]][t, 'a']] - w[t]* V0[t, H[[m]][t, 's']])

  for(t in (horizon-1):1){
    V_WDR[t] <- V_WDR[t+1] + w[t+1]*H[[m]][t, 'r'] - (w[t+1]*Q0[t, H[[m]][t, 's'], H[[m]][t, 'a']] - w[t]* V0[t, H[[m]][t, 's']])
  }
  V_WDR_summands <- c(V_WDR_summands, V_WDR[1])
 
}

MC_V_IS <- mean(V_IS_summands)
MC_V_WIS <- mean(V_WIS_summands)
MC_V_DR <- mean(V_DR_summands)
MC_V_WDR <- mean(V_WDR_summands)

cat('MC_step_IPS: ', MC_step_IPS(state_transition_matrix,
                                 behavior_action_matrix,
                                 transition_based_rewards, 
                                 evaluation_action_matrix,
                                 horizon, M=M), '\n')
cat('MC_direct_evaluation: ', MC_direct_evaluation(state_transition_matrix,
                                                   transition_based_rewards, 
                                                   evaluation_action_matrix,
                                                   horizon, M=M), '\n')
cat('MC_DR: ', MC_V_DR, '\n')
cat('MC_IS: ', MC_V_IS, '\n')
cat('MC_WIS: ', MC_V_WIS, '\n')
cat('MC_WDR: ', MC_V_WDR, '\n')
cat('V0: ', V0[1,1], '\n')
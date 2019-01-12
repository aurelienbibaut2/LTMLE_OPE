#Simple empirical model-based learning
library(plyr)

est_state_trans<-function(D,s=NA,a=NA,horizon=NA){
  
  if(is.na(s)){
    s<-length(unique(as.vector(D[,,1])))
  }
  if(is.na(a)){
    a<-length(unique(as.vector(D[,,2])))
  }
  if(is.na(horizon)){
    horizon<-dim(D[1,,])[1]
  }

  state_transition_matrix <- array(data=NA, dim=c(s, a, s))
  transition_based_rewards <- array(data=NA, dim=c(s, s))
  exp<-expand.grid(s=1:s,a=1:a,s_new=1:s)
  exp_rew<-expand.grid(s=1:s,s_new=1:s)
  
  prob<-NULL
  for(i in 1:(horizon-1)){
    new<-cbind(D[,i,][,1:3],s_new=D[,i+1,][,1])
    prob<-rbind(prob,new)
  }
  
  reward<-prob[,c(1,4,3)]
  prob<-prob[,-3]
  
  p<-merge(plyr::count(prob, 1:ncol(prob)),plyr::count(prob[,1:2], 1:ncol(prob[,1:2])),by=c("x.s","x.a"))
  names(p)<-c("s","a","s_new","freq.x","freq.y")
  p$prob<-p$freq.x/p$freq.y
  
  #If transition not observed in the data, assign 0
  all<-merge(exp, p, by=c("s","a","s_new"), all=TRUE)
  all[is.na(all$prob),"prob"]<-0

  for(i in 1:s){
    state_transition_matrix[i,,] <- t(matrix(all[all$s==i,"prob"], ncol = a, nrow=s))
  }
  
  #Approximate reward
  r<-merge(plyr::count(reward, 1:ncol(reward)),plyr::count(reward[,1:2], 1:ncol(reward[,1:2])),by=c("x.s","x.s_new"))
  names(r)<-c("s","s_new","r","freq.x","freq.y")
  r$reward<-r$r*r$freq.x/r$freq.y
  
  #If reward not observed in the data, assign 0
  all_r<-merge(exp_rew, r, by=c("s","s_new"), all=TRUE)
  all_r[is.na(all_r$reward),"reward"]<-0
  
  for(i in 1:s){
    transition_based_rewards[i,] <- all_r[all_r$s==i,"reward"]
  }

  return(list(state_transition_matrix=state_transition_matrix, 
         transition_based_rewards=transition_based_rewards))
}

compute_V_and_Q <- function(est_state_transition,
                            est_transition_rewards,
                            evaluation_action_matrix, horizon, gamma=1){
  
  s<-dim(est_transition_rewards)[1]
  new_s<-dim(est_transition_rewards)[2]
  a<-dim(evaluation_action_matrix)[2]

  state_transition_tensor <- as.tensor(est_state_transition, dims=c(s=s, a=a, new_s=s))
  rewards_tensor <- as.tensor(est_transition_rewards, dims=c(s=s, new_s=s))
  pi_tensor <- as.tensor(evaluation_action_matrix, dims=c(s=s, a=a))
  
  # expected reward given s and a
  ER <- reorder(mul.tensor(state_transition_tensor, i='new_s', 
                           rewards_tensor, j='new_s', by=c('s', 'a')), i='s') 
  
  V0 <- matrix(NA, nrow=horizon, ncol=s) 
  Q0 <- array(NA, dim=c(horizon, s, a))
  Q0[horizon, ,] <- ER
  V0[horizon, ] <- mul.tensor(ER, i='a', pi_tensor, j='a', by='s')
  
  for(t in (horizon-1):1){
    V0[t, ] = (mul.tensor(ER, i='a', pi_tensor, j='a', by='s')
               + gamma * mul.tensor(X=pi_tensor, i='a',
                                    Y=mul.tensor(state_transition_tensor, i='new_s', 
                                                 as.tensor(V0[t+1, ], dims=c('new_s'=s)), j='new_s',
                                                 by=c('s', 'a')),
                                    j='a',
                                    by='s'))
    Q0[t, , ] <- (ER + gamma * mul.tensor(state_transition_tensor, i='new_s',
                                          as.tensor(V0[t+1, ], dims=c('new_s'=s)), j='new_s',
                                          by=c('s', 'a')))
  }
  list(V0=V0, Q0=Q0)
}

#library(here)
#setwd(here("R"))

#Test
#
source('MDP_modelWin.R')
n=100; horizon=10
D <- generate_discrete_MDP_dataset(n=n, 1, state_transition_matrix,behavior_action_matrix,
                                   transition_based_rewards,horizon=horizon)
est<-est_state_trans(D)
est_V_and_Q<-compute_V_and_Q(est_state_transition=est$state_transition_matrix,
                est_transition_rewards=est$transition_based_rewards,
                evaluation_action_matrix=evaluation_action_matrix, horizon=horizon, gamma=1)

V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,transition_based_rewards,
                                  evaluation_action_matrix, horizon, gamma=1)

#source('MDP_modelFail.R')
#n=100; horizon=5
#D <- generate_discrete_MDP_dataset(n=n, 1, state_transition_matrix, behavior_action_matrix,
#                                   transition_based_rewards, horizon=horizon)
#est<-est_state_trans(D)

#source('MDP_Gridworld.R')
#n=1000; horizon=100
#D <- generate_discrete_MDP_dataset(n=n, 1, state_transition_matrix, behavior_action_matrix,
#                                   transition_based_rewards, horizon=horizon)
#est<-est_state_trans(D)


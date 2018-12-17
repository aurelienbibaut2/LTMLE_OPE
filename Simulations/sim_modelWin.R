library(tensorA)
library(here)
source(here("R/MDP_modelWin.R"))
source(here("R/utils.R"))

#Set simulaton parameters:
M=100
MC=1000
horizon=20
noiseb=0
noisee=0
noiseQ=0
noiseV=0
res<-list()

#Get the true Q0 and V0:
V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix, 
                                  transition_based_rewards, 
                                  evaluation_action_matrix, horizon)
V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0

for(i in 1:MC){
  
  # Generate M trajectories with specified 
  # state_transition_matrix, behavior_action_matrix and transition_based_rewards
  # with specified horizon 
  H <- gen_traj(M=M, state_transition_matrix,
                behavior_action_matrix,
                transition_based_rewards,
                horizon)
  
  # Run all the estimators (add ltmle)
  res[[i]] <- RL_sim(H=H, state_transition_matrix, behavior_action_matrix, transition_based_rewards,
                horizon, Q0, V0, noiseb, noisee, noiseQ, noiseV)
  
  
  
}




  

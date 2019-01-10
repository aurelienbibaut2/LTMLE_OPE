library(tensorA)
# Generate a discrete MDP trajectory
generate_discrete_MDP_trajectory <- function(s0, state_transition_matrix,
                                             behavior_action_matrix,
                                             transition_based_rewards,
                                             horizon){
  states_trajectory <- c(); actions_trajectory <- c(); rewards_trajectory <- c()
  s <- s0
  for(t in 1:horizon){
    if(s==16){
      #So all trajectories are the same length 
      a <- which(rmultinom(n=1, size=1, prob=behavior_action_matrix[s, ]) == 1)
      new_s <- 16
      r <- 0
      states_trajectory <- c(states_trajectory, 16)
      actions_trajectory <- c(actions_trajectory, a)
      rewards_trajectory <- c(rewards_trajectory, r)
      s <- new_s
    }else{
      a <- which(rmultinom(n=1, size=1, prob=behavior_action_matrix[s, ]) == 1)
      new_s <- which(rmultinom(n=1, size=1, prob=state_transition_matrix[s, a, ]) == 1)
      r <- transition_based_rewards[s, new_s]
      states_trajectory <- c(states_trajectory, s)
      actions_trajectory <- c(actions_trajectory, a)
      rewards_trajectory <- c(rewards_trajectory, r)
      s <- new_s
    }
    
  }
  cbind(s=states_trajectory, a=actions_trajectory, r=rewards_trajectory)
}

# Generate discrete trajectory with probabilities of transitions under both
# the behavior and the exploration policy
generate_discrete_MDP_trajectory_with_pi_a_pi_b <- function(s0, state_transition_matrix,
                                                            behavior_action_matrix,
                                                            transition_based_rewards,
                                                            horizon, gamma=1){
  H <- generate_discrete_MDP_trajectory(s0, state_transition_matrix,
                                        behavior_action_matrix,
                                        transition_based_rewards,
                                        horizon)
  pi_b <- apply(H, 1, function(x) behavior_action_matrix[x[1], x[2]]) # propensity score 
  pi_a <- apply(H, 1, function(x) evaluation_action_matrix[x[1], x[2]]) # probabilities of actions taken in H under the evaluation policy
  rho_t <- cumprod(pi_a/pi_b)
  rtg <- rev(cumsum(rev(gamma^(1:nrow(H)) * H[, 'r']))) / (gamma^(1:nrow(H)))
  gamma_t <- gamma^(1:nrow(H))
  cbind(H, pi_a=pi_a, pi_b=pi_b, rho_t=rho_t, rtg=rtg, gamma_t=gamma_t)
}

# Generate dataset, that is a bundle of trajectories
generate_discrete_MDP_dataset <- function(n, s0, state_transition_matrix,
                                          behavior_action_matrix,
                                          transition_based_rewards,
                                          horizon, gamma=1){
  aperm(replicate(n, generate_discrete_MDP_trajectory_with_pi_a_pi_b(s0, state_transition_matrix,
                                                                     behavior_action_matrix,
                                                                     transition_based_rewards,
                                                                     horizon, gamma)), c(3,1,2))
}

# Monte-Carlo estimation of the value of the evaluation policy.
# This version uses IS, which is suboptimal in terms of variance if
# we just care about finding the true value of the policy
MC_step_IPS <- function(state_transition_matrix,
                        behavior_action_matrix,
                        transition_based_rewards, 
                        evaluation_action_matrix,
                        horizon, M){
  V_step_IS <- c()
  for(m in 1:M){
    H <- generate_discrete_MDP_trajectory(1, state_transition_matrix,
                                          behavior_action_matrix,
                                          transition_based_rewards,
                                          horizon)
    pi_b <- apply(H, 1, function(x) behavior_action_matrix[x[1], x[2]]) # propensity score 
    pi_a <- apply(H, 1, function(x) evaluation_action_matrix[x[1], x[2]]) # probabilities of actions taken in H under the evaluation policy
    V_step_IS <- c(V_step_IS,
                   sum(H[, 'r'] * cumprod(pi_a / pi_b)))
    
  }
  mean(V_step_IS)
}
# Monte Carlo estimation of the value of the evaluation policy
# This version samples directly from the evaluation policy
# More efficient than the above version that uses IS
MC_direct_evaluation <- function(state_transition_matrix,
                                 transition_based_rewards, 
                                 evaluation_action_matrix,
                                 horizon, M){
  total_rewards <- c()
  for(m in 1:M){
    H <- generate_discrete_MDP_trajectory(1, state_transition_matrix,
                                          behavior_action_matrix=evaluation_action_matrix,
                                          transition_based_rewards=transition_based_rewards,
                                          horizon=horizon)
    total_rewards <- c(total_rewards, sum(H[, 'r']))
    
  }
  mean(total_rewards)
}

# GridWorld environment specification 
# Define MDP transition matrices
# State transition matrix dimensions: current state x action x next state
state_transition_matrix <- array(data=NA, dim=c(16, 5, 16))
#s1;
state_transition_matrix[1, ,] <- rbind(c(0,1/2,0,0,1/2,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
#s2;
state_transition_matrix[2, ,] <- rbind(c(1/3,0,1/3,0,0,1/3,0,0,0,0,0,0,0,0,0,0),
                                       c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0))
#s3;
state_transition_matrix[3, ,] <- rbind(c(0,1/3,0,1/3,0,0,1/3,0,0,0,0,0,0,0,0,0),
                                       c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0))
#s4;
state_transition_matrix[4, ,] <- rbind(c(0,0,1/2,0,0,0,0,1/2,0,0,0,0,0,0,0,0),
                                       c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0),
                                       c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0),
                                       c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0),
                                       c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0))
#s5;
state_transition_matrix[5, ,] <- rbind(c(1/3,0,0,0,0,1/3,0,0,1/3,0,0,0,0,0,0,0),
                                       c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
#s6;
state_transition_matrix[6, ,] <- rbind(c(0,1/4,0,0,1/4,0,1/4,0,0,1/4,0,0,0,0,0,0),
                                       c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
#s7;
state_transition_matrix[7, ,] <- rbind(c(0,0,1/4,0,0,1/4,0,1/4,0,0,1/4,0,0,0,0,0),
                                       c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
#s8;
state_transition_matrix[8, ,] <- rbind(c(0,0,0,1/3,0,0,1/3,0,0,0,0,1/3,0,0,0,0),
                                       c(0,0,0,0,0,0,0,0.993,0,0,0,0.007,0,0,0,0),
                                       c(0,0,0,0,0,0,0,0.007,0,0,0,0.993,0,0,0,0),
                                       c(0,0,0,0,0,0,0,0.998,0,0,0,0.002,0,0,0,0),
                                       c(0,0,0,0,0,0,0,0.002,0,0,0,0.998,0,0,0,0))
#s9;
state_transition_matrix[9, ,] <- rbind(c(0,0,0,0,1/3,0,0,0,0,1/3,0,0,1/3,0,0,0),
                                       c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                       c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
#s10;
state_transition_matrix[10, ,] <- rbind(c(0,0,0,0,0,1/4,0,0,1/4,0,1/4,0,0,1/4,0,0),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
#s11;
state_transition_matrix[11, ,] <- rbind(c(0,0,0,0,0,0,1/4,0,0,1/4,0,1/4,0,0,1/4,0),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
#s12;
state_transition_matrix[12, ,] <- rbind(c(0,0,0,0,0,0,0,1/3,0,0,1/3,0,0,0,0,1/3),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1))
#s13;
state_transition_matrix[13, ,] <- rbind(c(0,0,0,0,0,0,0,0,1/2,0,0,0,0,1/2,0,0),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
#s14;
state_transition_matrix[14, ,] <- rbind(c(0,0,0,0,0,0,0,0,0,1/3,0,0,1/3,0,1/3,0),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
#s15;
state_transition_matrix[15, ,] <- rbind(c(0,0,0,0,0,0,0,0,0,0,1/3,0,0,1/3,0,1/3),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
#s16;
state_transition_matrix[16, ,] <- rbind(c(0,0,0,0,0,0,0,0,0,0,0,1/2,0,0,1/2,0),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
                                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1))

# In GridWorld, the reward is a deterministic function of the transition.
# At entry s, new_s of this natrix is the reward correspoding to s -> new_s
transition_based_rewards <- rbind(c(-1, -1, -1, -1, -1, -10, -1, 1, -1, -1, -1, -1, -1, -1, -1, 10),
                                  c(-1, -1, -1, -1, -1, -10, -1, 1, -1, -1, -1, -1, -1, -1, -1, 10),
                                  c(-1, -1, -1, -1, -1, -10, -1, 1, -1, -1, -1, -1, -1, -1, -1, 10),
                                  c(-1, -1, -1, -1, -1, -10, -1, 1, -1, -1, -1, -1, -1, -1, -1, 10),
                                  c(-1, -1, -1, -1, -1, -10, -1, 1, -1, -1, -1, -1, -1, -1, -1, 10),
                                  c(-1, -1, -1, -1, -1, -10, -1, 1, -1, -1, -1, -1, -1, -1, -1, 10),
                                  c(-1, -1, -1, -1, -1, -10, -1, 1, -1, -1, -1, -1, -1, -1, -1, 10),
                                  c(-1, -1, -1, -1, -1, -10, -1, 1, -1, -1, -1, -1, -1, -1, -1, 10),
                                  c(-1, -1, -1, -1, -1, -10, -1, 1, -1, -1, -1, -1, -1, -1, -1, 10),
                                  c(-1, -1, -1, -1, -1, -10, -1, 1, -1, -1, -1, -1, -1, -1, -1, 10),
                                  c(-1, -1, -1, -1, -1, -10, -1, 1, -1, -1, -1, -1, -1, -1, -1, 10),
                                  c(-1, -1, -1, -1, -1, -10, -1, 1, -1, -1, -1, -1, -1, -1, -1, 10),
                                  c(-1, -1, -1, -1, -1, -10, -1, 1, -1, -1, -1, -1, -1, -1, -1, 10),
                                  c(-1, -1, -1, -1, -1, -10, -1, 1, -1, -1, -1, -1, -1, -1, -1, 10),
                                  c(-1, -1, -1, -1, -1, -10, -1, 1, -1, -1, -1, -1, -1, -1, -1, 10),
                                  c(-1, -1, -1, -1, -1, -10, -1, 1, -1, -1, -1, -1, -1, -1, -1, 10))

# Define behavior (logging) policy: s x a

#In the gridworld problem, the behavior policy randomly selects actions with equal probability regardless of which
#state it is in.  This results in an expected return of approximately −72.377.
#This policy selects each of the four actions with probability 0.25, regardless of the observation.
behavior_action_matrix<-  rbind(c(1, 0, 0, 0, 0),
                                c(1, 0, 0, 0, 0),
                                c(1, 0, 0, 0, 0),
                                c(1, 0, 0, 0, 0),
                                c(1, 0, 0, 0, 0),
                                c(1, 0, 0, 0, 0),
                                c(1, 0, 0, 0, 0),
                                c(1, 0, 0, 0, 0),
                                c(1, 0, 0, 0, 0),
                                c(1, 0, 0, 0, 0),
                                c(1, 0, 0, 0, 0),
                                c(1, 0, 0, 0, 0),
                                c(1, 0, 0, 0, 0),
                                c(1, 0, 0, 0, 0),
                                c(1, 0, 0, 0, 0),
                                c(1, 0, 0, 0, 0))
#The evaluation policy is a significantly better policy, although it is still far from optimal (it learns
#to reach the goal while avoiding the large penalty, but it  does not remain in the position (2,4) very long.
evaluation_action_matrix_p4 <-  rbind(c(0, 0.5, 0.5, 0, 0),
                                      c(0, 0.5, 0.5, 0, 0),
                                      c(0, 0.5, 0.5, 0, 0),
                                      c(0, 0.5, 0.5, 0, 0),
                                      c(0, 0.5, 0.5, 0, 0),
                                      c(0, 0.5, 0.5, 0, 0),
                                      c(0, 0.5, 0.5, 0, 0),
                                      c(0, 1, 0, 0, 0),
                                      c(0, 0.5, 0.5, 0, 0),
                                      c(0, 0.5, 0.5, 0, 0),
                                      c(0, 0.5, 0.5, 0, 0),
                                      c(0, 0.5, 0.5, 0, 0),
                                      c(0, 0.5, 0.5, 0, 0),
                                      c(0, 0.5, 0.5, 0, 0),
                                      c(0, 0.5, 0.5, 0, 0),
                                      c(0, 0.5, 0.5, 0, 0))

#This policy is a hand coded near-optimal policy that moves quickly to the position s8 without hitting s6.  
#It then usually takes the action to move down, and occasionally takes the action to move right.  
#Once it is in s8, it moves almost deterministically to s16. (near-optimal, pi5) 
evaluation_action_matrix_p5 <-  rbind(c(0, 0, 0, 0.5, 0.5),
                                      c(0, 0, 0, 0.5, 0.5),
                                      c(0, 0, 0, 0.5, 0.5),
                                      c(0, 0, 0, 0.5, 0.5),
                                      c(0, 0, 0, 0.5, 0.5),
                                      c(0, 0, 0, 0.5, 0.5),
                                      c(0, 0, 0, 0.5, 0.5),
                                      c(0, 0, 0, 1, 0),
                                      c(0, 0, 0, 0.5, 0.5),
                                      c(0, 0, 0, 0.5, 0.5),
                                      c(0, 0, 0, 0.5, 0.5),
                                      c(0, 0, 0, 0.5, 0.5),
                                      c(0, 0, 0, 0.5, 0.5),
                                      c(0, 0, 0, 0.5, 0.5),
                                      c(0, 0, 0, 0.5, 0.5),
                                      c(0, 0, 0, 0.5, 0.5))

#Check V0 and Q0:
horizon <- 100

MCres_eval_p4 <- MC_direct_evaluation(state_transition_matrix, transition_based_rewards, 
                                   evaluation_action_matrix=evaluation_action_matrix_p4, horizon, M=10000)
MCres_eval_p5 <- MC_direct_evaluation(state_transition_matrix, transition_based_rewards, 
                                      evaluation_action_matrix=evaluation_action_matrix_p5, horizon, M=10000)
MCres_beha <- MC_direct_evaluation(state_transition_matrix, transition_based_rewards, 
                                   behavior_action_matrix, horizon, M=10000)


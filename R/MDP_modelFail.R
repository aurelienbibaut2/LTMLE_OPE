library(tensorA)
# Generate a discrete MDP trajectory
generate_discrete_MDP_trajectory <- function(s0, state_transition_matrix,
                                             behavior_action_matrix,
                                             transition_based_rewards,
                                             horizon){
  states_trajectory <- c(); actions_trajectory <- c(); rewards_trajectory <- c()
  s <- s0
  for(t in 1:horizon){
    a <- which(rmultinom(n=1, size=1, prob=behavior_action_matrix[s, ]) == 1)
    new_s <- which(rmultinom(n=1, size=1, prob=state_transition_matrix[s, a, ]) == 1)
    r <- transition_based_rewards[s, new_s]
    states_trajectory <- c(states_trajectory, s)
    actions_trajectory <- c(actions_trajectory, a)
    rewards_trajectory <- c(rewards_trajectory, r)
    s <- new_s
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

# Dynamic programming based computation of the value function
compute_true_V_and_Q <- function(state_transition_matrix,
                                 transition_based_rewards,
                                 evaluation_action_matrix, horizon, gamma=1){
  # True value-to-go under policy pi
  # Compute the true value-to-go under pi with dynamic programming
  state_transition_tensor <- as.tensor(state_transition_matrix, dims=c(s=4, a=2, new_s=4))
  rewards_tensor <- as.tensor(transition_based_rewards, dims=c(s=4, new_s=4))
  pi_tensor <- as.tensor(evaluation_action_matrix, dims=c(s=4, a=2))
  ER <- reorder(mul.tensor(state_transition_tensor, i='new_s', 
                           rewards_tensor, j='new_s', by=c('s', 'a')), i='s') # expected reward given s and a
  V0 <- matrix(NA, nrow=horizon, ncol=4) # matrix of values to go: s X t
  Q0 <- array(NA, dim=c(horizon, 4, 2))
  Q0[horizon, ,] <- ER
  V0[horizon, ] <- mul.tensor(ER, i='a', pi_tensor, j='a', by='s')
  for(t in (horizon-1):1){
    V0[t, ] = (mul.tensor(ER, i='a', pi_tensor, j='a', by='s')
               + gamma * mul.tensor(X=pi_tensor, i='a',
                                    Y=mul.tensor(state_transition_tensor, i='new_s', 
                                                 as.tensor(V0[t+1, ], dims=c('new_s'=4)), j='new_s',
                                                 by=c('s', 'a')),
                                    j='a',
                                    by='s')
    )
    Q0[t, , ] <- (ER + gamma * mul.tensor(state_transition_tensor, i='new_s',
                                          as.tensor(V0[t+1, ], dims=c('new_s'=4)), j='new_s',
                                          by=c('s', 'a')))
  }
  list(V0=V0, Q0=Q0)
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

# ModelFail environment specification (see page 23 of Thomas and Brunskill)
# Define MDP transition matrices
# State transition matrix dimensions: current state x action x next state
state_transition_matrix <- array(data=NA, dim=c(4, 2, 4))
state_transition_matrix[, , 1] <- cbind(c(0, 0, 0, 0),
                                        c(0, 0, 0, 0))
state_transition_matrix[, , 2] <- cbind(c(1, 0, 0, 0),
                                        c(0, 0, 0, 0))
state_transition_matrix[, , 3] <- cbind(c(0, 0, 0, 0),
                                        c(1, 0, 0, 0))
state_transition_matrix[, , 4] <- cbind(c(0, 1, 0, 1),
                                        c(0, 0, 1, 1))
# Rewards:
transition_based_rewards <- cbind(c(0, 0, 0, 0),
                                  c(0, 0, 0, 0),
                                  c(0, 0, 0, 0),
                                  c(0, 1, -1, 0))
# Define behavior (logging) policy: s x a
behavior_action_matrix <- rbind(c(0.88, 0.12),
                                c(1, 0),
                                c(0, 1),
                                c(0.5, 0.5))
evaluation_action_matrix <-  rbind(c(0.12, 0.88),
                                   c(1, 0),
                                   c(0, 1),
                                   c(0.5, 0.5))

#Check V0 and Q0:
horizon <- 5
V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
                                  transition_based_rewards,
                                  evaluation_action_matrix, 
                                  horizon)
V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0

MCres_eval <- MC_direct_evaluation(state_transition_matrix, transition_based_rewards, 
                                   evaluation_action_matrix, horizon, M=1000)
MCres_beha <- MC_direct_evaluation(state_transition_matrix, transition_based_rewards, 
                                   behavior_action_matrix, horizon, M=1000)
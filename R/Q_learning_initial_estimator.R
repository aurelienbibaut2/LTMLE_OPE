source('MDP_modelWin.R')
source('Estimators.R')

set.seed(1)
horizon <- 10; n <- 1000
V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
                                  transition_based_rewards,
                                  evaluation_action_matrix, horizon, gamma=0.8)
V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0

Q_hat <- Q0; V_hat <- V0
D <- generate_discrete_MDP_dataset(n, 1, state_transition_matrix,
                                   behavior_action_matrix,
                                   transition_based_rewards,
                                   horizon)

# Make dataset of transitions with next state and done indicator
# apply(D[1:3, , ], 1, function(d) cbind(d[1:5, c('s', 'a', 'r')], new_s = c(d[2:5, c('s')], NA))
add_new_s <- function(d, horizon){
  cbind( d[1:horizon, c('s', 'a', 'r')], new_s = c(d[2:horizon, 's'], NA) , done = c(rep(0, horizon-1), 1) )
}

transitions_dataset <- do.call('rbind', lapply(1:dim(D)[1], function(i) add_new_s(D[i, , ], horizon=dim(D)[2])))


# Initialize Q_hat and other things to be initialized
bellman_dataset <- cbind(transitions_dataset, bellman_target=0)
n_states <- 3; n_actions <- 2; gamma <- 0.9
Q_hat <- array(0, dim=c(n_states, n_actions))
pi_tensor <- as.tensor(evaluation_action_matrix, dims=c(s=3, a=2))

# Define one helper function
na_to_zero <- Vectorize(function(x){
  if(is.na(x)){
    return(0)
  }else{
    return(x)
  }
})


L2_errors <- c()
# Bellman iterations
for(m in 1:10){
  # Compute Bellman target
  V_hat <- mul.tensor(as.tensor(Q_hat, dims=c(s=3, a=2)), i='a',
                      pi_tensor, j='a', by='s')
  # Communicate a bit
  print(as.vector(V_hat))
  print(as.vector(V0[1, ]))
  
  # MSE
  L2_errors <- c(L2_errors, sum( (as.vector(V_hat) - as.vector(V0[1, ]))^2 ) )
  
  # Regress Q_hat on Bellman target 
  bellman_dataset[, 'bellman_target'] <- apply(transitions_dataset, 1, function(x) x['r'] + gamma * na_to_zero(as.numeric(V_hat[x['new_s']])) )
  Q_hat_lm_fit <- lm(bellman_target ~ (factor(s) + factor(a))^2, data=data.frame(bellman_dataset))
  # Update Q_hat
  Q_hat <- outer(1:n_states, 1:n_actions, FUN= function(s, a) predict(Q_hat_lm_fit, newdata=data.frame(cbind(s=s, a=a))) )
  
}

plot(1:length(L2_errors), L2_errors)

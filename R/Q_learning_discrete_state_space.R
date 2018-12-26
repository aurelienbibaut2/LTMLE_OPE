# Make dataset of transitions with next state and done indicator
make_transitions_dataset <- function(D){
  add_new_s <- function(d, horizon){
    cbind( d[1:horizon, c('s', 'a', 'r')], new_s = c(d[2:horizon, 's'], NA) , done = c(rep(0, horizon-1), 1) )
  }
  transitions_dataset <- do.call('rbind', lapply(1:dim(D)[1], function(i) add_new_s(D[i, , ], horizon=dim(D)[2])))
}

# Define one helper function
na_to_zero <- Vectorize(function(x){
  if(is.na(x)){
    return(0)
  }else{
    return(x)
  }
})

# Bellman iterations
bellman_iterations <- function(transitions_dataset, evaluation_action_matrix, gamma, n_iterations, 
                               Q_hat_init=NULL){
  # Initialize Q_hat and other things to be initialized
  bellman_dataset <- cbind(transitions_dataset, bellman_target=0)
  n_states <- dim(evaluation_action_matrix)[1]; n_actions <- dim(evaluation_action_matrix)[2]
  
  if(!is.null(Q_hat_init)){
    Q_hat <- Q_hat_init
  }else{
    Q_hat <- array(0, dim=c(n_states, n_actions))
  }
  pi_tensor <- as.tensor(evaluation_action_matrix, dims=c(s=n_states, a=n_actions))
  
  # Bellman iterations
  for(m in 1:n_iterations){
    # Compute Bellman target
    V_hat <- mul.tensor(as.tensor(Q_hat, dims=c(s=n_states, a=n_actions)), i='a',
                        pi_tensor, j='a', by='s')
    # Regress Q_hat on Bellman target 
    bellman_dataset[, 'bellman_target'] <- apply(transitions_dataset, 1, 
                                                 function(x) x['r'] + gamma * na_to_zero(as.numeric(V_hat[x['new_s']])) )
    Q_hat_lm_fit <- lm(bellman_target ~ (factor(s) + factor(a))^2, data=data.frame(bellman_dataset))
    # Update Q_hat
    Q_hat <- outer(1:n_states, 1:n_actions, FUN= function(s, a) predict(Q_hat_lm_fit, newdata=data.frame(cbind(s=s, a=a))) )
  }
  return(Q_hat)
}
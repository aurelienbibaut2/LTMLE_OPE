library(gym)
library(keras)
remote_base <- "http://127.0.0.1:5000"
client <- create_GymClient(remote_base)
print(client)

# Create environment
env_id <- "CartPole-v0"
instance_id <- env_create(client, env_id)
print(instance_id)

# List all environments
all_envs <- env_list_all(client)
print(all_envs)


model_pi_b <- load_model_hdf5("cartPole_25_policy.h5")

# Behavior policy
pi_b <- function(state){
  #model_pi_b$predict(x=matrix(state, nrow=1))
  c(0.1, 0.9)
}

model_pi_e <- load_model_hdf5("~/cartPole_84_policy.h5")
# Target/evaluation policy
pi_e <- function(state){
  model_pi_e$predict(x=matrix(state, nrow=1))
}

generate_trajectory <- function(client, instance_id, pi_b, pi_e, max_steps=80){
  states <- c(); actions <- c(); rewards <- c()
  s <- unlist(env_reset(client, instance_id))
  done <- F
  step <- 0
  omega_t_vals <- c()
  while(!done & step < max_steps){
    a <- rbinom(1, 1, pi_b(s)[2])
    
    step_results <- env_step(client, instance_id, a, render = F)
    s2 <- unlist(step_results$observation)
    r <- step_results$reward
    done <- step_results$done
    
    omega_t <- (pi_e(s)/ pi_b(s))[a+1]
    
    states <- rbind(states, s)
    actions <- c(actions, a)
    rewards <- c(rewards, r)
    omega_t_vals <- c(omega_t_vals, omega_t)
    
    
    s <- s2
    step <- step + 1
    # print(step)
  }
  
  states <- rep(1, max_steps)
  actions <- actions + 1
  rho_t <- cumprod(omega_t_vals)
  rewards <- 2*(rewards - 1/2)
  
  if(step < max_steps){
    actions <- c(actions, rep(1, max_steps-step))
    rho_t <- c(rho_t, rep(rho_t[step], max_steps-step))
    rewards <- c(rewards, rep(-1, max_steps-step))
  }
  
  rtg <- rev(cumsum(rev(rewards)))
  cbind(s=states, a=actions, r=rewards, rtg=rtg, rho_t=rho_t, gamma_t = rep(1, max_steps))
}

generate_CartPole_dataset <- function(n, max_steps, client, instance_id, pi_b, pi_e){
  aperm(replicate(n, generate_trajectory(client, instance_id, pi_b, pi_e, max_steps)), c(3,1,2))
}


MC_true_value <- function(n, max_steps, pi_e){
  D <- generate_CartPole_dataset(n, max_steps, client, instance_id, pi_e, pi_e)
  estimate <- mean(D[, 1, 'rtg'])
  stddev <- sd(D[, 1, 'rtg'])
  list(estimate, LB=estimate - 1.96*stddev / sqrt(n), UB=estimate + 1.96*stddev / sqrt(n))
}

# Ridiculously biased initial estimators
get_Vhat_Qhat <- function(horizon){
  V0 <- matrix(0, nrow=horizon, ncol=1) # matrix of values to go: s X t
  Q0 <- array(0, dim=c(horizon, 1, 2))
  for(t in horizon:1){
    Q0[t, , ] <- horizon - t
    V0[t, ] <- horizon - t
  }
  list(V0=V0, Q0=Q0)
}

evaluation_action_matrix <- matrix(c(0.5, 0.5), nrow=1, ncol=2)
# Debugging experiment ----------------------------------------------------
source('Estimators.R')
source('partial_LTMLE.R')
horizon <- 40
Q_hat_and_V_hat <- get_Vhat_Qhat(horizon=horizon)
Q_hat <- Q_hat_and_V_hat$Q0
V_hat <- Q_hat_and_V_hat$V0

D <- generate_CartPole_dataset(n=50, max_steps=horizon, client, instance_id, pi_b, pi_e)
MC_res <- MC_true_value(200, horizon, pi_e)
V0 <- MC_res[[1]]
partial_LTMLE_estimator(D, Q_hat, V_hat, evaluation_action_matrix = evaluation_action_matrix, gamma=1, alpha=1, j=horizon)$estimate
WDR_estimator_TB_softened(D, Q_hat, V_hat, gamma=1, alpha=1, j=NULL)
D[1, , 'rho_t']

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

generate_trajectory <- function(client, instance_id, model){
  states <- c(); actions <- c(); rewards <- c()
  s <- unlist(env_reset(client, instance_id))
  done <- F
  while(!done){
    a <- rbinom(1, 1, model$predict(x=matrix(s, nrow=1))[2])
    
    step_results <- env_step(client, instance_id, a, render = F)
    s2 <- unlist(step_results$observation)
    r <- step_results$reward
    done <- step_results$done
    
    states <- rbind(states, s)
    actions <- c(actions, a)
    rewards <- c(rewards, r)
    
    s <- s2
  }
  rtg <- rev(cumsum(rev(rewards)))
  list(states=states, actions=to_categorical(actions), rtg=matrix(rtg))
}

# Define net
inputs <- layer_input(shape=4)
pi <- inputs
for(i in 1:2){
  pi <- layer_dense(units=32, activation="relu")(pi)
}
pi <- layer_dense(units=2, activation = 'softmax')(pi)

model <- keras_model(inputs=inputs, outputs=pi)

# Define loss
action_one_hot_ph <- k_placeholder(shape=list(NULL, 2))
action_prob <- model$output
rtg_ph <- k_placeholder(shape=list(NULL))
loss <- -k_mean(k_log(k_sum(action_prob * action_one_hot_ph, axis=1)) * rtg_ph)

adam <- optimizer_adam(lr=0.0001)

updates <- adam$get_updates(params=model$trainable_weights,
                            loss=loss)

train_fn <- k_function(inputs=list(model$input,
                                   action_one_hot_ph,
                                   rtg_ph),
                       outputs=list(),
                       updates=updates)
stored_total_rewards <- c()
model_25_saved <- F; model_30_saved <- F; model_35_saved <- F
for(i in 1:1000){
  H <- generate_trajectory(client, instance_id, model)
  train_fn(list(H$states, H$actions, H$rtg) )
  stored_total_rewards <- c(stored_total_rewards, H$rtg[1])
  if(i %% 10 == 0){
    print(mean(stored_total_rewards))
    
    if((!model_25_saved) && mean(stored_total_rewards) > 25){
      save_model_hdf5(model, filepath='cartPole_25_policy.h5', 
                      overwrite = TRUE,
                      include_optimizer = TRUE)
      model_25_saved <- T
    }
    if((!model_30_saved) && mean(stored_total_rewards) > 30){
      save_model_hdf5(model, filepath='cartPole_30_policy.h5', 
                      overwrite = TRUE,
                      include_optimizer = TRUE)
      model_30_saved <- T
    }
    if((!model_35_saved) && mean(stored_total_rewards) > 35){
      save_model_hdf5(model, filepath='cartPole_35_policy.h5', 
                      overwrite = TRUE,
                      include_optimizer = TRUE)
      model_35_saved <- T
      break
    }
    stored_total_rewards <- c()
  }
}
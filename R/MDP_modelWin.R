source(here("R/mdp_utils.R"))

# ModelWin environment specification (see page 24 of Thomas and Brunskill)
# Define MDP transition matrices
# State transition matrix dimensions: current state x action x next state
state_transition_matrix <- array(data=NA, dim=c(3, 2, 3))
state_transition_matrix[1, ,] <- rbind(c(0, 0.4, 0.6),
                                       c(0, 0.6, 0.4))
state_transition_matrix[2, ,] <- rbind(c(1, 0, 0),
                                       c(1, 0, 0))
state_transition_matrix[3, ,] <- rbind(c(1, 0, 0),
                                       c(1, 0, 0))
# In ModelWin, the reward is a deterministic function of the transition.
# At entry s, new_s of this natrix is the reward correspoding to s -> new_s
transition_based_rewards <- rbind(c(0, 1, -1),
                                  c(0, 0, 0),
                                  c(0, 0, 0))
# Define behavior (logging) policy: s x a
behavior_action_matrix <- rbind(c(0.73, 0.27),
                                c(0.5, 0.5),
                                c(0.5, 0.5))
evaluation_action_matrix <-  rbind(c(0.27, 0.73),
                                   c(0.5, 0.5),
                                   c(0.5, 0.5))

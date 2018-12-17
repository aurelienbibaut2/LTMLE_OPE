source(here("R/mdp_utils.R"))

# ModelFails environment specification (see page 24 of Thomas and Brunskill)
# Define MDP transition matrices
# State transition matrix dimensions: current state x action x next state
state_transition_matrix <- array(data=NA, dim=c(4, 2, 4))
state_transition_matrix[1, ,] <- rbind(c(0, 0.5, 0.5, 0),
                                       c(0, 0.5, 0.5, 0))
state_transition_matrix[2, ,] <- rbind(c(0, 0, 0, 1),
                                       c(0, 0, 0, 1))
state_transition_matrix[3, ,] <- rbind(c(0, 0, 0, 1),
                                       c(0, 0, 0, 1))
state_transition_matrix[4, ,] <- rbind(c(0, 0, 0, 0),
                                       c(0, 0, 0, 0))
# At entry s, new_s of this natrix is the reward correspoding to s -> new_s
transition_based_rewards <- rbind(c(0, 0, 0, 0),
                                  c(0, 0, 0, 1),
                                  c(0, 0, 0, -1),
                                  c(0, 0, 0, 0))
# Define behavior (logging) policy: s x a
behavior_action_matrix <- rbind(c(0.88, 0.12),
                                c(0.5, 0.5),
                                c(0.5, 0.5),
                                c(0.5, 0.5))
evaluation_action_matrix <-  rbind(c(0.12, 0.88),
                                   c(0.5, 0.5),
                                   c(0.5, 0.5),
                                   c(0.5, 0.5))

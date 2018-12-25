import numpy as np
"""The modelWin environment is a discrete MDP where r_t is a
deterministic function of (s_t, s_{t+1}). This is to be used
in conjuction with the discrete_MDP_transition_based_rewards class."""

# Transition matrix of the states: state * action * new_state
state_transition_matrix = np.array([[[0, 0.4, 0.6], [0., 0.6, 0.4]],
                                    [[1, 0, 0], [1, 0, 0]],
                                    [[1, 0, 0], [1, 0, 0]]])
# Transition matrix of giving the deterministic rewards
# given the transition (s, new_s)
# Dimensions: state * new_s
transition_based_rewards = np.array([[0, 1, -1],
                                     [0, 0, 0],
                                     [0, 0, 0]])
# Behavior policy is specified by its action matrix, that
# is P(a|s). Dimensions: s * a
behavior_action_matrix = np.array([[0.73, 0.27],
                                   [0.5, 0.5],
                                   [0.5, 0.5]])
# Same for evaluation policy
evaluation_policy_matrix = np.array([[0.27, 0.73],
                                     [0.5, 0.5],
                                     [0.5, 0.5]])
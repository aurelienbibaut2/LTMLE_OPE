import numpy as np
from numpy.random import multinomial

class discrete_MDP_transition_based_rewards(object):
    """ MDP class for discrete states, actions and rewards,
    where the rewards at t is a deterministic function of the
    transition (s_t, s_{t+1})."""
    def __init__(self, state_transition_matrix, transition_based_rewards):
        self.state_transition_matrix = state_transition_matrix
        self.transition_based_rewards = transition_based_rewards

    def reset(self, initial_state):
        self.state = initial_state

    def step(self, action):
        vector_state = multinomial(n=1, pvals=self.state_transition_matrix[self.state, action, :], size=1)
        new_state = np.sum(vector_state * np.arange(0, self.state_transition_matrix.shape[0]))
        reward = self.transition_based_rewards[self.state, new_state]
        self.state = new_state
        return reward

    def get_true_Q_and_V_functions(self, action_matrix, gamma, horizon):
        """Use dynamic programming to compute the true value and action-value
        function under a policy, from the transition matrices of the MDP and
        of the policy"""
        n_states = self.state_transition_matrix.shape[0]
        n_actions = self.state_transition_matrix.shape[1]

        V0 = np.zeros((horizon, n_states))
        Q0 = np.zeros((horizon, n_states, n_actions))
        # In einsum, z will stand for the new state
        ER = np.einsum('saz,sz->sa', self.state_transition_matrix, self.transition_based_rewards)
        Q0[horizon - 1, :, :] = ER
        V0[horizon - 1, :] = np.einsum('sa,sa->s', action_matrix, Q0[horizon - 1, :])
        for t in np.arange(0, horizon - 1)[::-1]:
            Q0[t, :] = ER + gamma * np.einsum('saz,z->sa', self.state_transition_matrix, V0[t + 1, :])
            V0[t, :] = np.einsum('sa,sa->s', action_matrix, Q0[t, :])

        return Q0, V0
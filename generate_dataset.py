import numpy as np
from numpy.random import multinomial

def generate_trajectory(env, s0, behavior_action_matrix, horizon):
    env.reset(s0)
    sar_triples = [] # State, action, reward triples
    for t in range(horizon):
        triple = {}
        s = env.state
        # Sample action
        a = np.sum(np.arange(behavior_action_matrix.shape[1]) *
                   multinomial(n=1, pvals=behavior_action_matrix[env.state, :], size=1)[0])
        r = env.step(a)
        sar_triples.append([s, a, r])

    return np.array(sar_triples)

def generate_dataset(env, s0, horizon, n,
                    behavior_action_matrix, evaluation_policy_matrix):
    trajectories = []
    for i in range(n):
        trajectory =  generate_trajectory(env, s0, behavior_action_matrix, horizon)
        pi_b = np.apply_along_axis(lambda x: behavior_action_matrix[x[0], x[1]], axis=1, arr=trajectory)
        pi_e = np.apply_along_axis(lambda x: evaluation_policy_matrix[x[0], x[1]], axis=1, arr=trajectory)
        rho_t = np.cumprod(pi_e / pi_b)
        trajectories.append(np.hstack((trajectory,
                                       pi_b[:, np.newaxis],
                                       pi_e[:, np.newaxis],
                                       rho_t[:, np.newaxis])))
    # Dataset has dimension n x horizon x sar triple info
    return np.array(trajectories)
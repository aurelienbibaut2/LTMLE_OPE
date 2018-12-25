import numpy as np

def IPS(D):
    horizon = D.shape[1]
    return np.mean(D[:, horizon-1, 5] * np.sum(D[:, :, 2], axis=1))

def stepIPS(D):
    return np.mean(np.sum(D[:, :, 5] * D[:, :, 2], axis=1) )

def WIPS(D):
    horizon = D.shape[1]
    w_H = np.mean(D[:, horizon-1, 5])
    return np.mean(D[:, horizon-1, 5] / w_H * np.sum(D[:, :, 2], axis=1))

def stepWIPS(D):
    n = D.shape[0]
    horizon = D.shape[1]
    w_t = np.mean(D[:, :, 5], axis=0)
    return np.mean(np.sum(D[:, :, 5] / (np.ones((n,1)) * w_t.T) * D[:, :, 2], axis=1))

def DR(D, Q_hat, V_hat):
    #import pdb; pdb.set_trace()
    n = D.shape[0]
    horizon = D.shape[1]
    D_star = np.zeros((n, horizon))

    t = horizon-1
    D_star[:, t] = D[:, t, 5] * (D[:, t, 2] - Q_hat[t, D[:, t, 0].astype(int), D[:, t, 1].astype(int)])
    for t in np.arange(0, horizon-1)[::-1]:
        D_star[:, t] = D[:, t, 5] * (D[:, t, 2] + V_hat[t+1, D[:, t+1, 0].astype(int)]
                                      - Q_hat[t, D[:, t, 0].astype(int), D[:, t, 1].astype(int)])
    return np.mean(V_hat[t, D[:, t, 0].astype(int)] + np.sum(D_star, axis=1))


def WDR(D, Q_hat, V_hat):
    n = D.shape[0]
    horizon = D.shape[1]
    D_star = np.zeros((n, horizon))
    w_t = np.mean(D[:, :, 5], axis=0)

    t = horizon - 1
    D_star[:, t] = D[:, t, 5] / w_t[t] * (D[:, t, 2] - Q_hat[t, D[:, t, 0].astype(int), D[:, t, 1].astype(int)])
    for t in np.arange(0, horizon - 1)[::-1]:
        D_star[:, t] = D[:, t, 5] / w_t[t] * (D[:, t, 2] + V_hat[t + 1, D[:, t + 1, 0].astype(int)]
                                              - Q_hat[t, D[:, t, 0].astype(int), D[:, t, 1].astype(int)])
    return np.mean(V_hat[t, D[:, t, 0].astype(int)] + np.sum(D_star, axis=1))


# Define the helper functions for LTMLE
def logit(x):
    return np.log(x) - np.log(1 - x)
def expit(x):
    return 1 / (1 + np.exp(-x))

# Function to fit the targeting step of the logistic regression
def fit_epsilon(y, Q, weights, epsilon0, max_it):
    """Fit epsilon in the targeting step, using IRLS.
    The argument Q here is the \tilde{Q}(A^(i)_t, S^(i)_t)
    from the write-up"""
    epsilon = epsilon0
    for it in range(max_it):
        ll = -np.sum(weights * (y * np.log(expit(logit(Q) + epsilon))
                                + (1 - y) * np.log(1 - expit(logit(Q) + epsilon))))
        grad = -np.sum(weights * (y - expit(logit(Q) + epsilon)))
        hessian = np.sum(weights * expit(logit(Q) + epsilon)
                         * (1 - expit(logit(Q) + epsilon)))
        epsilon = epsilon - grad / hessian

    return epsilon

def LTMLE(D, Q_hat, V_hat, evaluation_policy_matrix):
    horizon = D.shape[1]; n = D.shape[0]

    V_evaluated = np.zeros(n)

    epsilons = []
    for t in np.arange(0, horizon)[::-1]:
        Delta_t = horizon - t
        R = D[:, t, 2]
        U_tilde = (R + V_evaluated + Delta_t) / (2 * Delta_t)
        Q_t_evaluated = Q_hat[t, D[:, t, 0].astype(int), D[:, t, 1].astype(int)]
        Q_tilde_evaluated = (Q_t_evaluated + Delta_t) / (2 * Delta_t)
        # Targeting step
        epsilons.append(fit_epsilon(y=U_tilde, Q=Q_tilde_evaluated, weights=D[:, t, 5],
                                   epsilon0=0, max_it=5))
        # Evaluate Q_tilde(s_t, a) for a_t = 1, a_t = 2
        Q_tilde_t_star = expit(logit((Q_hat[t, :, :] + Delta_t) / (2 * Delta_t)) + epsilons[-1])
        # Then set V_tilde(s_t) = \sum_{a} Q_tilde(s_t, a) pi_a(a|s_t)
        V_tilde = np.einsum('sa,sa->s', Q_tilde_t_star, evaluation_policy_matrix)
        # Compute V = 2 * Delta_t * (V_tilde - 1)
        V = 2 * Delta_t * (V_tilde - 1/2)
        # Evaluate V
        V_evaluated = V[D[:, t, 0].astype(int)]
    return np.mean(V_evaluated)
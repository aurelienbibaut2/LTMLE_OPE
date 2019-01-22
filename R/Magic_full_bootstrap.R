library(boot)
source('Estimators.R')
source('utils.R')
source('partial_LTMLE.R')
source('Magic_estimator.R')

MAGIC_new_bootstrap_LTMLE <- function(D, Q_hat, V_hat, gamma, evaluation_action_matrix, force_PD=T){
  n <- dim(D)[1]
  horizon <- dim(D)[2]
  
  # # Get g^(j)s
  n_ids <- 10
  alphas <- c(rep(0, n_ids/2), seq(0, 1, length.out = n_ids/2))
  lambdas <- c(rev(seq(0, 5e-4, length.out = n_ids/2)), rep(0, n_ids/2))
  js <- c(seq(1, horizon, length.out = n_ids/2), rep(horizon, n_ids/2))
  
  R1 <- min(horizon * 8, 40)
  R2 <- min(horizon * 4, 20)
  g_js <- sapply(1:n_ids, function(j) partial_LTMLE_estimator(D, Q_hat=Q_hat, V_hat=V_hat, gamma=gamma,
                                                              evaluation_action_matrix = evaluation_action_matrix, 
                                                              alpha=alphas[j], j=js[j], lambda=lambdas[j])$estimate)
  
  # # Get bias by bootstrapping g^(horizon)
  bootstrap_CI <- quantile(boot(data=D, statistic=function(data, indices) partial_LTMLE_estimator(D[indices, , ], Q_hat=Q_hat, V_hat=V_hat, gamma=gamma, 
                                                                                                  evaluation_action_matrix = evaluation_action_matrix, 
                                                                                                  alpha=1, j=horizon)$estimate,
                                R = R1)$t, probs = c(0.1 / 2, 1 - 0.1 / 2))
  b_n <- sapply(g_js, Vectorize(function(g_j) distance_to_interval(bootstrap_CI, g_j)) )
  
  # Get covariance matrix by bootstrapping all g_js
  X <- boot(data=D, 
            statistic=function(data, indices) sapply(1:n_ids, 
                                                     function(j) partial_LTMLE_estimator(D[indices, , ], 
                                                                                         Q_hat=Q_hat, V_hat=V_hat, gamma=gamma, 
                                                                                         evaluation_action_matrix = evaluation_action_matrix, 
                                                                                         alpha=alphas[j], j=js[j], lambda=lambdas[j])$estimate),
            R = R2)$t
  Omega_n <- t((X - rep(1, R2) %*% t(apply(X, 2, mean)))) %*% (X - rep(1, R2) %*% t(apply(X, 2, mean))) / R2
  
  # Define and solve QP
  Dmat <- Omega_n + b_n %*% t(b_n)
  if(force_PD)
    Dmat <- Matrix::nearPD(Dmat, eig.tol=1e-10)$mat
  Amat <- t(rbind(rep(1, n_ids), diag(n_ids)) )
  dvec <- rep(0, n_ids)
  b0 <- c(1, rep(0, n_ids))
  x_star <- solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=b0, meq=1)$solution
  
  # Compute the MAGIC estimate as the weighted sum of the g^(j)'s, that is x_star^\top b_n[2:horizon]
  estimate <- t(x_star) %*% g_js
  # Output
  list(estimate=estimate, x_star=x_star, g_js=g_js, b_n=b_n, Omega_n=Omega_n, bootstrap_CI=bootstrap_CI)
}

# # MAGIC full library debugging experiments -------------------------------------
# source('MDP_modelWin.R')
# horizon <- 5; gamma <- 1; n_states <- 3; n_actions <- 2
# V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
#                                   transition_based_rewards,
#                                   evaluation_action_matrix, horizon, gamma = gamma)
# V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0
# Q_hat <- array(dim=dim(Q0)); V_hat <- array(dim=dim(V0))
# Q_hat <- Q0; V_hat <- V0
# b <- 5e-2 * rnorm(1)
# Delta_t <- 0
# n <- 200; gamma <- 1
# 
# cat(detectCores(), 'cores detected\n')
# cl <- makeCluster(getOption("cl.cores", detectCores()-1), outfile = '')
# registerDoParallel(cl)
# x_stars <- foreach(i=1:63, .combine = rbind,
#                    .packages = c('boot', 'quadprog'),
#                    .verbose = T, .inorder = T) %dopar% {
#                      D <- generate_discrete_MDP_dataset(n, 1, state_transition_matrix,
#                                                         behavior_action_matrix,
#                                                         transition_based_rewards,
#                                                         horizon)
#                      MAGIC_new_bootstrap_LTMLE(D, Q_hat, V_hat, gamma, evaluation_action_matrix)$x_star
#                    }
# stopCluster(cl)
# 
# mean_x_star <- apply(x_stars, 2, mean)
# x_star_plot <-ggplot(data=data.frame(id=1:10, mean_x_star=mean_x_star), aes(x=id, y=mean_x_star)) +
#   geom_bar(stat="identity")
# 
# library(gridExtra)
# grid.arrange(MSE_plot, x_star_plot, nrow=2)
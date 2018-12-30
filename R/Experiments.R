library(tensorA)
# library(here)
# source(here("R/MDP_modelWin.R"))
# source(here("R/Estimators.R"))
source('MDP_modelWin.R')
# source('MDP_modelFail.R')
source('Estimators.R')
source('C_TMLE.R')
source('Magic_estimator.R')
source('Q_learning_discrete_state_space.R')
source('penalized_LTMLE.R')

# Simulations -------------------------------------------------------------
horizon <- 20; gamma <- 0.95; n_states <- 3; n_actions <- 2
V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
                                  transition_based_rewards,
                                  evaluation_action_matrix, horizon, gamma = gamma)
V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0


# Specify jobs ------------------------------------------------------------
library(foreach); library(doParallel)
nb_repeats <- (parallel::detectCores() - 1) * 1
# ns <- c(50, 100, 200, 500, 1000, 5000, 10000)
ns <- c(100, 500, 1000)
jobs <- expand.grid(n = ns, repeat.id = 1:nb_repeats)


# Q-learning of the initial estimator -------------------------------------
# Generate data
# Q_hats <- list(); V_hats <- list()
# for(n in ns){
#   D <- generate_discrete_MDP_dataset(n, 1, state_transition_matrix,
#                                      behavior_action_matrix,
#                                      transition_based_rewards,
#                                      horizon)
#   # Fit Q_function
#   cat('Start Bellman iterations\n')
#   transitions_dataset <- make_transitions_dataset(D)
#   Q_learning_results <- bellman_iterations(transitions_dataset, evaluation_action_matrix,
#                                            gamma=gamma, max_it=100, relative_tol=1e-4, V0=V0[1, ],
#                                            verbose=T, start_new_plot = (n==ns[1]) )
#   cat('Done with Bellman iterations\n')
#   # Replicate accross time points the Bellman iterations based Q-function
#   # That would be correct under infinite horizon.
#   # Here's it's only approximately correct. Better for early time points than late time points.
#   
#   Q_hat <- array(0, dim=c(horizon, n_states, n_actions)); V_hat <- array(0, dim=c(horizon, n_states))
#   for(t in 1:horizon){
#     Q_hat[t, ,] <-  Q_learning_results$Q_hat
#     V_hat[t, ] <-  Q_learning_results$V_hat
#   }
# 
#   Q_hats[[as.character(n)]] <- Q_hat; V_hats[[as.character(n)]] <- V_hat
# }


# Monte Carlo simulation --------------------------------------------------
cat(detectCores(), 'cores detected\n')
cl <- makeCluster(getOption("cl.cores", detectCores()-1), outfile = '')
registerDoParallel(cl)

results <- foreach(i=1:nrow(jobs), .combine = rbind,
                   .packages = c('boot', 'quadprog'),
                   .verbose = T, .inorder = T) %dopar% {
                     D <- generate_discrete_MDP_dataset(jobs[i, ]$n, 1, state_transition_matrix,
                                                        behavior_action_matrix,
                                                        transition_based_rewards,
                                                        horizon)
                     b <- 3e-1 * rnorm(1)
                     Q_hat <- Q0 + b; V_hat <- V0 + b
                     #Q_hat <- Q_hats[[as.character(jobs[i, ]$n)]]
                     #V_hat <- V_hats[[as.character(jobs[i, ]$n)]]
                     
                     rbind(
                       #c(n=jobs[i, ]$n, estimator='IS', estimate=IS_estimator(D)),
                           # c(n=jobs[i, ]$n, estimator='WIS', estimate=WIS_estimator(D)),
                           # c(n=jobs[i, ]$n, estimator='stepIS', estimate=stepIS_estimator(D)),
                           # c(n=jobs[i, ]$n, estimator='stepWIS', estimate=stepWIS_estimator(D)),
                           # c(n=jobs[i, ]$n, estimator='DR',  estimate=DR_estimator_JL(D, Q_hat=Q0, V_hat=V0)),
                           # c(n=jobs[i, ]$n, estimator='DR_TB',  estimate=DR_estimator_TB(D, Q_hat=Q0, V_hat=V0)),
                           c(n=jobs[i, ]$n, estimator='WDR', estimate=try(WDR_estimator_TB(D, Q_hat=Q_hat, V_hat=V_hat,
                                                                                        gamma = gamma, j = horizon)$g_js[horizon+1])),
                           c(n=jobs[i, ]$n, estimator='MAGIC', estimate=MAGIC_estimator(D, Q_hat, V_hat, gamma = gamma,
                                                                                        horizon = horizon, n_bootstrap = 100,
                                                                                        force_PD = T)$estimate),
                           c(n=jobs[i, ]$n, estimator='LTMLE', estimate=try(LTMLE_estimator(D, Q_hat, V_hat, 
                                                                                               evaluation_action_matrix, gamma, alpha=1)$estimate)),
                           c(n=jobs[i, ]$n, estimator='LTMLE_0.7', estimate=try(LTMLE_estimator(D, Q_hat, V_hat, 
                                                                                                evaluation_action_matrix, gamma, alpha=0.7)$estimate)),
                           c(n=jobs[i, ]$n, estimator='LTMLE_0.1', estimate=try(LTMLE_estimator(D, Q_hat, V_hat, 
                                                                                       evaluation_action_matrix, gamma, alpha=0.1)$estimate)),
                           c(n=jobs[i, ]$n, estimator='penalizedLTMLE', estimate=try(penalized_LTMLE_estimator(D, Q_hat, V_hat, 
                                                                                                evaluation_action_matrix, gamma, alpha=0.1)$estimate)),
                           c(n=jobs[i, ]$n, estimator='C-TMLE', estimate=C_LTMLE_softening(D, Q_hat, V_hat, 
                                                                                           evaluation_action_matrix, gamma)$estimate)
                     )
                   }
stopCluster(cl)
# Compute MSE from results matrix
results_df <- transform(as.data.frame(results), 
                        n=as.numeric(as.character(n)),
                        estimate=as.numeric(as.character(estimate)))
results_df$squared_error <- (results_df$estimate - V0[1,1])^2

MSE_table <- aggregate(results_df$squared_error, 
                       list(estimator=results_df$estimator, n=results_df$n), 
                       mean)
colnames(MSE_table)[3] <- 'MSE'

# Plot nMSE against n
library(ggplot2)
MSE_plot <- ggplot(data=MSE_table, aes(x=log10(n), y=log10(n*MSE), color=estimator, shape=estimator, size=estimator)) + 
  # scale_shape_manual( values=c('LTMLE'=15, 'DR'=19, 'IS'=19, 'stepIS'=19, 'stepWIS'=19, 'WDR'=19, 'WIS'=19)) +
  # scale_size_manual( values=c('LTMLE'=6, 'DR'=2, 'IS'=2, 'stepIS'=2, 'stepWIS'=2, 'WDR'=2, 'WIS'=2)) +
  geom_line(size=1) + geom_point() + ggtitle(paste('ModelWin, horizon=', horizon, ', number of draws per point=', nb_repeats))
print(MSE_plot)
library(tensorA)
# library(here)
# source(here("R/MDP_modelWin.R"))
# source(here("R/Estimators.R"))
source('Estimators.R')
source('C_TMLE.R')
source('Magic_estimator.R')
source('Q_learning_discrete_state_space.R')
source('penalized_LTMLE.R')
source('MAGIC-LTMLE_Estimator.R')
source('single_epsilon_LTMLE.R')
source('partial_LTMLE.R')
source('MDP_modelWin.R')
# source('MDP_modelFail.R')

# Simulations -------------------------------------------------------------
# ModelFail parameters
# horizon <- 3; gamma <- 1; n_states <- 4; n_actions <- 2
# V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
#                                   transition_based_rewards,
#                                   evaluation_action_matrix, horizon, gamma = gamma)
# V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0
# Q_hat <- array(0.38, dim=c(3, 4, 2))
# V_hat <- array(0.38, dim=c(3, 4))
# Q_hat[3, , ] <- 0
# V_hat[3, ] <- 0

# ModelWin parameters
horizon <- 10; gamma <- 1; n_states <- 3; n_actions <- 2
V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
                                  transition_based_rewards,
                                  evaluation_action_matrix, horizon, gamma = gamma)
V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0

# Specify jobs ------------------------------------------------------------
library(foreach); library(doParallel)
nb_repeats <- (parallel::detectCores() - 1)  * 1
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
#   plot(Q_learning_results$l2_errors)
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
                     # Apply bias while respecting model's constraints
                     Q_hat <- array(dim=dim(Q0)); V_hat <- array(dim=dim(V0))
                     b <- 5e-2 * rnorm(1)
                     Delta_t <- 0
                     for(t in horizon:1){
                       Delta_t <- 1 + gamma * Delta_t
                       Q_hat[t, ,] <- 2 * Delta_t * (expit(logit( (Q0[t, ,] + Delta_t) / (2*Delta_t) ) + b) - 1/2)
                       V_hat[t, ] <- 2 * Delta_t * (expit(logit( (V0[t, ] + Delta_t) / (2 * Delta_t) ) + b) - 1/2)
                     }
                     #Q_hat <- Q_hats[[as.character(jobs[i, ]$n)]]
                     #V_hat <- V_hats[[as.character(jobs[i, ]$n)]]
                     
                     # D_large <- generate_discrete_MDP_dataset(1e4, 1, state_transition_matrix,
                     #                                          behavior_action_matrix,
                     #                                          transition_based_rewards,
                     #                                          horizon)
                     
                     C_LTMLE_result <- try(C_LTMLE_softening(D, Q_hat, V_hat,
                                                             evaluation_action_matrix,
                                                             gamma, D_large=NULL, V=3, greedy=T))
                     
                     # one_step_result <- try(one_step_LTMLE(D, Q_hat, V_hat, evaluation_action_matrix, 
                     #                                   gamma=gamma, verbose=F, max_it=20, step_size=1e-1))
                     LTMLE_1.0_result <- try(LTMLE_estimator(D, Q_hat, V_hat,
                                                             evaluation_action_matrix, gamma, alpha=1))
                     
                     
                     rbind(
                       #c(n=jobs[i, ]$n, estimator='IS', estimate=IS_estimator(D)),
                           # c(n=jobs[i, ]$n, estimator='WIS', estimate=WIS_estimator(D)),
                           # c(n=jobs[i, ]$n, estimator='stepIS', estimate=stepIS_estimator(D)),
                           # c(n=jobs[i, ]$n, estimator='stepWIS', estimate=stepWIS_estimator(D)),
                           # c(n=jobs[i, ]$n, estimator='DR',  estimate=DR_estimator_JL(D, Q_hat=Q0, V_hat=V0)),
                           # c(n=jobs[i, ]$n, estimator='DR_TB',  estimate=DR_estimator_TB(D, Q_hat=Q0, V_hat=V0)),
                           c(n=jobs[i, ]$n, estimator='WDR', estimate=try(WDR_estimator_TB(D, Q_hat=Q_hat, V_hat=V_hat,
                                                                                        gamma = gamma, j = horizon)$g_js[horizon+1]), 
                             base_est_id=NA, epsilon=NA, score_eq=NA)
                           , c(n=jobs[i, ]$n, estimator='MAGIC', estimate=try(MAGIC_estimator(D, Q_hat, V_hat, gamma = gamma,
                                                                                        horizon = horizon, n_bootstrap = 1000,
                                                                                        force_PD = T)$estimate),
                               base_est_id=NA, epsilon=NA, score_eq=NA)
                           # , c(n=jobs[i, ]$n, estimator='MAGIC-LTMLE', estimate=try(MAGIC_LTMLE_estimator_hacky(D, Q_hat, V_hat, evaluation_action_matrix, 
                           #                                                          gamma, n_bootstrap=1000) ))
                           , c(n=jobs[i, ]$n, estimator='LTMLE_1.0', estimate=try(LTMLE_1.0_result$estimate),
                               base_est_id=NA, epsilon=LTMLE_1.0_result$epsilons[1], score_eq=NA)
                           , c(n=jobs[i, ]$n, estimator='LTMLE_0.0', estimate=try(LTMLE_estimator(D, Q_hat, V_hat,
                                                                                                  evaluation_action_matrix, gamma, alpha=0)$estimate),
                               base_est_id=NA, epsilon=LTMLE_1.0_result$epsilons[1], score_eq=NA)
                           , c(n=jobs[i, ]$n, estimator='partial_LTMLE_1.0', estimate=try(partial_LTMLE_estimator(D, Q_hat, V_hat,
                                                                                                          evaluation_action_matrix, gamma, alpha=1,j=1)$estimate),
                               base_est_id=NA, epsilon=NA, score_eq=NA)
                           , c(n=jobs[i, ]$n, estimator='partial_LTMLE_0.3', estimate=try(partial_LTMLE_estimator(D, Q_hat, V_hat,
                                                                                                                  evaluation_action_matrix, gamma, alpha=0.3,j=1)$estimate),
                               base_est_id=NA, epsilon=NA, score_eq=NA)
                           # , c(n=jobs[i, ]$n, estimator='LTMLE_0.7', estimate=try(LTMLE_estimator(D, Q_hat, V_hat, 
                           #                                                                        evaluation_action_matrix, gamma, alpha=0.7)$estimate))
                           # , c(n=jobs[i, ]$n, estimator='LTMLE_0.5', estimate=try(LTMLE_estimator(D, Q_hat, V_hat, 
                           #                                                                       evaluation_action_matrix, gamma, alpha=0.5)$estimate))
                           # , c(n=jobs[i, ]$n, estimator='LTMLE_0.1', estimate=try(LTMLE_estimator(D, Q_hat, V_hat, 
                           #                                                             evaluation_action_matrix, gamma, alpha=0.1)$estimate))
                           # , c(n=jobs[i, ]$n, estimator='LTMLE_0.0', estimate=try(LTMLE_estimator(D, Q_hat, V_hat, 
                           #                                                                        evaluation_action_matrix, gamma, alpha=0)$estimate))
                           # #, c(n=jobs[i, ]$n, estimator='penalizedLTMLE', estimate=try(penalized_LTMLE_estimator(D, Q_hat, V_hat, 
                           #                                                                                      evaluation_action_matrix, gamma, alpha=1,
                           #                                                                                      penalty = 1)$estimate))
                           # # c(n=jobs[i, ]$n, estimator='penalizedLTMLE_3', estimate=try(penalized_LTMLE_estimator(D, Q_hat, V_hat, 
                           #                                                                                     # evaluation_action_matrix, gamma, alpha=1, 
                           #                                                                                     # penalty = 0.1)$estimate)),
                           # , c(n=jobs[i, ]$n, estimator='C-TMLE', estimate=try(C_LTMLE_softening(D, Q_hat, V_hat,
                           #                                                                       evaluation_action_matrix,
                           #                                                                       gamma, D_large=NULL, V=3, greedy=F)$estimate) )
                           , c(n=jobs[i, ]$n, estimator='MAGIC_LTMLE', estimate=try(MAGIC_bootstrap_LTMLE(D, Q_hat, V_hat, gamma,
                                                                                                          evaluation_action_matrix, force_PD=T)$estimate),
                               base_est_id=NA, epsilon=NA, score_eq=NA )
                           , c(n=jobs[i, ]$n, estimator='C-TMLE-sftning', estimate=try(C_LTMLE_result$estimate),
                               base_est_id=try(C_LTMLE_result$softening_coeff), epsilon=NA, score_eq=NA )
                           # , c(n=jobs[i, ]$n, estimator='1step_LTMLE', estimate=try(one_step_result$estimate),
                           #     base_est_id=NA, epsilon=one_step_result$epsilon, score_eq=try(one_step_result$score_eq))
                           #, c(n=jobs[i, ]$n, estimator='C-TMLE-pnlty', estimate=try(C_LTMLE_penalization(D, Q_hat, V_hat, 
                           #                                                                             evaluation_action_matrix,
                           #                                                                             gamma, V=3, plot_risk=F, greedy=T)$estimate) )
                           # c(n=jobs[i, ]$n, estimator='C-TMLE_true_risk', estimate=try(C_TMLE_results$true_risk_estimate))
                     )
                   }
stopCluster(cl)
# Compute MSE from results matrix
results_df <- transform(as.data.frame(results), 
                        n=as.numeric(as.character(n)),
                        estimate=as.numeric(as.character(estimate)))
estimators <- c('C-TMLE-sftning', 'MAGIC', 'MAGIC_LTMLE')
results_df <- subset(results_df, estimator %in% estimators)

results_df$squared_error <- (results_df$estimate - V0[1,1])^2

MSE_table <- aggregate(results_df$squared_error, 
                       list(estimator=results_df$estimator, n=results_df$n), 
                       mean)
colnames(MSE_table)[3] <- 'MSE'
var_table <- aggregate(results_df$estimate, 
                       list(estimator=results_df$estimator, n=results_df$n), 
                       var)
colnames(var_table)[3] <- 'var'
bias_table <- aggregate(results_df$estimate - V0[1,1], 
                        list(estimator=results_df$estimator, n=results_df$n), 
                        mean)
colnames(bias_table)[3] <- 'bias'

summary_table <- transform(cbind(MSE_table, var=var_table$var, bias=bias_table$bias),
                           MSE=round(MSE, 5), var=round(var, 5), bias=round(bias, 5))
print(summary_table)

# Base estimator id:
base_est_id_df <- subset(data.frame(results), subset=estimator=='C-TMLE-sftning' )
print(table(base_est_id_df$n, base_est_id_df$base_est_id))
# Plot nMSE against n
library(ggplot2)
MSE_plot <- ggplot(data=MSE_table, aes(x=log10(n), y=log10(n*MSE), color=estimator, shape=estimator)) + 
  scale_shape_manual( values=c('MAGIC'=15,
                               'C-TMLE-sftning'=15, '1step_LTMLE'=15, 'MAGIC_LTMLE'=15, 'partial_LTMLE_1.0'=15, 'partial_LTMLE_0.3'=15,
                               'LTMLE_1.0'=19, 'LTMLE_0.7'=19, 'LTMLE_0.5'=19, 'LTMLE_0.1'=19, 'LTMLE_0.0'=19, 
                               'WDR'=18) ) +
  scale_size_manual( values=c('MAGIC'=8,
                              'C-TMLE-sftning'=8, '1step_LTMLE'=8, 'MAGIC_LTMLE'=8, 'partial_LTMLE_1.0'=8, 'partial_LTMLE_0.3'=8,
                              'LTMLE_1.0'=4, 'LTMLE_0.7'=4, 'LTMLE_0.5'=4, 'LTMLE_0.1'=4, 'LTMLE_0.0'=4, 
                              'WDR'=4)) +
  geom_line(size=1) + geom_point(aes(size=estimator)) + ggtitle(paste('ModelWin, horizon=', horizon, ', number of draws per point=', nb_repeats))
print(MSE_plot)

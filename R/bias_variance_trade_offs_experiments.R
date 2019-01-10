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
horizon <- 20; gamma <- 1; n_states <- 3; n_actions <- 2
V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
                                  transition_based_rewards,
                                  evaluation_action_matrix, horizon, gamma = gamma)
V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0

# Specify jobs ------------------------------------------------------------
library(foreach); library(doParallel)
nb_repeats <- (parallel::detectCores() - 1) * 10
# ns <- c(50, 100, 200, 500, 1000, 5000, 10000)
ns <- c(1000)
n_ids <- 10
alphas <- seq(0, 1, length.out = n_ids)
js <- ceiling(seq(1, horizon, length.out=n_ids))
jobs <- expand.grid(n = ns, id = 1:n_ids,repeat.id = 1:nb_repeats)


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
                     b <- 2e-1 * rnorm(1)
                     Delta_t <- 0
                     for(t in horizon:1){
                       Delta_t <- 1 + gamma * Delta_t
                       Q_hat[t, ,] <- 2 * Delta_t * (expit(logit( (Q0[t, ,] + Delta_t) / (2*Delta_t) ) + b) - 1/2)
                       V_hat[t, ] <- 2 * Delta_t * (expit(logit( (V0[t, ] + Delta_t) / (2 * Delta_t) ) + b) - 1/2)
                     }
                     
                     
                     rbind(c(n=jobs[i, ]$n, id=jobs[i, ]$id, estimator='WDR',
                             estimate=WDR_estimator_TB_old(D, Q_hat, V_hat, gamma=1, j=js[jobs[i, ]$id])),
                           
                           c(n=jobs[i, ]$n, id=jobs[i, ]$id, estimator='partial_LTMLE',
                             estimate=partial_LTMLE_estimator(D, Q_hat, V_hat,
                                                              evaluation_action_matrix, gamma,
                                                              alpha=alphas[jobs[i, ]$id], js[jobs[i, ]$id])$estimate)
                           )
                   }
stopCluster(cl)
# Compute MSE from results matrix
results_df <- transform(as.data.frame(results), 
                        n=as.numeric(as.character(n)),
                        id=as.numeric(as.character(id)),
                        estimate=as.numeric(as.character(estimate)))

results_df$squared_error <- (results_df$estimate - V0[1,1])^2

MSE_table <- aggregate(results_df$squared_error, 
                       list(estimator=results_df$estimator, n=results_df$n, id=results_df$id), 
                       mean)
colnames(MSE_table)[4] <- 'MSE'
var_table <- aggregate(results_df$estimate, 
                       list(estimator=results_df$estimator, n=results_df$n, id=results_df$id), 
                       var)
colnames(var_table)[4] <- 'var'
bias_table <- aggregate(results_df$estimate - V0[1,1], 
                        list(estimator=results_df$estimator, n=results_df$n, id=results_df$id), 
                        mean)
colnames(bias_table)[4] <- 'bias'

summary_table <- transform(cbind(MSE_table, var=var_table$var, bias=bias_table$bias),
                           MSE=round(MSE, 5), var=round(var, 5), bias=round(bias, 5))
print(summary_table)

MSE_plot <- ggplot(data=MSE_table, aes(x=id, y=log10(n*MSE), color=estimator)) + geom_line() + geom_point()
print(MSE_plot)
# 
# # Base estimator id:
# base_est_id_df <- subset(data.frame(results), subset=estimator=='C-TMLE-sftning' )
# print(table(base_est_id_df$n, base_est_id_df$base_est_id))
# # Plot nMSE against n
# library(ggplot2)
# MSE_plot <- ggplot(data=MSE_table, aes(x=log10(n), y=log10(n*MSE), color=estimator, shape=estimator)) + 
#   scale_shape_manual( values=c('MAGIC'=15,
#                                'C-TMLE-sftning'=15, '1step_LTMLE'=15, 'MAGIC_LTMLE'=15, 'partial_LTMLE_1.0'=15, 'partial_LTMLE_0.3'=15,
#                                'LTMLE_1.0'=19, 'LTMLE_0.7'=19, 'LTMLE_0.5'=19, 'LTMLE_0.1'=19, 'LTMLE_0.0'=19, 
#                                'WDR'=18) ) +
#   scale_size_manual( values=c('MAGIC'=8,
#                               'C-TMLE-sftning'=8, '1step_LTMLE'=8, 'MAGIC_LTMLE'=8, 'partial_LTMLE_1.0'=8, 'partial_LTMLE_0.3'=8,
#                               'LTMLE_1.0'=4, 'LTMLE_0.7'=4, 'LTMLE_0.5'=4, 'LTMLE_0.1'=4, 'LTMLE_0.0'=4, 
#                               'WDR'=4)) +
#   geom_line(size=1) + geom_point(aes(size=estimator)) + ggtitle(paste('ModelWin, horizon=', horizon, ', number of draws per point=', nb_repeats))
# print(MSE_plot)

library(tensorA)
library(here)
source(here("R/MDP_modelWin.R"))
source(here("R/Estimators.R"))

# Simulations -------------------------------------------------------------
horizon <- 4
V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
                                  transition_based_rewards,
                                  evaluation_action_matrix, horizon)
V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0
# Specify jobs ------------------------------------------------------------
library(foreach); library(doParallel)
nb_repeats <- (parallel::detectCores() - 1) * 10
ns <- c(50, 100, 200, 500, 1000, 5000, 10000)
jobs <- expand.grid(n = ns, repeat.id = 1:nb_repeats)


cat(detectCores(), 'cores detected\n')
cl <- makeCluster(getOption("cl.cores", detectCores()-1), outfile = '')
registerDoParallel(cl)

results <- foreach(i=1:nrow(jobs), .combine = rbind,
                   #.packages = c('speedglm'),
                   #.export = c('TMLE_EY1_speedglm', 'CV_AIPTW_EY1', 'expit', 'logit', 'g_to_g_delta'),
                   .verbose = T, .inorder = T) %dopar% {
                     D <- generate_discrete_MDP_dataset(jobs[i, ]$n, 1, state_transition_matrix,
                                                        behavior_action_matrix,
                                                        transition_based_rewards,
                                                        horizon)
                     rbind(c(n=jobs[i, ]$n, estimator='IS', estimate=IS_estimator(D)),
                           c(n=jobs[i, ]$n, estimator='WIS', estimate=WIS_estimator(D)),
                           c(n=jobs[i, ]$n, estimator='stepIS', estimate=stepIS_estimator(D)),
                           c(n=jobs[i, ]$n, estimator='stepWIS', estimate=stepWIS_estimator(D)),
                           c(n=jobs[i, ]$n, estimator='DR',  estimate=DR_estimator_JL(D, Q_hat=Q0, V_hat=V0)),
                           # c(n=jobs[i, ]$n, estimator='DR_TB',  estimate=DR_estimator_TB(D, Q_hat=Q0, V_hat=V0)),
                           c(n=jobs[i, ]$n, estimator='WDR',  estimate=WDR_estimator_TB(D, Q_hat=Q0, V_hat=V0)),
                           c(n=jobs[i, ]$n, estimator='LTMLE', estimate=LTMLE_estimator(D, Q_hat=Q0, V_hat=V0, 
                                                                                        evaluation_action_matrix=evaluation_action_matrix))
                     )
                   }

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
  scale_shape_manual( values=c('LTMLE'=15, 'DR'=19, 'IS'=19, 'stepIS'=19, 'stepWIS'=19, 'WDR'=19, 'WIS'=19)) +
  scale_size_manual( values=c('LTMLE'=6, 'DR'=2, 'IS'=2, 'stepIS'=2, 'stepWIS'=2, 'WDR'=2, 'WIS'=2)) +
  geom_line(size=1) + geom_point() + ggtitle(paste('ModelWin, horizon=', horizon, ', number of draws per point=', nb_repeats))
print(MSE_plot)
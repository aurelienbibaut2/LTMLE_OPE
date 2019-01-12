source('Estimators.R')
source('penalized_LTMLE.R')
source('utils.R')
source('evaluate_EIC.R')
source('partial_LTMLE.R')

# Collaborative TMLE that uses a sequence of decreasingly softened weights
C_LTMLE_softening <- function(D, Q_hat, V_hat, evaluation_action_matrix, gamma, V=3, plot_risk=F, D_large=NULL, greedy=T){
  # Split the dataset. Compute sequence of epsilons for each softening coeff. Pick the one that minimizes a cross-validated risk
  # D_split1 <- D[1:floor(0.5 * n), ,]; D_split2 <- D[(floor(0.5 * n) + 1):n, ,]
  n <- dim(D)[1]
  softening_coeffs <- seq(1e-2, 1, length.out = 10)
  js <- (ceiling(seq(1, horizon, length.out=10))) 
  
  CV_doubly_roubust_risks <- rep(Inf, length(softening_coeffs))
  true_risks <- rep(0, length(softening_coeffs))
  for(i in 1:length(softening_coeffs)){
    CV_doubly_roubust_risks[i] <- 0
    for(v in 1:V){
      test_ids <- (floor( (v-1) / V * n) + 1) : floor( v / V * n)
      epsilons <- partial_LTMLE_estimator(D[setdiff(1:n, test_ids), ,], Q_hat, V_hat, evaluation_action_matrix, gamma, 
                                          alpha=softening_coeffs[i], j=js[i])$epsilons
      CV_doubly_roubust_risks[i] <- (CV_doubly_roubust_risks[i] +
                                       var(evaluate_EIC(D[test_ids, ,], epsilons, Q_hat, V_hat, evaluation_action_matrix, gamma))
                                     )
      if(!is.null(D_large)) true_risks[i] <- true_risks[i] + var(evaluate_EIC(D_large, epsilons, Q_hat, V_hat, evaluation_action_matrix, gamma))
      
    }
    if(greedy & (i > 1) && (CV_doubly_roubust_risks[i] > CV_doubly_roubust_risks[i-1]) ) break # Enforce greedy search
  }
  if(plot_risk){
    
    print(CV_doubly_roubust_risks)
    if(!is.null(D_large)){
      plot(softening_coeffs, CV_doubly_roubust_risks, ylim=range(c(true_risks, CV_doubly_roubust_risks)))
      points(softening_coeffs, true_risks, col='green')
      plot(softening_coeffs, true_risks, col='green')
    }else{
      plot(softening_coeffs, CV_doubly_roubust_risks)
    }
  }
  # Output the estimate computed on the full dataset that uses the softening coeff that minimizes the CV risk
  list(estimate=partial_LTMLE_estimator(D, Q_hat, V_hat, evaluation_action_matrix, gamma, 
                                alpha=softening_coeffs[which.min(CV_doubly_roubust_risks)], j=js[which.min(CV_doubly_roubust_risks)]
                                )$estimate,
       true_risk_estimate=partial_LTMLE_estimator(D, Q_hat, V_hat, evaluation_action_matrix, gamma, 
                  alpha=softening_coeffs[which.min(true_risks)], j=js[which.min(true_risks)]
                  )$estimate,
       softening_coeff=softening_coeffs[which.min(CV_doubly_roubust_risks)])
}


# Collaborative TMLE that uses a sequence of decreasingly softened weights
C_LTMLE_penalization <- function(D, Q_hat, V_hat, evaluation_action_matrix, gamma, V=3, plot_risk=F,  greedy=T){
  # Split the dataset. Compute sequence of epsilons for each softening coeff. Pick the one that minimizes a cross-validated risk
  n <- dim(D)[1]
  penalty_coeffs <- c(10^c(seq(0, -5, length.out = 9)), 0)
  CV_doubly_roubust_risks <- rep(Inf, length(penalty_coeffs))
  for(i in 1:length(penalty_coeffs)){
    CV_doubly_roubust_risks[i] <- 0
    for(v in 1:V){
      test_ids <- (floor( (v-1) / V * n) + 1) : floor( v / V * n)
      epsilons <- penalized_LTMLE_estimator(D[setdiff(1:n, test_ids), ,], Q_hat, V_hat, evaluation_action_matrix, gamma, alpha=1,
                                            penalty=penalty_coeffs[i])$epsilons
      CV_doubly_roubust_risks[i] <- (CV_doubly_roubust_risks[i] +
                                       var(evaluate_EIC(D[test_ids, ,], epsilons, Q_hat, V_hat, evaluation_action_matrix, gamma))
      )
    }
    if(greedy & (i > 2) &&
       (CV_doubly_roubust_risks[i] > CV_doubly_roubust_risks[i-1] 
        + 0.1 * abs(CV_doubly_roubust_risks[i-1] > CV_doubly_roubust_risks[i-2]) ) ) break # Enforce greedy search
  }
  if(plot_risk){
    print(CV_doubly_roubust_risks)
    plot(CV_doubly_roubust_risks)
  }
  # Output the estimate computed on the full dataset that uses the softening coeff that minimizes the CV risk
  list(estimate = penalized_LTMLE_estimator(D, Q_hat, V_hat, evaluation_action_matrix, gamma, alpha=1,
                                  penalty=penalty_coeffs[which.min(CV_doubly_roubust_risks)])$estimate,
       penalty = penalty_coeffs[ which.min(CV_doubly_roubust_risks) ])
}


# Debugging experiments ---------------------------------------------------
# Set DGP parameters
# horizon <- 20; n <- 1e2; gamma <- 1
# V0_and_Q0 <- compute_true_V_and_Q(state_transition_matrix,
#                                   transition_based_rewards,
#                                   evaluation_action_matrix, horizon, gamma)
# V0 <- V0_and_Q0$V0; Q0 <- V0_and_Q0$Q0
# b <- 0.1e-1*rnorm(1)
# Q_hat <- Q0 +  b
# V_hat <-  V0 +  b
# D <- generate_discrete_MDP_dataset(n, 1, state_transition_matrix,
#                                    behavior_action_matrix,
#                                    transition_based_rewards,
#                                    horizon)
# # debug(C_LTMLE_softening)
# # Check risk
# # D_large <- generate_discrete_MDP_dataset(1e4, 1, state_transition_matrix,
# #                                          behavior_action_matrix,
# #                                          transition_based_rewards,
# #                                          horizon)
# par(mfrow=c(2,1))
# C_LTMLE_softening(D, Q_hat, V_hat, evaluation_action_matrix, gamma, V=3, plot_risk=T, D_large=NULL)
# C_LTMLE_penalization(D, Q_hat, V_hat, evaluation_action_matrix, gamma, V=3, plot_risk=T, greedy=T)
# 



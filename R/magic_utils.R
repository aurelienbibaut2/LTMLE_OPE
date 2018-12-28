bootstrap_bias <- function(D, Q_hat, V_hat, horizon, gamma, number_bootstrap = 1e3, alpha = 0.1) {
  n <- dim(D)[1]
  wdr_bootstrap <- c()
  for (i in 1:number_bootstrap) {
    index_bootstrap <- sample(x = 1:n, size = n, replace = TRUE)
    D_bootstrap <- D[index_bootstrap, , ]
    wdr_bootstrap <- c(wdr_bootstrap, WDR_estimator_TB(D=D_bootstrap, Q_hat=Q_hat, V_hat=V_hat, 
                                                       gamma=gamma, j=horizon, compute_covariance=F)$g_js[horizon+1])
  }
  quantile(wdr_bootstrap, probs = c(alpha / 2, 1 - alpha / 2))
}


bootstrap_bias <- function(D, Q0, V0, number_bootstrap = 1e3) {
  n <- dim(D)[1]
  wdr_bootstrap <- c()
  for (i in 1:number_bootstrap) {
    index_bootstrap <- sample(x = 1:n, size = n, replace = TRUE)
    D_bootstrap <- D[index_bootstrap, , ]
    wdr_once <- WDR_estimator_TB(D = D_bootstrap, Q_hat = Q0, V_hat = V0)
    wdr_bootstrap <- c(wdr_bootstrap, wdr_once)
  }
  wdr_quantiles <- quantile(wdr_bootstrap, probs = c(0.05, 0.95))
  return(wdr_quantiles)  
}
# bootstrap_bias(D, Q0, V0)

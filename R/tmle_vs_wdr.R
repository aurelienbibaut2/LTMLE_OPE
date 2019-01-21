source('utils.R')
generate_data <- function(n, use_pi_e=F){
  W1 <- (2*rbinom(n, 1, 0.5) - 1) * rexp(n)
  W2 <- rnorm(n)
  pi_b <- expit(0.5 * W1)
  pi_e <- rep(0.5, n)
  if(use_pi_e){
    A <- rbinom(n, 1, prob = pi_e)
  }else{
    A <- rbinom(n, 1, prob = pi_b)
  }
  Y <- rbinom(n, 1, prob = expit((2*A - 1) + W2 + 2 * W1 * W2))
  
  cbind(W1=W1, W2=W2, A=A, Y=Y, pi_b=pi_b, pi_e=pi_e)
}

EIC <- function(D, Q_hat, Q_hat_a1, Q_hat_a2, epsilon){
  IC_term1 <- D[, 'pi_e'] / D[, 'pi_b' ] * (D[, 'Y'] - expit(logit(Q_hat) + epsilon) ) 
  IC_term2a <- 0.5 * expit(logit(Q_hat_a1) + epsilon) + 0.5 * expit(logit(Q_hat_a2) + epsilon)
  
  IC_term1 + IC_term2a - mean(IC_term2a)
}

targeting_risk <- function(D, Q_hat, epsilon, alpha, lambda){
  omega <- soften(D[, 'pi_e'] / D[, 'pi_b' ], alpha)
  (-mean(omega * ( D[, 'Y'] * log(expit(logit(Q_hat) + epsilon)) + 
                     ( 1- D[, 'Y']) * log(1 - expit(logit(Q_hat) + epsilon) ) )
         ) 
  + lambda  * abs(epsilon)
  )
}

TMLE <- function(D, Q_hat, Q_hat_a1, Q_hat_a2, alpha, lambda){
  omega <- soften(D[, 'pi_e'] / D[, 'pi_b' ], alpha)
  #epsilon <- glm(D[, 'Y'] ~ offset(logit(Q_hat)) + 1, family = quasibinomial(), weights = omega)$coefficients[1]
  
  epsilon <- optimize(function(x) targeting_risk(D, Q_hat, x, alpha, lambda), lower=-6, upper=6)$minimum
  
  mean(0.5 * expit(logit(Q_hat_a1) + epsilon) + 0.5 * expit(logit(Q_hat_a2) + epsilon))
}

# debug(TMLE)

IPS <- function(D, Q_hat){
  mean(D[, 'pi_e'] / D[, 'pi_b'] * D[, 'Y'])
}

WDR <- function(D, Q_hat, V_hat, alpha){
  omega <- n * soften(D[, 'pi_e'] / D[, 'pi_b' ], alpha)
  mean(omega / mean(omega) * (D[, 'Y'] - Q_hat)  + V_hat)
}

# Get truth
mean(generate_data(5e3, use_pi_e = T)[, 'Y'])

MSE_TMLE <- c()
MSE_WDR <- c()
n_reg_params <- 10
alphas <- seq(0, 1, length.out = n_reg_params)
lambdas <- rev(seq(0, 1e-4, length.out = n_reg_params))
for(k in 1:n_reg_params){
  TMLE_estimates <- c()
  WDR_estimates <- c()
  for(m in 1:100){
    n <- 1000
    D <- generate_data(n=n)
    Q_hat_model <- glm(Y ~ A + W1 , family = binomial(), data=data.frame(D) )
    Q_hat <- predict(Q_hat_model, newdata = data.frame(D), type='response')
    Da1 <- D; Da1[, 'A'] <- 0
    Da2 <- D; Da2[, 'A'] <- 1
    Q_hat_a1 <- predict(Q_hat_model, newdata = data.frame(Da1), type='response')
    Q_hat_a2<- predict(Q_hat_model, newdata = data.frame(Da2), type='response')
    V_hat <- (Q_hat_a1 + Q_hat_a2) /  2
    TMLE_estimates <- c(TMLE_estimates, TMLE(D, Q_hat, Q_hat_a1, Q_hat_a2, alphas[k], lambdas[k]))
    WDR_estimates <- c(WDR_estimates, WDR(D, Q_hat, V_hat, alphas[k]))
  }
  MSE_TMLE <- c(MSE_TMLE, mean((TMLE_estimates - 0.5)^2))
  MSE_WDR <- c(MSE_WDR, mean((WDR_estimates - 0.5)^2))
}

plot(alphas, log10(MSE_TMLE), col='green', type='b', lty=2, ylim=range(log10(c(MSE_TMLE, MSE_WDR))))
lines(alphas, log10(MSE_WDR), col='red', type='b', lty=2)
legend(0.01, 0.04, legend = c('TMLE', 'WDR'), col=c('green', 'red'), lty=c(2,2))
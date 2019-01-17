source('utils.R')
generate_data <- function(n, use_pi_e=F){
  W1 <- (2*rbinom(n, 1, 0.5) - 1) * rexp(n)
  W2 <- rnorm(n)
  pi_b <- expit(4 * W1)
  pi_e <- rep(0.5, n)
  if(use_pi_e){
    A <- rbinom(n, 1, prob = pi_e)
  }else{
    A <- rbinom(n, 1, prob = pi_b)
  }
  Y <- rbinom(n, 1, prob = expit((2*A - 1) + W2))
  
  cbind(W1=W1, W2=W2, A=A, Y=Y, pi_b=pi_b, pi_e=pi_e)
}

TMLE <- function(D, Q_hat){
  omega <- D[, 'pi_e'] / D[, 'pi_b' ]
  epsilon <- glm(D[, 'Y'] ~ offset(logit(Q_hat)) + 1, family = quasibinomial(), weights = omega)$coefficients[1]
  mean(expit(logit(Q_hat) + epsilon))
}

IPS <- function(D, Q_hat){
  mean(D[, 'pi_e'] / D[, 'pi_b']* D[, 'Y'])
}

WDR <- function(D, Q_hat){
  omega <- D[, 'pi_e'] / D[, 'pi_b' ]
  mean(omega / mean(omega) * (D[, 'Y'] - 0.3)  + 0.3)
}


# Get truth
mean(generate_data(5e3, use_pi_e = T)[, 'Y'])

TMLE_estimates <- c()
WDR_estimates <- c()
for(m in 1:1e3){
  n <- 100
  D <- generate_data(n=n)
  Q_hat_model <- glm(Y ~ A + W2, family = binomial(), data=data.frame(D) )
  Q_hat <- predict(Q_hat_model, newdata = data.frame(D), type='response')
  TMLE_estimates <- c(TMLE_estimates, TMLE(D, Q_hat))
  WDR_estimates <- c(WDR_estimates, WDR(D, Q_hat))
}

print(mean((TMLE_estimates - 0.5)^2))
print(mean((WDR_estimates - 0.5)^2))
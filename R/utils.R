########################################
## Run multiple estimators on the same
## set of trajectories H
##
## Possibly add noise to p_b,p_e,Q,V
##

RL_sim = function(H, state_transition_matrix, behavior_action_matrix, transition_based_rewards,
                  horizon, Q0=Q0, V0=V0, noiseb, noisee, noiseQ=noiseQ, noiseV=noiseV){
  
  eIS <- est_IS(H_all=H,
                behavior_action_matrix=behavior_action_matrix,
                evaluation_action_matrix=evaluation_action_matrix,
                horizon=horizon,
                noiseb=noiseb,
                noisee=noisee)
  
  eWIS <- est_WIS(H_all=H,
                  behavior_action_matrix=behavior_action_matrix,
                  evaluation_action_matrix=evaluation_action_matrix,
                  horizon=horizon,
                  noiseb=noiseb,
                  noisee=noisee)
  
  eDR <- est_DR(H_all=H,
                behavior_action_matrix=behavior_action_matrix,
                evaluation_action_matrix=evaluation_action_matrix,
                horizon=horizon,
                Q0=Q0,
                V0=V0,
                noiseQ=noiseQ,
                noiseV=noiseV,
                noiseb=noiseb,
                noisee=noisee)
  
  eWDR <- est_WDR(H_all=H,
                  behavior_action_matrix=behavior_action_matrix,
                  evaluation_action_matrix=evaluation_action_matrix,
                  horizon=horizon,
                  Q0=Q0,
                  V0=V0,
                  noiseQ=noiseQ,
                  noiseV=noiseV,
                  noiseb=noiseb,
                  noisee=noisee)
  
  return(list(eIS=eIS, eWIS=eWIS, eDR=eDR, eWDR=eWDR))
}

########################################
## Generate many trajectores 

gen_traj = function(M, 
                    state_transition_matrix,
                    behavior_action_matrix,
                    transition_based_rewards,
                    horizon){
  
  H <- list()
  for(m in 1:M){
    H[[m]]<- generate_discrete_MDP_trajectory(1, state_transition_matrix,
                                              behavior_action_matrix,
                                              transition_based_rewards,
                                              horizon)}
  
  return(H)
}

########################################
## Step-wise Importance Sampling. 
## Following paper 2015 by Nan Jiang and Nihong Li

est_IS <- function(H_all,behavior_action_matrix,evaluation_action_matrix,horizon,
                   noiseb=0,noisee=0){
  
  M<-length(H_all)
  V_IS_summands <- c()
  
  #Possibly add some noise to pi_b and pi_e
  behavior_action_matrix <- behavior_action_matrix + 
    rnorm(n=length(behavior_action_matrix), mean = 0, sd = noiseb)
  evaluation_action_matrix <- evaluation_action_matrix + 
    rnorm(n=length(evaluation_action_matrix), mean = 0, sd = noisee)
  
  for(m in 1:M){
    
    H <- H_all[[m]]
    
    pi_b <- apply(H, 1, function(x) behavior_action_matrix[x[1], x[2]])
    pi_a <- apply(H, 1, function(x) evaluation_action_matrix[x[1], x[2]])
    
    V_IS <- rep(NA, horizon)
    t <- horizon
    V_IS[t] <- pi_a[t] / pi_b[t] * (H[t, 'r'])
    
    for(t in (horizon-1):1){
      V_IS[t] <- V_IS[t+1] + pi_a[t] / pi_b[t] * (H[t, 'r'])
    }
    V_IS_summands <- c(V_IS_summands, V_IS[1])
  }
  
  MC_V_IS <- mean(V_IS_summands)
  return(MC_V_IS)
}

########################################
## Step-wise Weighted Importance Sampling. 
## Following paper 2015 by Nan Jiang and Nihong Li

est_WIS <- function(H_all,behavior_action_matrix,evaluation_action_matrix,horizon,
                    noiseb=0,noisee=0){
  
  M<-length(H_all)
  V_WIS_summands <- c()
  
  #Possibly add some noise to pi_b and pi_e
  behavior_action_matrix <- behavior_action_matrix + 
    rnorm(n=length(behavior_action_matrix), mean = 0, sd = noiseb)
  evaluation_action_matrix <- evaluation_action_matrix + 
    rnorm(n=length(evaluation_action_matrix), mean = 0, sd = noisee)
  
  pi_b <- do.call(rbind, lapply(H_all, function(h) 
    {apply(h, 1, function(x) behavior_action_matrix[x[1], x[2]])}))
  pi_a <- do.call(rbind, lapply(H_all, function(h) 
    {apply(h, 1, function(x) evaluation_action_matrix[x[1], x[2]])}))
  
  for(m in 1:M){
    
    t <- horizon
    w <- rep(NA, horizon)
    w[t]<-sum((pi_a[1:M,t] / pi_b[1:M,t]) / M)
    
    V_WIS <- rep(NA, horizon)
    
    V_WIS[t] <- (pi_a[m,t] / pi_b[m,t])/w[t] * (H_all[[m]][t, 'r'])
    
    for(t in (horizon-1):1){
      w[t]<-sum((pi_a[1:M,t] / pi_b[1:M,t]) / M)
      V_WIS[t] <- V_WIS[t+1] + (pi_a[m,t] / pi_b[m,t])/w[t] * (H_all[[m]][t, 'r'])
    }
    V_WIS_summands <- c(V_WIS_summands, V_WIS[1])
  }
  
  MC_V_WIS <- mean(V_WIS_summands)
  return(MC_V_WIS)
  
}

########################################
## Double Robust estimator. 
## Following 2015 paper by Nan Jiang and Nihong Li

est_DR <- function(H_all,behavior_action_matrix,evaluation_action_matrix,horizon,
                   Q0,V0,noiseQ=0,noiseV=0,noiseb=0,noisee=0){
  
  M<-length(H_all)
  V_DR_summands <- c()
  
  #Possibly add some noise to pi_b and pi_e
  behavior_action_matrix <- behavior_action_matrix + 
    rnorm(n=length(behavior_action_matrix), mean = 0, sd = noiseb)
  evaluation_action_matrix <- evaluation_action_matrix + 
    rnorm(n=length(evaluation_action_matrix), mean = 0, sd = noisee)
  
  #Possibly add some noise to Q and V ("estimated")
  Q0 <- Q0 + rnorm(n=length(Q0), mean = 0, sd = noiseQ)
  V0 <- V0 + rnorm(n=length(V0), mean = 0, sd = noiseV)
  
  for(m in 1:M){
    
    H <- H_all[[m]]
    t <- horizon
    
    pi_b <- apply(H, 1, function(x) behavior_action_matrix[x[1], x[2]])
    pi_a <- apply(H, 1, function(x) evaluation_action_matrix[x[1], x[2]])
  
    V_DR <- rep(NA, horizon)
    
    V_DR[t] <- V0[t, H[t, 's']] +  pi_a[t] / pi_b[t] * (H[t, 'r'] - Q0[t, H[t, 's'], H[t, 'a']])
    
    for(t in (horizon-1):1){
      V_DR[t] <- V0[t, H[t, 's']] +  pi_a[t] / pi_b[t] * (H[t, 'r'] + V_DR[t+1] - Q0[t, H[t, 's'], H[t, 'a']])
    }
    
    V_DR_summands <- c(V_DR_summands, V_DR[1])
  }
  
  MC_V_DR <- mean(V_DR_summands)
  return(MC_V_DR)
}

########################################
## Weighted Double Robust estimator. 
## Following 2016 paper by Philip Thomas and Emma Brunskill

est_WDR <- function(H_all,behavior_action_matrix,evaluation_action_matrix,horizon,
                    Q0,V0,noiseQ=0,noiseV=0,noiseb=0,noisee=0){
  
  M<-length(H_all)
  V_WDR_summands <- c()
  
  #Possibly add some noise to pi_b and pi_e
  behavior_action_matrix <- behavior_action_matrix + 
    rnorm(n=length(behavior_action_matrix), mean = 0, sd = noiseb)
  evaluation_action_matrix <- evaluation_action_matrix + 
    rnorm(n=length(evaluation_action_matrix), mean = 0, sd = noisee)
  
  pi_b <- do.call(rbind, lapply(H_all, function(h) 
    {apply(h, 1, function(x) behavior_action_matrix[x[1], x[2]])}))
  pi_a <- do.call(rbind, lapply(H_all, function(h) 
    {apply(h, 1, function(x) evaluation_action_matrix[x[1], x[2]])}))

  for(m in 1:M){
    
    t <- horizon
    w <- c(1,rep(NA, horizon))
    w[t+1]<-(pi_a[m,t]/pi_b[m,t]) / sum((pi_a[1:M,t] / pi_b[1:M,t])^2)
    
    for(i in (horizon-1):1){
      w[i+1]<-(pi_a[m,i]/pi_b[m,i]) / sum((pi_a[1:M,i] / pi_b[1:M,i])^2)
    }
    
    V_WDR <- rep(NA, horizon)
    
    #Possibly add some noise to Q and V ("estimated")
    Q0 <- Q0 + rnorm(n=length(Q0), mean = 0, sd = noiseQ)
    V0 <- V0 + rnorm(n=length(V0), mean = 0, sd = noiseV)
    
    V_WDR[t] <- w[t+1]*H_all[[m]][t, 'r'] - 
      (w[t+1]*Q0[t, H_all[[m]][t, 's'], H_all[[m]][t, 'a']] - w[t]* V0[t, H_all[[m]][t, 's']])
    
    for(t in (horizon-1):1){
      V_WDR[t] <- V_WDR[t+1] + w[t+1]*H_all[[m]][t, 'r'] - 
        (w[t+1]*Q0[t, H_all[[m]][t, 's'], H_all[[m]][t, 'a']] - w[t]* V0[t, H_all[[m]][t, 's']])
    }
    
    V_WDR_summands <- c(V_WDR_summands, V_WDR[1])
    
  }
  
  MC_V_WDR <- mean(V_WDR_summands)
  return(MC_V_WDR)
  
}

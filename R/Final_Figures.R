#Create prettier plots:
library(grid)
library(gridExtra)
library(here)
setwd(here("results"))

MSE_table_be_lb<-read.csv("Final/ModelWin_base_estimators_b0=5e-3_nr=71.csv")
MSE_table_be_hb<-read.csv("Final/ModelWin_base_estimators_b0=5e-2_nr=71.csv")
MSE_table_be_lb$estimator<-as.character(MSE_table_be_lb$estimator)
MSE_table_be_hb$estimator<-as.character(MSE_table_be_hb$estimator)
MSE_table_be_lb$estimator[MSE_table_be_lb$estimator == 'softened LTMLE'] <- 's LTMLE'
MSE_table_be_hb$estimator[MSE_table_be_hb$estimator == 'softened LTMLE'] <- 's LTMLE'

library(ggplot2)

p1 <- ggplot(data=MSE_table_be_lb, aes(x=id, y=log10(n*MSE), color=estimator, shape=estimator)) + geom_line() + 
  geom_point(size=2.5) + theme_bw() +
  ggtitle(paste("(a) ModelWin low bias")) +
  theme(legend.position="none") + 
  #ggtitle(paste(env_name,' horizon=', horizon, ', number of draws per point=', nb_repeats,
  #              '\nbias=', b0, '*rnorm(1), n=', ns[1])) +
  scale_shape_manual( values=c('WDR'=19, 'softened_WDR'=19, 'ps WDR'=19, 
                               'partial LTMLE'=15, 's LTMLE'=15, 'ps LTMLE'=15, 'psp LTMLE'=15) )

p2 <- ggplot(data=MSE_table_be_hb, aes(x=id, y=log10(n*MSE), color=estimator, shape=estimator)) + geom_line() + 
  geom_point(size=2.5) + theme_bw() +
  ggtitle(paste("(a) ModelWin high bias")) +
  #ggtitle(paste(env_name,' horizon=', horizon, ', number of draws per point=', nb_repeats,
  #              '\nbias=', b0, '*rnorm(1), n=', ns[1])) +
  scale_shape_manual( values=c('WDR'=19, 'softened_WDR'=19, 'ps WDR'=19, 
                               'partial LTMLE'=15, 's LTMLE'=15, 'ps LTMLE'=15, 'psp LTMLE'=15) )

grid.arrange(p1, p2, nrow = 1)

MSE_table_GW<-read.csv("GridWorld-1e2-1e3-horizon=100b0=0.005-71draws.csv")
MSE_table_MWlb<-read.csv("ModelWin-1e2-1e4-horizon=10-b0=5e-3-63draws.csv")
MSE_table_MWhb<-read.csv("ModelWin-1e2-1e4-horizon=10-b0=5e-2-63draws.csv")
MSE_table_MF<-read.csv("ModelFail-1e2-1e2.7-horizon=3b0=0.005-71draws.csv")
MSE_table_MWhb$estimator<-as.character(MSE_table_MWhb$estimator)
MSE_table_MWhb$estimator[MSE_table_MWhb$estimator == 'MAGIC_LTMLE'] <- 'RLTMLE 1'
MSE_table_MWhb$estimator[MSE_table_MWhb$estimator == 'MAGIC_bootstrap'] <- 'RLTMLE 2'

p1<-ggplot(data=MSE_table_GW, aes(x=log10(n), y=log10(n*MSE), color=estimator, shape=estimator)) + 
  scale_shape_manual( values=c('MAGIC'=15, 'MAGIC_full_library'=15, 'MAGIC_bootstrap'=15, 'RLTMLE 1'=15, 'RLTMLE 2'=15,
                               'C-TMLE-sftning'=15, '1step_LTMLE'=15, 'MAGIC_LTMLE'=15, 'partial_LTMLE_1.0'=15, 'partial_LTMLE_0.3'=15,
                               'LTMLE_1.0'=19, 'LTMLE_0.7'=19, 'LTMLE_0.5'=19, 'LTMLE_0.1'=19, 'LTMLE_0.0'=19,
                               'WDR'=18) ) +
  scale_size_manual( values=c('MAGIC'=8, 'MAGIC_full_library'=8, 'MAGIC_bootstrap'=8, 'RLTMLE 1'=8, 'RLTMLE 2'=8,
                              'C-TMLE-sftning'=8, '1step_LTMLE'=8, 'MAGIC_LTMLE'=8, 'partial_LTMLE_1.0'=8, 'partial_LTMLE_0.3'=8,
                              'LTMLE_1.0'=4, 'LTMLE_0.7'=4, 'LTMLE_0.5'=4, 'LTMLE_0.1'=4, 'LTMLE_0.0'=4, 
                              'WDR'=4)) +
  theme_bw() +
  geom_line(size=1) + 
  theme(legend.position="none") +
  geom_point(aes(size=estimator)) +
  ggtitle(paste("(a) GridWorld"))

p2<-ggplot(data=MSE_table_MF, aes(x=log10(n), y=log10(n*MSE), color=estimator, shape=estimator)) + 
  scale_shape_manual( values=c('MAGIC'=15, 'MAGIC_full_library'=15, 'MAGIC_bootstrap'=15, 'RLTMLE 1'=15, 'RLTMLE 2'=15,
                               'C-TMLE-sftning'=15, '1step_LTMLE'=15, 'MAGIC_LTMLE'=15, 'partial_LTMLE_1.0'=15, 'partial_LTMLE_0.3'=15,
                               'LTMLE_1.0'=19, 'LTMLE_0.7'=19, 'LTMLE_0.5'=19, 'LTMLE_0.1'=19, 'LTMLE_0.0'=19,
                               'WDR'=18) ) +
  scale_size_manual( values=c('MAGIC'=8, 'MAGIC_full_library'=8, 'MAGIC_bootstrap'=8, 'RLTMLE 1'=8, 'RLTMLE 2'=8,
                              'C-TMLE-sftning'=8, '1step_LTMLE'=8, 'MAGIC_LTMLE'=8, 'partial_LTMLE_1.0'=8, 'partial_LTMLE_0.3'=8,
                              'LTMLE_1.0'=4, 'LTMLE_0.7'=4, 'LTMLE_0.5'=4, 'LTMLE_0.1'=4, 'LTMLE_0.0'=4, 
                              'WDR'=4)) +
  theme_bw() +
  theme(legend.position="none") +
  geom_line(size=1) + geom_point(aes(size=estimator)) +
  ggtitle(paste("(b) ModelFail"))

p3<-ggplot(data=MSE_table_MWlb, aes(x=log10(n), y=log10(n*MSE), color=estimator, shape=estimator)) + 
  scale_shape_manual( values=c('MAGIC'=15, 'MAGIC_full_library'=15, 'MAGIC_bootstrap'=15, 'RLTMLE 1'=15, 'RLTMLE 2'=15,
                               'C-TMLE-sftning'=15, '1step_LTMLE'=15, 'MAGIC_LTMLE'=15, 'partial_LTMLE_1.0'=15, 'partial_LTMLE_0.3'=15,
                               'LTMLE_1.0'=19, 'LTMLE_0.7'=19, 'LTMLE_0.5'=19, 'LTMLE_0.1'=19, 'LTMLE_0.0'=19,
                               'WDR'=18) ) +
  scale_size_manual( values=c('MAGIC'=8, 'MAGIC_full_library'=8, 'MAGIC_bootstrap'=8, 'RLTMLE 1'=8, 'RLTMLE 2'=8,
                              'C-TMLE-sftning'=8, '1step_LTMLE'=8, 'MAGIC_LTMLE'=8, 'partial_LTMLE_1.0'=8, 'partial_LTMLE_0.3'=8,
                              'LTMLE_1.0'=4, 'LTMLE_0.7'=4, 'LTMLE_0.5'=4, 'LTMLE_0.1'=4, 'LTMLE_0.0'=4, 
                              'WDR'=4)) +
  theme_bw() +
  geom_line(size=1) + geom_point(aes(size=estimator)) +
  ggtitle(paste("(c) ModelWin low bias"))

p4<-ggplot(data=MSE_table_MWhb, aes(x=log10(n), y=log10(n*MSE), color=estimator, shape=estimator)) + 
  scale_shape_manual( values=c('MAGIC'=15, 'MAGIC_full_library'=15, 'MAGIC_bootstrap'=15, 'RLTMLE 1'=15, 'RLTMLE 2'=15,
                               'C-TMLE-sftning'=15, '1step_LTMLE'=15, 'MAGIC_LTMLE'=15, 'partial_LTMLE_1.0'=15, 'partial_LTMLE_0.3'=15,
                               'LTMLE_1.0'=19, 'LTMLE_0.7'=19, 'LTMLE_0.5'=19, 'LTMLE_0.1'=19, 'LTMLE_0.0'=19,
                               'WDR'=18) ) +
  scale_size_manual( values=c('MAGIC'=8, 'MAGIC_full_library'=8, 'MAGIC_bootstrap'=8, 'RLTMLE 1'=8, 'RLTMLE 2'=8,
                              'C-TMLE-sftning'=8, '1step_LTMLE'=8, 'MAGIC_LTMLE'=8, 'partial_LTMLE_1.0'=8, 'partial_LTMLE_0.3'=8,
                              'LTMLE_1.0'=4, 'LTMLE_0.7'=4, 'LTMLE_0.5'=4, 'LTMLE_0.1'=4, 'LTMLE_0.0'=4, 
                              'WDR'=4)) +
  theme_bw() +
  geom_line(size=1) + geom_point(aes(size=estimator)) +
  ggtitle(paste("(d) ModelWin high bias"))

grid.arrange(p1, p3, p2, p4, nrow = 2)

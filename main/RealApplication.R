pkgs <- c('openxlsx','copula','R2OpenBUGS','R2jags','coda','MCMCvis','readxl','survivalROC',
          'DescTools','dplyr','ggplot2','gridExtra','grid','devtools','latex2exp', 'risksetROC')
lapply(pkgs, require, character.only = TRUE)
load_all(".")
f.path='C:/Users/farea/Documents/RWorkspace'

# Storing data
library(survivalROC)

data(mayo)
mayo$time <- mayo$time/365 #convert to year
mayo <- mayo[,-4]
colnames(mayo) <- c('t','status','x')

set.seed(123456)

# Finding the best model to use for PH model
p1.jags <- c('lambda.T','beta','mu.X','Deviance')
i1.jags <- function(){list(lambda.T=runif(1),beta=.1, mu.X=runif(1))}
MCMCestimate_PH('exp_exp', mayo, f.path=f.path, p.jags=p1.jags, i.jags=i1.jags)
par.est <- read.table(file.path(f.path,"Simulation_result","MCMC_ph_expo_expo_result.txt"))
params.est <- list(mu.x=par.est[5,1],lambda.t=par.est[4,1],beta=par.est[1,1], dev=par.est[3,1])

# calculating log-likehood & Dhat
ph.const <- exp(params.est$beta*mayo$x)
logHaz <- log(params.est$lambda.t * ph.const)
logSurv <- log(exp(-params.est$lambda.t*mayo$t)^ph.const)
logfX <- log(params.est$mu.x*exp(-params.est$mu.x*mayo$x))
L <- mayo$status * (logHaz + logfX) + logSurv
Dhat <- sum(-2*L)
DIC.1 <- 2*params.est$dev - Dhat

p2.jags <- c('lambda.T','beta','mu.X','sigma.X','Deviance')
i2.jags <- function(){list(lambda.T=runif(1),beta=.1, mu.X=runif(1), tau=runif(1))}
MCMCestimate_PH('norm_exp', mayo, f.path=f.path, p.jags=p2.jags, i.jags=i2.jags)
par.est <- read.table(file.path(f.path,"Simulation_result","MCMC_ph_norm_expo_result.txt"))
params.est <- list(mu.x=par.est[5,1],lambda.t=par.est[4,1],beta=par.est[1,1],sigma.x=par.est[6,1], dev=par.est[3,1])

# calculating log-likehood & Dhat
ph.const <- exp(params.est$beta*mayo$x)
logHaz <- log(params.est$lambda.t * ph.const)
logSurv <- -params.est$lambda.t *mayo$t*ph.const
logfX <- -(mayo$x-params.est$mu.x)^2/(2*(params.est$sigma.x^2)) - log(params.est$sigma.x*(2*(22/7))^(1/2))
L <- mayo$status * (logHaz + logfX) + logSurv
Dhat <- sum(-2*L)
DIC.2 <- 2*params.est$dev - Dhat


p3.jags <- c('shape.T','scale.T','beta','mu.X','sigma.X','Deviance')
i3.jags <- function(){list(shape.T=runif(1),scale.T=runif(1),beta=.1, mu.X=runif(1), tau=runif(1))}
MCMCestimate_PH('norm_weib', mayo, f.path=f.path, p.jags=p3.jags, i.jags=i3.jags)
par.est <- read.table(file.path(f.path,"Simulation_result","MCMC_ph_norm_weib_result.txt"))
params.est <- list(mu.x=par.est[4,1],scale.t=par.est[5,1],beta=par.est[1,1],shape.t=par.est[6,1],sigma.x=par.est[7,1], dev=par.est[3,1])

# calculating log-likehood & Dhat
ph.const <- exp(params.est$beta*mayo$x)
logHaz <- log(params.est$shape.t) + log(params.est$scale.t) + (params.est$shape.t-1)*log(params.est$scale.t*mayo$t) + log(ph.const)
logSurv <- -(params.est$scale.t*mayo$t)^params.est$shape.t*ph.const
logfX <- -(mayo$x-params.est$mu.x)^2/(2*(params.est$sigma.x^2)) - log(params.est$sigma.x*(2*(22/7))^(1/2))
L <- mayo$status * (logHaz + logfX) + logSurv
Dhat <- sum(-2*L)
DIC.3 <- 2*params.est$dev - Dhat

# Finding the best model to use for Copula function
p1.jags <- c('lambda.T','theta','mu.X','Deviance')
i1.jags <- function(){list(lambda.T=runif(1),theta=1, mu.X=runif(1))}
MCMCestimate_COP('exp_exp', mayo, f.path=f.path, p.jags=p1.jags, i.jags=i1.jags)
par.est <- read.table(file.path(f.path,"Simulation_result","MCMC_cop_expo_expo_result.txt"))
params.est <- list(mu.x=par.est[4,1],lambda.t=par.est[3,1],theta=par.est[5,1], dev=par.est[2,1])

# calculating log-likehood & Dhat
u <- pexp(mayo$x, rate=params.est$mu.x)
v <- pexp(mayo$t, rate=params.est$lambda.t)
log.c <- log(1+params.est$theta) + (-1/params.est$theta-2)*log((1-u)^-params.est$theta + v^-params.est$theta - 1) + (-params.est$theta-1)*(log(1-u) + log(v))
log.C1 <- (-1/params.est$theta-1)*log((1-u)^-params.est$theta + v^-params.est$theta - 1) + (-params.est$theta-1)*log(1-u)
log.fx <- dexp(mayo$x, rate=params.est$mu.x, log=T)
log.ft <- dexp(mayo$t, rate=params.est$lambda.t, log=T)
L <- mayo$status * (log.c+log.fx+log.ft) + (1-mayo$status)*(log.C1+log.fx)
Dhat <- sum(-2*L)
DIC.4 <- 2*params.est$dev - Dhat

p2.jags <- c('lambda.T','theta','mu.X','sigma.X','Deviance')
i2.jags <- function(){list(lambda.T=runif(1),theta=1, mu.X=6, tau=runif(1))}
MCMCestimate_COP('norm_exp', mayo, f.path=f.path, p.jags=p2.jags, i.jags=i2.jags)
par.est <- read.table(file.path(f.path,"Simulation_result","MCMC_cop_norm_expo_result.txt"))
params.est <- list(mu.x=par.est[4,1],lambda.t=par.est[3,1],theta=par.est[6,1],sigma.x=par.est[5,1], dev=par.est[2,1])

# calculating log-likehood & Dhat
u <- pnorm(mayo$x, mean=params.est$mu.x, sd=params.est$sigma.x)
v <- pexp(mayo$t, rate=params.est$lambda.t)
log.c <- log(1+params.est$theta) + (-1/params.est$theta-2)*log((1-u)^-params.est$theta + v^-params.est$theta - 1) + (-params.est$theta-1)*(log(1-u) + log(v))
log.C1 <- (-1/params.est$theta-1)*log((1-u)^-params.est$theta + v^-params.est$theta - 1) + (-params.est$theta-1)*log(1-u)
log.fx <- dnorm(mayo$x, mean=params.est$mu.x, sd=params.est$sigma.x, log=T)
log.ft <- dexp(mayo$t, rate=params.est$lambda.t, log=T)
L <- mayo$status * (log.c+log.fx+log.ft) + (1-mayo$status)*(log.C1+log.fx)
Dhat <- sum(-2*L)
DIC.5 <- 2*params.est$dev - Dhat

p3.jags <- c('shape.T','scale.T','theta','mu.X','sigma.X','Deviance')
i3.jags <- function(){list(shape.T=1.5,k=0.2,theta=1, mu.X=6, tau=runif(1))}
MCMCestimate_COP('norm_weib', mayo, f.path=f.path, p.jags=p3.jags, i.jags=i3.jags)
par.est <- read.table(file.path(f.path,"Simulation_result","MCMC_cop_norm_weib_result.txt"))
params.est <- list(mu.x=par.est[5,1],shape.t=par.est[3,1],theta=par.est[7,1],sigma.x=par.est[6,1],scale.t=par.est[4,1], dev=par.est[2,1])

# calculating log-likehood & Dhat
u <- pnorm(mayo$x, mean=params.est$mu.x, sd=params.est$sigma.x)
v <- pweibull(mayo$t, shape=params.est$shape.t, scale=1/params.est$scale.t)
log.c <- log(1+params.est$theta) + (-1/params.est$theta-2)*log((1-u)^-params.est$theta + v^-params.est$theta - 1) + (-params.est$theta-1)*(log(1-u) + log(v))
log.C1 <- (-1/params.est$theta-1)*log((1-u)^-params.est$theta + v^-params.est$theta - 1) + (-params.est$theta-1)*log(1-u)
log.fx <- dnorm(mayo$x, mean=params.est$mu.x, sd=params.est$sigma.x, log=T)
log.ft <- dweibull(mayo$t, shape=params.est$shape.t, scale=1/params.est$scale.t, log=T)
L <- mayo$status * (log.c+log.fx+log.ft) + (1-mayo$status)*(log.C1+log.fx)
Dhat <- sum(-2*L)
DIC.6 <- 2*params.est$dev - Dhat

# ===================================================================================
# Both MCMC process for PH model and Copula function suggests to use the joint
# model from Normal-Weibull scenario since its DIC value is the lowest compare to
# the other 2 scenario
# ===================================================================================

#showing the ROC and AUC
n.quant=21; monte.sample=1000000
t <- quantile(mayo$t,probs=c(0.1,0.25,0.4,0.6,0.75,0.9))

set.seed(123456)
par.est <- read.table(file.path(f.path,"Simulation_result","MCMC_ph_norm_weib_result.txt"))
PH.est <- list(mu.x=par.est[4,1],scale.t=par.est[5,1],beta=par.est[1,1],shape.t=par.est[6,1],sigma.x=par.est[7,1], dev=par.est[3,1])
PH.roc <- roc_PH(max.x=max(mayo$x),time.t=t,mod='norm_weib',
                      params=PH.est,n.quant=n.quant,n.sample=monte.sample)
PH.auc <- PH.roc %>% group_by(time) %>% summarize(AUC=AUC(1-specificity,sensitivity)) %>% data.frame()
PH.auc$model <- 'PH'

par.est <- read.table(file.path(f.path,"Simulation_result","MCMC_cop_norm_weib_result.txt"))
COP.est <- list(mu.x=par.est[5,1],shape.t=par.est[3,1],theta=par.est[7,1],sigma.x=par.est[6,1],scale.t=par.est[4,1], dev=par.est[2,1])
COP.roc <- roc_COP(max.x=max(mayo$x), time.t=t, mod='norm_weib', params=COP.est, n.quant=n.quant)
COP.auc <- COP.roc %>% group_by(time) %>% summarize(AUC=AUC(1-specificity,sensitivity)) %>% data.frame()
COP.auc$model <- 'COP'

KM.auc <- data.frame()
for (i in t){
  res <- survivalROC(Stime=mayo$t, status=mayo$status, marker=mayo$x, predict.time=i, method='KM')$AUC
  KM.auc <- rbind(KM.auc, c(i, res))
}
colnames(KM.auc) <- c('time','AUC')
KM.auc$model <- 'KM'

RiskROC.auc <- data.frame()
for (i in t){
  roc <- risksetROC(Stime=mayo$t, status=mayo$status, marker=mayo$x, predict.time=i)$AUC
  RiskROC.auc <- rbind(RiskROC.auc, c(i,roc))
}
colnames(RiskROC.auc) <- c('time','AUC')
RiskROC.auc$model <- 'Risk Set'

# ROC curve
p1 <- ggplot(PH.roc, aes(x=1-specificity, y=sensitivity, group=time, color=time)) +
  geom_line(linewidth=0.8) +
  geom_point(shape=21) +
  scale_colour_gradient(high = "#132B43", low = "#56B1F7") +
  xlim(0,1) +
  ylim(0,1) +
  geom_abline(slope=1, linetype='dashed') +
  labs(title='Proportional Hazard Model') +
  theme(legend.position='none')

p2 <- ggplot(COP.roc, aes(x=1-specificity, y=sensitivity, group=time, color=time)) +
  geom_line(linewidth=0.8) +
  geom_point(shape=21) +
  scale_colour_gradient(high = "#132B43", low = "#56B1F7") +
  xlim(0,1) +
  ylim(0,1) +
  geom_abline(slope=1, linetype='dashed') +
  labs(y=element_blank(), title='Copula Function', color='Time (years)') +
  theme(axis.text.y=element_blank(),legend.position=c(.7,0.2), legend.direction='horizontal')

grid.arrange(p1,p2, nrow=1, ncol=2,
             top = textGrob("Time-dependent ROC curve", gp=gpar(fontsize=20, font=3)))

# AUC trend
AUC.res <- rbind(PH.auc,COP.auc,KM.auc, RiskROC.auc)
ggplot(AUC.res, aes(x=time, y=AUC)) +
  geom_line(linewidth = 1, aes(linetype=model)) +
  scale_linetype_manual(values=c('dotted','solid','dashed', 'dotdash'))+
  ylim(0.5,1)+
  labs(x='Time (years)',title='Comparison for AUC trend') +
  theme(legend.position=c(.5,0.1), legend.direction='horizontal')

pkgs <- c('openxlsx','copula','R2OpenBUGS','R2jags','coda','MCMCvis','readxl',
          'DescTools','dplyr','ggplot2','gridExtra','grid','devtools','latex2exp')
lapply(pkgs, require, character.only = TRUE)
load_all(".")

# =============================================
# Scripts below is used to study the impact of sample size towards the time-dependent ROC curve produced.
# We will see whether degeneracy problem in the ROC curve exists or not with our model that is
# derived from PH and Copula method.
# =============================================

n.quant=21; monte.sample=100000; data.sample=c(30,60,120,200); c.rate=0.2; beta <- 1.2; tau <- 0.7; theta <- 2*tau/(1-tau)
f.path='C:/Users/farea/Documents/RWorkspace'
# =========================================
# Simulation study for Expo-Expo scenario
# =========================================
ROC1_res <- NA
store.param1 <- data.frame()
mu.x=1.4; lambda.t=0.2;

set.seed(123456)
for(i in 1:length(data.sample)){
  params=list(mu.x=mu.x,lambda.t=lambda.t,c.rate=c.rate,beta=beta)
  df.PH <- generate_dataPH('exp_exp', n.sample=data.sample[i],
                           params=params)
  p.jags <- c('lambda.T','beta','mu.X','Deviance')
  i.jags <- function(){list(lambda.T=lambda.t,beta=beta, mu.X=mu.x)}
  MCMCestimate_PH(model='exp_exp', df=df.PH, f.path=f.path, p.jags=p.jags, i.jags=i.jags)
  par.est <- read.table(file.path(f.path,"Simulation_result","MCMC_ph_expo_expo_result.txt"))
  params.est <- list(mu.x=par.est[5,1],lambda.t=par.est[4,1],beta=par.est[1,1])
  store.param1 <- rbind(store.param1,params.est)
  t <- quantile(df.PH$t,probs=0.75)
  roc.exp_exp <- roc_PH(max.x=max(df.PH$x),time.t=t,mod='exp_exp',
                        params=params.est,n.quant=n.quant,n.sample=monte.sample)
  roc.exp_exp$n.sam <- rep(data.sample[i],n.quant)
  ROC1_res <- rbind(ROC1_res,roc.exp_exp)
}

ROC1_resCOP <- NA
store.paramCOP1 <- data.frame()
mu.x <- 1.4; lambda.t <- 0.2;
set.seed(123456)
for(i in 1:length(data.sample)){
  params=list(mu.x=mu.x,lambda.t=lambda.t,c.rate=c.rate,theta=theta)
  df.COP <- generate_dataCOP('exp_exp', n.sample=data.sample[i],params=params)
  p.jags <- c('lambda.T','theta','mu.X','Deviance')
  i.jags <- function(){list(lambda.T=lambda.t,theta=theta, mu.X=mu.x)}
  MCMCestimate_COP(model='exp_exp', df=df.COP, f.path=f.path, p.jags=p.jags, i.jags=i.jags)
  par.est <- read.table(file.path(f.path,"Simulation_result","MCMC_cop_expo_expo_result.txt"))
  params.est <- list(mu.x=par.est[4,1],lambda.t=par.est[3,1],theta=par.est[5,1])
  store.paramCOP1 <- rbind(store.paramCOP1,params.est)
  t <- quantile(df.COP$t,probs=0.75)
  roc.exp_exp <- roc_COP(max.x=max(df.COP$x), time.t=t, mod='exp_exp', params=params.est, n.quant=n.quant)
  roc.exp_exp$n.sam <- rep(data.sample[i],n.quant)
  ROC1_resCOP <- rbind(ROC1_resCOP,roc.exp_exp)
}


# =========================================
# Simulation study for Norm-Expo scenario
# =========================================
ROC2_res <- NA
store.param2 <- data.frame()
mu.x=5; sigma.x=0.8; lambda.t=0.2;

set.seed(123456)
for(i in 1:length(data.sample)){
  params=list(mu.x=mu.x,lambda.t=lambda.t,c.rate=c.rate,beta=beta, sigma.x=sigma.x)
  df.PH <- generate_dataPH('norm_exp', n.sample=data.sample[i], params=params)
  p.jags <- c('lambda.T','beta','mu.X','sigma.X','Deviance')
  i.jags <- function(){list(lambda.T=lambda.t,beta=beta, mu.X=mu.x, tau=runif(1))}
  MCMCestimate_PH(model='norm_exp', df=df.PH, f.path=f.path, p.jags=p.jags, i.jags=i.jags)
  par.est <- read.table(file.path(f.path,"Simulation_result","MCMC_ph_norm_expo_result.txt"))
  params.est <- list(mu.x=par.est[5,1],lambda.t=par.est[4,1],beta=par.est[1,1],sigma.x=par.est[6,1])
  store.param2 <- rbind(store.param2,params.est)
  t <- quantile(df.PH$t,probs=0.75)
  roc.norm_exp <- roc_PH(max.x=max(df.PH$x), time.t=t, mod='norm_exp',
                         params=params.est, n.quant=n.quant, n.sample=monte.sample)
  roc.norm_exp$n.sam <- rep(data.sample[i],n.quant)
  ROC2_res <- rbind(ROC2_res,roc.norm_exp)
}

ROC2_resCOP <- NA
store.paramCOP2 <- data.frame()
mu.x=5; sigma.x=0.8; lambda.t=0.2;
set.seed(123456)
for(i in 1:length(data.sample)){
  params=list(mu.x=mu.x,lambda.t=lambda.t,c.rate=c.rate,theta=theta, sigma.x=sigma.x)
  df.COP <- generate_dataCOP('norm_exp', n.sample=data.sample[i], params=params)
  p.jags <- c('lambda.T','theta','mu.X','sigma.X','Deviance')
  i.jags <- function(){list(lambda.T=lambda.t,theta=theta, mu.X=mu.x, tau=runif(1))}
  MCMCestimate_COP(model='norm_exp', df=df.COP, f.path=f.path, p.jags=p.jags, i.jags=i.jags)
  par.est <- read.table(file.path(f.path,"Simulation_result","MCMC_cop_norm_expo_result.txt"))
  params.est <- list(mu.x=par.est[4,1],lambda.t=par.est[3,1],theta=par.est[6,1],sigma.x=par.est[5,1])
  store.paramCOP2 <- rbind(store.paramCOP2,params.est)
  t <- quantile(df.COP$t,probs=0.75)
  roc.norm_exp <- roc_COP(max.x=max(df.COP$x), time.t=t, mod='norm_exp', params=params.est,n.quant=n.quant)
  roc.norm_exp$n.sam <- rep(data.sample[i],n.quant)
  ROC2_resCOP <- rbind(ROC2_resCOP,roc.norm_exp)
}


# =========================================
# Simulation study for Norm-Weib scenario
# =========================================
ROC3_res <- NA
store.param3 <- data.frame()
mu.x=5; sigma.x=0.8; shape.t=1.2; scale.t=1.4;

set.seed(123456)
for(i in 1:length(data.sample)){
  params=list(mu.x=mu.x, shape.t=shape.t, scale.t=scale.t, c.rate=c.rate, beta=beta, sigma.x=sigma.x)
  df.PH <- generate_dataPH('norm_weib', n.sample=data.sample[i], params=params)
  p.jags <- c('shape.T','scale.T','beta','mu.X','sigma.X','Deviance')
  i.jags <- function(){list(shape.T=shape.t,scale.T=scale.t,beta=beta, mu.X=mu.x, tau=runif(1))}
  MCMCestimate_PH(model='norm_weib', df=df.PH, f.path=f.path, p.jags=p.jags, i.jags=i.jags)
  par.est <- read.table(file.path(f.path,"Simulation_result","MCMC_ph_norm_weib_result.txt"))
  params.est <- list(mu.x=par.est[4,1],scale.t=par.est[5,1],beta=par.est[1,1],shape.t=par.est[6,1],sigma.x=par.est[7,1])
  store.param3 <- rbind(store.param3,params.est)
  t <- quantile(df.PH$t,probs=0.75)
  roc.norm_weib <- roc_PH(max.x=max(df.PH$x), time.t=t, mod='norm_weib',
                          params=params.est, n.quant=n.quant, n.sample=monte.sample)
  roc.norm_weib$n.sam <- rep(data.sample[i],n.quant)
  ROC3_res <- rbind(ROC3_res,roc.norm_weib)
}

ROC3_resCOP <- NA
store.paramCOP3 <- data.frame()
mu.x=5; sigma.x=0.8; shape.t=1.2; scale.t=1/1.4;
set.seed(123456)
for(i in 1:length(data.sample)){
  params=list(mu.x=mu.x, sigma.x=sigma.x, shape.t=shape.t,
              scale.t=scale.t, c.rate=c.rate, theta=theta)
  df.COP <- generate_dataCOP('norm_weib', n.sample=data.sample[i],params=params)
  p.jags <- c('scale.T','shape.T','theta','mu.X','sigma.X','Deviance')
  i.jags <- function(){list(k=scale.t,shape.T=shape.t,theta=theta, mu.X=mu.x, tau=runif(1))}
  MCMCestimate_COP(model='norm_weib', df=df.COP, f.path=f.path, p.jags=p.jags, i.jags=i.jags)
  par.est <- read.table(file.path(f.path,"Simulation_result","MCMC_cop_norm_weib_result.txt"))
  params.est <- list(mu.x=par.est[3,1],shape.t=par.est[5,1],theta=par.est[7,1],sigma.x=par.est[6,1],scale.t=par.est[4,1])
  store.paramCOP3 <- rbind(store.paramCOP3,params.est)
  t <- quantile(df.COP$t,probs=0.75)
  roc.norm_weib <- roc_COP(max.x=max(df.COP$x), time.t=t,mod='norm_weib', params=params.est,n.quant=n.quant)
  roc.norm_weib$n.sam <- rep(data.sample[i],n.quant)
  ROC3_resCOP <- rbind(ROC3_resCOP,roc.norm_weib)
}

# ================
# ROC PH model
# ================
ROC1_res <- ROC1_res[-1,]
ROC2_res <- ROC2_res[-1,]
ROC3_res <- ROC3_res[-1,]

p1 <- ggplot(ROC1_res, aes(x=1-specificity, y=sensitivity)) +
  geom_line(linewidth=0.8, aes(linetype=as.factor(n.sam))) +
  scale_linetype_manual(values=c("twodash","solid","dotdash","dotted")) +
  xlim(0,1) +
  ylim(0,1) +
  geom_abline(slope=1, linetype='dashed') +
  labs(x=element_blank(),y=element_blank(),title='Exponential-Exponential') +
  theme(legend.position='none')

p2 <- ggplot(ROC2_res, aes(x=1-specificity, y=sensitivity)) +
  geom_line(linewidth=0.8, aes(linetype=as.factor(n.sam))) +
  scale_linetype_manual(values=c("twodash","solid","dotdash","dotted")) +
  xlim(0,1) +
  ylim(0,1) +
  geom_abline(slope=1, linetype='dashed') +
  labs(x=element_blank(),y=element_blank(), title='Normal-Exponential') +
  theme(axis.text.y=element_blank(), legend.position='none')

p3 <- ggplot(ROC3_res, aes(x=1-specificity, y=sensitivity)) +
  geom_line(linewidth=0.8, aes(linetype=as.factor(n.sam))) +
  scale_linetype_manual(values=c("twodash","solid","dotdash","dotted")) +
  xlim(0,1) +
  ylim(0,1) +
  geom_abline(slope=1, linetype='dashed') +
  labs(x=element_blank(),y=element_blank(), title='Normal-Weibull', lty="Sample size:") +
  theme(axis.text.y=element_blank(),legend.position=c(.85,0.2), legend.direction='vertical')

grid.arrange(p1,p2,p3,nrow=1,ncol=3,
             top = textGrob("Time-dependent ROC with various sample size:\n PH model",gp=gpar(fontsize=20,font=3)),
             bottom=textGrob("1-Specificity",gp=gpar(fontsize=20,font=3)),
             left=textGrob("Sensitivity",gp=gpar(fontsize=20,font=3), rot=90))

# # ================
# # ROC Copula model
# # ================
ROC1_resCOP <- ROC1_resCOP[-1,]
ROC2_resCOP <- ROC2_resCOP[-1,]
ROC3_resCOP <- ROC3_resCOP[-1,]

p1 <- ggplot(ROC1_resCOP, aes(x=1-specificity, y=sensitivity)) +
  geom_line(linewidth=0.8, aes(linetype=as.factor(n.sam))) +
  scale_linetype_manual(values=c("twodash","solid","dotdash","dotted")) +
  xlim(0,1) +
  ylim(0,1) +
  geom_abline(slope=1, linetype='dashed') +
  labs(x=element_blank(),y=element_blank(),title='Exponential-Exponential') +
  theme(legend.position='none')

p2 <- ggplot(ROC2_resCOP, aes(x=1-specificity, y=sensitivity)) +
  geom_line(linewidth=0.8, aes(linetype=as.factor(n.sam))) +
  scale_linetype_manual(values=c("twodash","solid","dotdash","dotted")) +
  xlim(0,1) +
  ylim(0,1) +
  geom_abline(slope=1, linetype='dashed') +
  labs(x=element_blank(),y=element_blank(), title='Normal-Exponential') +
  theme(axis.text.y=element_blank(), legend.position='none')

p3 <- ggplot(ROC3_resCOP, aes(x=1-specificity, y=sensitivity)) +
  geom_line(linewidth=0.8, aes(linetype=as.factor(n.sam))) +
  scale_linetype_manual(values=c("twodash","solid","dotdash","dotted")) +
  xlim(0,1) +
  ylim(0,1) +
  geom_abline(slope=1, linetype='dashed') +
  labs(x=element_blank(),y=element_blank(), title='Normal-Weibull') +
  theme(axis.text.y=element_blank(),legend.position=c(.85,0.2), legend.direction='vertical')

grid.arrange(p1,p2,p3,nrow=1,ncol=3,
             top = textGrob("Time-dependent ROC with various sample size:\n Copula Function",gp=gpar(fontsize=20,font=3)),
             bottom=textGrob("1-Specificity",gp=gpar(fontsize=20,font=3)),
             left=textGrob("Sensitivity",gp=gpar(fontsize=20,font=3), rot=90))

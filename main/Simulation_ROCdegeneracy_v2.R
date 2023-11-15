pkgs <- c('openxlsx','copula','R2OpenBUGS','R2jags','coda','MCMCvis','readxl',
          'DescTools','dplyr','ggplot2','gridExtra','grid','devtools','latex2exp',
          'survivalROC','risksetROC')
lapply(pkgs, require, character.only = TRUE)
load_all(".")

# =============================================
# Scripts below is used to study the impact of sample size towards the time-dependent ROC curve produced.
# We will see whether degeneracy problem in the ROC curve exists or not with our model that is
# derived from PH and Copula method.
# =============================================

n.quant=21; monte.sample=100000; data.sample=200; c.rate=0.2; beta <- 1.2; tau <- 0.7; theta <- 2*tau/(1-tau)
# =========================================
# Simulation study for Expo-Expo scenario
# =========================================
mu.x=1.4; lambda.t=0.2;
set.seed(123456)
params=list(mu.x=mu.x,lambda.t=lambda.t,c.rate=c.rate,beta=beta)
df.PH <- generate_dataPH('exp_exp', n.sample=data.sample,
                           params=params)
t <- quantile(df.PH$t,probs=c(0.1,0.25,0.4,0.6,0.75,0.9))
ROC1_res <- roc_PH(max.x=max(df.PH$x),time.t=t,mod='exp_exp',
                      params=params,n.quant=n.quant,n.sample=monte.sample)
KM1.auc <- data.frame()
for (i in t){
  res <- survivalROC(Stime=df.PH$t, status=df.PH$status, marker=df.PH$x, predict.time=i, method='KM')$AUC
  KM1.auc <- rbind(KM1.auc, c(i, res))
}
RiskROC1.auc <- data.frame()
for (i in t){
  roc <- risksetROC(Stime=df.PH$t, status=df.PH$status, marker=df.PH$x, predict.time=i)$AUC
  RiskROC1.auc <- rbind(RiskROC1.auc, c(i,roc))
}

mu.x <- 1.4; lambda.t <- 0.2;
set.seed(123456)
params=list(mu.x=mu.x,lambda.t=lambda.t,c.rate=c.rate,theta=theta)
df.COP <- generate_dataCOP('exp_exp', n.sample=data.sample,params=params)
t <- quantile(df.COP$t,probs=c(0.1,0.25,0.4,0.6,0.75,0.9))
ROC1_resCOP <- roc_COP(max.x=max(df.COP$x), time.t=t, mod='exp_exp', params=params, n.quant=n.quant)
KM1.resauc <- data.frame()
for (i in t){
  res <- survivalROC(Stime=df.COP$t, status=df.COP$status, marker=df.COP$x, predict.time=i, method='KM')$AUC
  KM1.resauc <- rbind(KM1.resauc, c(i, res))
}
RiskROC1.resauc <- data.frame()
for (i in t){
  roc <- risksetROC(Stime=df.COP$t, status=df.COP$status, marker=df.COP$x, predict.time=i)$AUC
  RiskROC1.resauc <- rbind(RiskROC1.resauc, c(i,roc))
}
# =========================================
# Simulation study for Norm-Expo scenario
# =========================================
mu.x=5; sigma.x=0.8; lambda.t=0.2;
set.seed(123456)
params=list(mu.x=mu.x,lambda.t=lambda.t,c.rate=c.rate,beta=beta, sigma.x=sigma.x)
df.PH <- generate_dataPH('norm_exp', n.sample=data.sample, params=params)
t <- quantile(df.PH$t,probs=c(0.1,0.25,0.4,0.6,0.75,0.9))
ROC2_res <- roc_PH(max.x=max(df.PH$x), time.t=t, mod='norm_exp',
                         params=params, n.quant=n.quant, n.sample=monte.sample)
KM2.auc <- data.frame()
for (i in t){
  res <- survivalROC(Stime=df.PH$t, status=df.PH$status, marker=df.PH$x, predict.time=i, method='KM')$AUC
  KM2.auc <- rbind(KM2.auc, c(i, res))
}
RiskROC2.auc <- data.frame()
for (i in t){
  roc <- risksetROC(Stime=df.PH$t, status=df.PH$status, marker=df.PH$x, predict.time=i)$AUC
  RiskROC2.auc <- rbind(RiskROC2.auc, c(i,roc))
}
mu.x=5; sigma.x=0.8; lambda.t=0.2;
set.seed(123456)
params=list(mu.x=mu.x,lambda.t=lambda.t,c.rate=c.rate,theta=theta, sigma.x=sigma.x)
df.COP <- generate_dataCOP('norm_exp', n.sample=data.sample, params=params)
t <- quantile(df.COP$t,probs=c(0.1,0.25,0.4,0.6,0.75,0.9))
ROC2_resCOP <- roc_COP(max.x=max(df.COP$x), time.t=t, mod='norm_exp', params=params,n.quant=n.quant)
KM2.resauc <- data.frame()
for (i in t){
  res <- survivalROC(Stime=df.COP$t, status=df.COP$status, marker=df.COP$x, predict.time=i, method='KM')$AUC
  KM2.resauc <- rbind(KM2.resauc, c(i, res))
}
RiskROC2.resauc <- data.frame()
for (i in t){
  roc <- risksetROC(Stime=df.COP$t, status=df.COP$status, marker=df.COP$x, predict.time=i)$AUC
  RiskROC2.resauc <- rbind(RiskROC2.resauc, c(i,roc))
}

# =========================================
# Simulation study for Norm-Weib scenario
# =========================================
mu.x=5; sigma.x=0.8; shape.t=1.2; scale.t=1.4;
set.seed(123456)
params=list(mu.x=mu.x, shape.t=shape.t, scale.t=scale.t, c.rate=c.rate, beta=beta, sigma.x=sigma.x)
df.PH <- generate_dataPH('norm_weib', n.sample=data.sample, params=params)
t <- quantile(df.PH$t,probs=c(0.1,0.25,0.4,0.6,0.75,0.9))
ROC3_res <- roc_PH(max.x=max(df.PH$x), time.t=t, mod='norm_weib',
                          params=params, n.quant=n.quant, n.sample=monte.sample)
KM3.auc <- data.frame()
for (i in t){
  res <- survivalROC(Stime=df.PH$t, status=df.PH$status, marker=df.PH$x, predict.time=i, method='KM')$AUC
  KM3.auc <- rbind(KM3.auc, c(i, res))
}
RiskROC3.auc <- data.frame()
for (i in t){
  roc <- risksetROC(Stime=df.PH$t, status=df.PH$status, marker=df.PH$x, predict.time=i)$AUC
  RiskROC3.auc <- rbind(RiskROC3.auc, c(i,roc))
}
mu.x=5; sigma.x=0.8; shape.t=1.2; scale.t=1/1.4;
set.seed(123456)
params=list(mu.x=mu.x, sigma.x=sigma.x, shape.t=shape.t,
            scale.t=scale.t, c.rate=c.rate, theta=theta)
df.COP <- generate_dataCOP('norm_weib', n.sample=data.sample,params=params)
t <- quantile(df.COP$t,probs=c(0.1,0.25,0.4,0.6,0.75,0.9))
ROC3_resCOP <- roc_COP(max.x=max(df.COP$x), time.t=t,mod='norm_weib', params=params,n.quant=n.quant)
KM3.resauc <- data.frame()
for (i in t){
  res <- survivalROC(Stime=df.COP$t, status=df.COP$status, marker=df.COP$x, predict.time=i, method='KM')$AUC
  KM3.resauc <- rbind(KM3.resauc, c(i, res))
}
RiskROC3.resauc <- data.frame()
for (i in t){
  roc <- risksetROC(Stime=df.COP$t, status=df.COP$status, marker=df.COP$x, predict.time=i)$AUC
  RiskROC3.resauc <- rbind(RiskROC3.resauc, c(i,roc))
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
  labs(x=element_blank(),y=element_blank(), title='Normal-Weibull', lty="Sample size:") +
  theme(axis.text.y=element_blank(),legend.position=c(.85,0.2), legend.direction='vertical')

grid.arrange(p1,p2,p3,nrow=1,ncol=3,
             top = textGrob("Time-dependent ROC with various sample size:\n Copula Function",gp=gpar(fontsize=20,font=3)),
             bottom=textGrob("1-Specificity",gp=gpar(fontsize=20,font=3)),
             left=textGrob("Sensitivity",gp=gpar(fontsize=20,font=3), rot=90))

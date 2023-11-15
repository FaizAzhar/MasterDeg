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
PH1.auc <- ROC1_res %>% group_by(time) %>%
  summarise(AUC = AUC(1-specificity, sensitivity))
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
colnames(RiskROC1.auc) <- colnames(KM1.auc) <- c('time','AUC')

mu.x <- 1.4; lambda.t <- 0.2;
set.seed(123456)
params=list(mu.x=mu.x,lambda.t=lambda.t,c.rate=c.rate,theta=theta)
df.COP <- generate_dataCOP('exp_exp', n.sample=data.sample,params=params)
t <- quantile(df.COP$t,probs=c(0.1,0.25,0.4,0.6,0.75,0.9))
ROC1_resCOP <- roc_COP(max.x=max(df.COP$x), time.t=t, mod='exp_exp', params=params, n.quant=n.quant)
COP1.auc <- ROC1_resCOP %>% group_by(time) %>%
  summarise(AUC = AUC(1-specificity, sensitivity))
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
colnames(RiskROC1.resauc) <- colnames(KM1.resauc) <- c('time','AUC')

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
PH2.auc <- ROC2_res %>% group_by(time) %>%
  summarise(AUC = AUC(1-specificity, sensitivity))
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
colnames(RiskROC2.auc) <- colnames(KM2.auc) <- c('time','AUC')


mu.x=5; sigma.x=0.8; lambda.t=0.2;
set.seed(123456)
params=list(mu.x=mu.x,lambda.t=lambda.t,c.rate=c.rate,theta=theta, sigma.x=sigma.x)
df.COP <- generate_dataCOP('norm_exp', n.sample=data.sample, params=params)
t <- quantile(df.COP$t,probs=c(0.1,0.25,0.4,0.6,0.75,0.9))
ROC2_resCOP <- roc_COP(max.x=max(df.COP$x), time.t=t, mod='norm_exp', params=params,n.quant=n.quant)
COP2.auc <- ROC2_resCOP %>% group_by(time) %>%
  summarise(AUC = AUC(1-specificity, sensitivity))
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
colnames(RiskROC2.resauc) <- colnames(KM2.resauc) <- c('time','AUC')

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
PH3.auc <- ROC3_res %>% group_by(time) %>%
  summarise(AUC = AUC(1-specificity, sensitivity))
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
colnames(RiskROC3.auc) <- colnames(KM3.auc) <- c('time','AUC')

mu.x=5; sigma.x=0.8; shape.t=1.2; scale.t=1/1.4;
set.seed(123456)
params=list(mu.x=mu.x, sigma.x=sigma.x, shape.t=shape.t,
            scale.t=scale.t, c.rate=c.rate, theta=theta)
df.COP <- generate_dataCOP('norm_weib', n.sample=data.sample,params=params)
t <- quantile(df.COP$t,probs=c(0.1,0.25,0.4,0.6,0.75,0.9))
ROC3_resCOP <- roc_COP(max.x=max(df.COP$x), time.t=t,mod='norm_weib', params=params,n.quant=n.quant)
COP3.auc <- ROC3_resCOP %>% group_by(time) %>%
  summarise(AUC = AUC(1-specificity, sensitivity))
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
colnames(RiskROC3.resauc) <- colnames(KM3.resauc) <- c('time','AUC')

# ================
# ROC PH model
# ================
PH3.auc$mode <- PH2.auc$mode <- PH1.auc$mode <- 'PH'
COP3.auc$mode <- COP2.auc$mode <- COP1.auc$mode <- 'COP'
KM3.resauc$mode <- KM2.resauc$mode <- KM1.resauc$mode <- KM3.auc$mode <- KM2.auc$mode <- KM1.auc$mode <- 'KM'
RiskROC3.resauc$mode <- RiskROC2.resauc$mode <- RiskROC1.resauc$mode <- RiskROC3.auc$mode <- RiskROC2.auc$mode <- RiskROC1.auc$mode <- 'Risk Set'

ph1 <- rbind(PH1.auc,KM1.auc,RiskROC1.auc)
ph2 <- rbind(PH2.auc,KM2.auc,RiskROC2.auc)
ph3 <- rbind(PH3.auc,KM3.auc,RiskROC3.auc)

p1 <- ggplot(ph1, aes(x=time, y=AUC)) +
  geom_line(linewidth = 1, aes(linetype=mode)) +
  scale_linetype_manual(values=c('solid','dashed', 'dotdash'))+
  ylim(0.5,1)+
  labs(x='Time (years)',title='Exponential-Exponential')

p2 <- ggplot(ph2, aes(x=time, y=AUC)) +
  geom_line(linewidth = 1, aes(linetype=mode)) +
  scale_linetype_manual(values=c('solid','dashed', 'dotdash'))+
  ylim(0.5,1)+
  labs(x='Time (years)',title='Normal-Exponential')

p3 <- ggplot(ph3, aes(x=time, y=AUC)) +
  geom_line(linewidth = 1, aes(linetype=mode)) +
  scale_linetype_manual(values=c('solid','dashed', 'dotdash'))+
  ylim(0.5,1)+
  labs(x='Time (years)',title='Normal-Weibull')

grid.arrange(p1,p2,p3,nrow=3,ncol=1,
             top = textGrob("Comparison of AUC trend:\n PH model",gp=gpar(fontsize=20,font=3)))

# # ================
# # ROC Copula model
# # ================
cop1 <- rbind(COP1.auc,KM1.resauc,RiskROC1.resauc)
cop2 <- rbind(COP2.auc,KM2.resauc,RiskROC2.resauc)
cop3 <- rbind(COP3.auc,KM3.resauc,RiskROC3.resauc)

p1 <- ggplot(cop1, aes(x=time, y=AUC)) +
  geom_line(linewidth = 1, aes(linetype=mode)) +
  scale_linetype_manual(values=c('dotted','solid', 'dotdash'))+
  ylim(0.5,1)+
  labs(x='Time (years)',title='Exponential-Exponential')

p2 <- ggplot(cop2, aes(x=time, y=AUC)) +
  geom_line(linewidth = 1, aes(linetype=mode)) +
  scale_linetype_manual(values=c('dotted','solid', 'dotdash'))+
  ylim(0.5,1)+
  labs(x='Time (years)',title='Normal-Exponential')

p3 <- ggplot(cop3, aes(x=time, y=AUC)) +
  geom_line(linewidth = 1, aes(linetype=mode)) +
  scale_linetype_manual(values=c('dotted','solid', 'dotdash'))+
  ylim(0.5,1)+
  labs(x='Time (years)',title='Normal-Weibull')

grid.arrange(p1,p2,p3,nrow=3,ncol=1,
             top = textGrob("Comparison of AUC trend:\n Copula",gp=gpar(fontsize=20,font=3)))

pkgs <- c('openxlsx','copula','R2OpenBUGS','R2jags','coda','MCMCvis','readxl',
          'DescTools','dplyr','ggplot2','gridExtra','grid','devtools','latex2exp')
lapply(pkgs, require, character.only = TRUE)
load_all(".")

# =============================================
# Scripts below is used to study the impact of dependence measurement in PH model (Beta) & Copula function (Tau)
# towards the trend of AUC produced in the time-dependent ROC analysis.
# =============================================

n.quant=21; monte.sample=100000; data.sample=100; c.rate=0.2; beta <- c(0.8,1.1,1.4,1.7,2,2.4,3,4.3,4.6,5,5.4,5.7)
tau <- c(0.1,0.2,0.3,0.4117647059,0.5,0.6,0.7142857143,0.8,0.9, 0.95);
theta <- 2*tau/(1-tau)
# =========================================
# Simulation study for Expo-Expo scenario
# =========================================
AUC1_res <- NA
mu.x=1.4; lambda.t=0.2;
set.seed(123456)
for(i in 1:length(beta)){
  params=list(mu.x=mu.x,lambda.t=lambda.t,c.rate=c.rate,beta=beta[i])
  df.PH <- generate_dataPH('exp_exp', n.sample=data.sample,
                           params=params)
  t <- quantile(df.PH$t,probs=c(0.1,0.25,0.4,0.6,0.75,0.9))
  roc.exp_exp <- roc_PH(max.x=max(df.PH$x),time.t=t,mod='exp_exp',
                        params=params,n.quant=n.quant,n.sample=monte.sample)
  res <- roc.exp_exp %>% group_by(time) %>%
    summarize(AUC=AUC(1-specificity,sensitivity)) %>% data.frame()
  res <- cbind(res,rep(beta[i],length(t)))
  AUC1_res <- rbind(AUC1_res,res)
}

AUC1_resCOP <- NA
mu.x <- 1.4; lambda.t <- 0.2;
set.seed(123456)
for(i in 1:length(tau)){
  params=list(mu.x=mu.x,lambda.t=lambda.t,c.rate=c.rate,theta=theta[i])
  df.COP <- generate_dataCOP('exp_exp', n.sample=data.sample,params=params)
  t <- quantile(df.COP$t,probs=c(0.1,0.25,0.4,0.6,0.75,0.9))
  roc.exp_exp <- roc_COP(max.x=max(df.COP$x), time.t=t, mod='exp_exp', params=params, n.quant=n.quant)
  res <- roc.exp_exp %>% group_by(time) %>%
    summarize(AUC=AUC(1-specificity,sensitivity)) %>% data.frame()
  res <- cbind(res,rep(tau[i],length(t)))
  AUC1_resCOP <- rbind(AUC1_resCOP,res)
}


# =========================================
# Simulation study for Norm-Expo scenario
# =========================================
AUC2_res <- NA
mu.x=5; sigma.x=0.8; lambda.t=0.2;
set.seed(123456)
for(i in 1:length(beta)){
  params=list(mu.x=mu.x,lambda.t=lambda.t,c.rate=c.rate,beta=beta[i], sigma.x=sigma.x)
  df.PH <- generate_dataPH('norm_exp', n.sample=data.sample, params=params)
  t <- quantile(df.PH$t,probs=c(0.1,0.25,0.4,0.6,0.75,0.9))
  roc.norm_exp <- roc_PH(max.x=max(df.PH$x), time.t=t, mod='norm_exp',
                         params=params, n.quant=n.quant, n.sample=monte.sample)
  res <- roc.norm_exp %>% group_by(time) %>%
    summarize(AUC=AUC(1-specificity,sensitivity)) %>% data.frame()
  res <- cbind(res,rep(beta[i],length(t)))
  AUC2_res <- rbind(AUC2_res,res)
}

AUC2_resCOP <- NA
mu.x=5; sigma.x=0.8; lambda.t=0.2;
set.seed(123456)
for(i in 1:length(tau)){
  params=list(mu.x=mu.x,lambda.t=lambda.t,c.rate=c.rate,theta=theta[i], sigma.x=sigma.x)
  df.COP <- generate_dataCOP('norm_exp', n.sample=data.sample, params=params)
  t <-  quantile(df.COP$t,probs=c(0.1,0.25,0.4,0.6,0.75,0.9))
  roc.norm_exp <- roc_COP(max.x=max(df.COP$x), time.t=t, mod='norm_exp', params=params,n.quant=n.quant)
  res <- roc.norm_exp %>% group_by(time) %>%
    summarize(AUC=AUC(1-specificity,sensitivity)) %>% data.frame()
  res <- cbind(res,rep(tau[i],length(t)))
  AUC2_resCOP <- rbind(AUC2_resCOP,res)
}


# =========================================
# Simulation study for Norm-Weib scenario
# =========================================
AUC3_res <- NA
mu.x=5; sigma.x=0.8; shape.t=1.2; scale.t=1.4;
set.seed(123456)
for(i in 1:length(beta)){
  params=list(mu.x=mu.x, shape.t=shape.t, scale.t=scale.t, c.rate=c.rate, beta=beta[i], sigma.x=sigma.x)
  df.PH <- generate_dataPH('norm_weib', n.sample=data.sample, params=params)
  t <- quantile(df.PH$t,probs=c(0.1,0.25,0.4,0.6,0.75,0.9))
  roc.norm_weib <- roc_PH(max.x=max(df.PH$x), time.t=t, mod='norm_weib',
                          params=params, n.quant=n.quant, n.sample=monte.sample)
  res <- roc.norm_weib %>% group_by(time) %>%
    summarize(AUC=AUC(1-specificity,sensitivity)) %>% data.frame()
  res <- cbind(res,rep(beta[i],length(t)))
  AUC3_res <- rbind(AUC3_res,res)
}

AUC3_resCOP <- NA
mu.x=5; sigma.x=0.8; shape.t=1.2; scale.t=1/1.4;
set.seed(123456)
for(i in 1:length(tau)){
  params=list(mu.x=mu.x, sigma.x=sigma.x, shape.t=shape.t,
              scale.t=scale.t, c.rate=c.rate, theta=theta[i])
  df.COP <- generate_dataCOP('norm_weib', n.sample=data.sample,params=params)
  t <-  quantile(df.COP$t,probs=c(0.1,0.25,0.4,0.6,0.75,0.9))
  roc.norm_weib <- roc_COP(max.x=max(df.COP$x), time.t=t,mod='norm_weib', params=params,n.quant=n.quant)
  res <- roc.norm_weib %>% group_by(time) %>%
    summarize(AUC=AUC(1-specificity,sensitivity)) %>% data.frame()
  res <- cbind(res,rep(tau[i],length(t)))
  AUC3_resCOP <- rbind(AUC3_resCOP,res)
}

# ===============
# AUC and beta
# ===============
AUC1_res <- AUC1_res[-1,]
AUC2_res <- AUC2_res[-1,]
AUC3_res <- AUC3_res[-1,]
AUC1_res[,1] <- rep(c("t1","t2","t3","t4","t5","t6"),12)
AUC2_res[,1] <- rep(c("t1","t2","t3","t4","t5","t6"),12)
AUC3_res[,1] <- rep(c("t1","t2","t3","t4","t5","t6"),12)
colnames(AUC1_res) <- c('time','AUC','beta')
colnames(AUC2_res) <- c('time','AUC','beta')
colnames(AUC3_res) <- c('time','AUC','beta')

p1.add <- AUC1_res[which(AUC1_res$beta == 1.4),]
p1 <- ggplot(AUC1_res,aes(x = time, y = AUC1_res$AUC, group=as.factor(beta), color=as.factor(beta))) +
  geom_line(linewidth = 0.5) +
  geom_line(data=p1.add, aes(x = time, y = p1.add$AUC), color='red', size=1) +
  scale_colour_viridis_d() +
  ylim(0.65,1)+
  labs(x=element_blank(), y=element_blank(), title='Exponential-Exponential', color='Beta') +
  scale_x_discrete(labels=c(parse(text = TeX('$T_{10}$')), parse(text = TeX('$T_{25}$')), parse(text = TeX('$T_{40}$')), parse(text = TeX('$T_{60}$')), parse(text = TeX('$T_{75}$')), parse(text = TeX('$T_{90}$')))) +
  theme(legend.position=c(.965,0.6))

p2.add <- AUC2_res[which(AUC2_res$beta == 5),]
p2 <- ggplot(AUC2_res,aes(x = time, y = AUC2_res$AUC, group=as.factor(beta), color=as.factor(beta))) +
  geom_line(linewidth = 0.5) +
  geom_line(data=p2.add, aes(x = time, y = p2.add$AUC), color='red', size=1) +
  scale_colour_viridis_d() +
  ylim(0.65,1)+
  labs(x=element_blank(), y=element_blank(), title='Normal-Exponential', color='Beta') +
  scale_x_discrete(labels=c(parse(text = TeX('$T_{10}$')), parse(text = TeX('$T_{25}$')), parse(text = TeX('$T_{40}$')), parse(text = TeX('$T_{60}$')), parse(text = TeX('$T_{75}$')), parse(text = TeX('$T_{90}$')))) +
  theme(legend.position=c(.965,0.6))

p3.add <- AUC3_res[which(AUC3_res$beta == 5),]
p3 <- ggplot(AUC3_res,aes(x = time, y = AUC3_res$AUC, group=as.factor(beta), color=as.factor(beta))) +
  geom_line(linewidth = 0.5) +
  geom_line(data=p3.add, aes(x = time, y = p3.add$AUC), color='red', size=1) +
  scale_colour_viridis_d() +
  ylim(0.55,1)+
  labs(x=element_blank(), y=element_blank(), title='Normal-Weibull', color='Beta') +
  scale_x_discrete(labels=c(parse(text = TeX('$T_{10}$')), parse(text = TeX('$T_{25}$')), parse(text = TeX('$T_{40}$')), parse(text = TeX('$T_{60}$')), parse(text = TeX('$T_{75}$')), parse(text = TeX('$T_{90}$')))) +
  theme(legend.position=c(.965,0.6))

grid.arrange(p1,p2,p3,nrow=3,ncol=1,
             top = textGrob("AUC across time with various Beta\n but fixed Scale parameter",gp=gpar(fontsize=20,font=3)),
             bottom=textGrob("Time",gp=gpar(fontsize=20,font=3)),
             left=textGrob("AUC",gp=gpar(fontsize=20,font=3), rot=90))

# ===============
# AUC and tau
# ===============
AUC1_resCOP <- AUC1_resCOP[-1,]
AUC2_resCOP <- AUC2_resCOP[-1,]
AUC3_resCOP <- AUC3_resCOP[-1,]
AUC1_resCOP[,1] <- rep(c("t1","t2","t3","t4","t5","t6"),10)
AUC2_resCOP[,1] <- rep(c("t1","t2","t3","t4","t5","t6"),10)
AUC3_resCOP[,1] <- rep(c("t1","t2","t3","t4","t5","t6"),10)
colnames(AUC1_resCOP) <- c('time','AUC','tau')
colnames(AUC2_resCOP) <- c('time','AUC','tau')
colnames(AUC3_resCOP) <- c('time','AUC','tau')

p1.add <- AUC1_resCOP[which(AUC1_resCOP$tau == 0.4117647059),]
p1 <- ggplot(AUC1_resCOP,aes(x = time, y = AUC, group=as.factor(round(tau,1)), color=as.factor(round(tau,1)))) +
  geom_line(linewidth = 0.5) +
  geom_line(data=p1.add, aes(x = time, y = AUC), color='red', size=1) +
  scale_colour_viridis_d() +
  ylim(0.5,1)+
  labs(x=element_blank(), y=element_blank(), title='Exponential-Exponential', color='Tau') +
  scale_x_discrete(labels=c(parse(text = TeX('$T_{10}$')), parse(text = TeX('$T_{25}$')), parse(text = TeX('$T_{40}$')), parse(text = TeX('$T_{60}$')), parse(text = TeX('$T_{75}$')), parse(text = TeX('$T_{90}$')))) +
  theme(legend.position=c(.965,0.6))

p2.add <- AUC2_resCOP[which(AUC2_resCOP$tau == 0.7142857143),]
p2 <- ggplot(AUC2_resCOP,aes(x = time, y = AUC, group=as.factor(round(tau,1)), color=as.factor(round(tau,1)))) +
  geom_line(linewidth = 0.5) +
  geom_line(data=p2.add, aes(x = time, y = AUC), color='red', size=1) +
  scale_colour_viridis_d() +
  ylim(0.5,1)+
  labs(x=element_blank(), y=element_blank(), title='Normal-Exponential', color='Tau') +
  scale_x_discrete(labels=c(parse(text = TeX('$T_{10}$')), parse(text = TeX('$T_{25}$')), parse(text = TeX('$T_{40}$')), parse(text = TeX('$T_{60}$')), parse(text = TeX('$T_{75}$')), parse(text = TeX('$T_{90}$')))) +
  theme(legend.position=c(.965,0.6))

p3.add <- AUC3_resCOP[which(AUC3_resCOP$tau == 0.7142857143),]
p3 <- ggplot(AUC3_resCOP,aes(x = time, y = AUC, group=as.factor(round(tau,1)), color=as.factor(round(tau,1)))) +
  geom_line(linewidth = 0.5) +
  geom_line(data=p3.add, aes(x = time, y = AUC), color='red', size=1) +
  scale_colour_viridis_d() +
  ylim(0.5,1)+
  labs(x=element_blank(), y=element_blank(), title='Normal-Weibull', color='Tau') +
  scale_x_discrete(labels=c(parse(text = TeX('$T_{10}$')), parse(text = TeX('$T_{25}$')), parse(text = TeX('$T_{40}$')), parse(text = TeX('$T_{60}$')), parse(text = TeX('$T_{75}$')), parse(text = TeX('$T_{90}$')))) +
  theme(legend.position=c(.965,0.6))

grid.arrange(p1,p2,p3,nrow=3,ncol=1,
             top = textGrob("AUC across time with various Kendall's Tau\n but fixed Scale parameter",gp=gpar(fontsize=20,font=3)),
             bottom=textGrob("Time",gp=gpar(fontsize=20,font=3)),
             left=textGrob("AUC",gp=gpar(fontsize=20,font=3), rot=90))

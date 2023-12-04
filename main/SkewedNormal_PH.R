pkgs <- c('survivalROC','cubature','DescTools','dplyr','ggplot2','plotly',
          'gridExtra','grid','devtools','parTimeROC', 'sn', 'rstan', 'smoothROCtime')
lapply(pkgs, require, character.only = TRUE)

load_all(".")

muX = 6; sigmaX = 1.5; alphaX = 1;
muT = 1.8; sigmaT = 4; alphaT = 5
beta <- c(0.1, 0.5, 0.9)
num_sample <- c(35, 60, 90, 100)
# num_sample <- 100

# ================================================
# Sample size impact towards time-dependent ROC
# ================================================

params <- list(mu.x=muX, sigma.x=sigmaX, skew.x=alphaX,
               mu.t=muT, sigma.t=sigmaT, skew.t=alphaT,
               beta=beta[4], c.rate=0.1)

bayes.result <- list()

for(sam.size in num_sample){
  set.seed(42)
  simdata <- generate_dataPH(mod="snorm_snorm", params=params, n.sample=sam.size)

 # estimate parameters using Bayesian inference
 MCMCestimate_PH("snorm_snorm",simdata, f.path="C:/Users/farea/Documents/Rworkspace",
                 n.chains=1, n.iter=3000, n.burnin=500)
 sn_ph <- readRDS("snorm_snorm.rds")
 bayes.result <- append(bayes.result,sn_ph)
 fit_summary <- summary(sn_ph)$summary

 tt <- quantile(simdata$t, probs=c(0.5))

  # parametric time-dependent ROC curve
  ph.roc <- roc_PH(max(simdata$x), time.t=tt, mod="snorm_snorm",
                   params=list(mu.x=fit_summary[1], sigma.x=fit_summary[2],skew.x=fit_summary[3],
                               beta=fit_summary[4],mu.t=fit_summary[5], sigma.t=fit_summary[6], skew.t=fit_summary[7]),
                   n.quant=100,
                   low.lim=c(-15,-15),upp.lim=c(20,20))

  # nonparametric
  km.roc <- survivalROC(simdata$t, simdata$status, simdata$x, predict.time=tt, method='KM')

  # semiparametric
  st.roc <- stRoc(data=simdata[,c(2,3,1)], t=tt, tcr="C")

  # plotting
  jpeg(file=paste0("sim1_size",sam.size,".jpeg"))
  plot(0:1, 0:1, type='n', main=paste0('50th percentile of Time-to-event (Sample = ', sam.size,')'),
       xlab='1-Specificity', ylab='Sensitivity')
  lines(1-ph.roc$specificity, ph.roc$sensitivity, lty=1)
  lines(km.roc$FP, km.roc$TP, lty=2)
  lines(st.roc$FP, st.roc$TP, lty=3)
  legend('bottomright', inset=0.05, c("PH","KM","smoothTime"), lty=c(1,2,3))
  dev.off()
}

saveRDS(bayes.result,"bayes_result2.rds")
bayes <- readRDS("bayes_result.rds")
bayes.result <- append(bayes,bayes.result)
saveRDS(bayes.result,"bayes_result.rds")

for(i in 1:4){
  traceplot(bayes[[i]])
  plot(bayes[[i]])
  print(bayes[[i]])
}

# ================================================
# Varying beta impact towards time-dependent ROC
# ================================================

sim2 <- data.frame()
bayes.result <- list()

for(varbeta in beta){
  params <- list(mu.x=muX, sigma.x=sigmaX, skew.x=alphaX,
               mu.t=muT, sigma.t=sigmaT, skew.t=alphaT,
               beta=varbeta, c.rate=0.1)
  set.seed(42)
  simdata <- generate_dataPH(mod="snorm_snorm", params=params, n.sample=100)
  simdata <- na.omit(simdata)
  tt <- quantile(simdata$t, probs=c(0.2,0.4,0.5,0.6,0.8,0.9))

  # estimate parameters using Bayesian inference
  MCMCestimate_PH("snorm_snorm",simdata, f.path="C:/Users/farea/Documents/Rworkspace",
                  n.chains=1, n.iter=3000, n.burnin=500)
  sn_ph <- readRDS("snorm_snorm.rds")
  bayes.result <- append(bayes.result,sn_ph)
  fit_summary <- summary(sn_ph)$summary

  # parametric time-dependent ROC curve
  ph.roc <- roc_PH(max(simdata$x), time.t=tt, mod="snorm_snorm",
                   params=list(mu.x=fit_summary[1], sigma.x=fit_summary[2],skew.x=fit_summary[3],
                   beta=fit_summary[4],mu.t=fit_summary[5], sigma.t=fit_summary[6], skew.t=fit_summary[7]),
                   n.quant=100,
                   low.lim=c(-15,-15),upp.lim=c(15,15))
  ph.roc$beta <- varbeta
  sim2 <- rbind(sim2, ph.roc)
}

auc.roc <- sim2 %>% group_by(beta,time) %>% summarise(.,auc=AUC(1-specificity,sensitivity)) %>% data.frame()
auc.roc$time <- rep(c(1:6),3)
ggplot(auc.roc, aes(x = time, y = auc, linetype = as.factor(beta))) +
  geom_line(linewidth = 0.7) +
  theme_classic() +
  labs(x='time (quantile)',title = 'Fixed biomarker location and varying association level', lty="beta") +
  scale_x_discrete(limits=c("0.1","0.25","0.4","0.6","0.75","0.9"))

# ================================================
# Real data application
# ================================================

# Clean dataset
data(mayo)
mayo$time <- mayo$time/365
df <- mayo[,c(1,2,3)]
colnames(df) <- c("t","status","x")

# Parameter estimation
f.path = "C://Users/farea/Documents/Rworkspace"
MCMCestimate_PH('snorm_snorm', df, f.path=f.path, n.chains=1,n.iter=3000, n.burnin=500)
sn_ph <- readRDS("snorm_snorm.rds")
print(sn_ph)
traceplot(sn_ph)
plot(sn_ph)
fit_summary <- summary(sn_ph)$summary

# Producing time-dependent ROC curve and AUC
tt <- quantile(mayo$time, probs=c(0.2,0.4,0.5,0.6,0.8,0.9))
ph.roc <- roc_PH(max(mayo$mayoscore5), time.t=tt, mod="snorm_snorm",
                 params=list(mu.x=fit_summary[1], sigma.x=fit_summary[2],skew.x=fit_summary[3], beta=fit_summary[4],mu.t=fit_summary[5], sigma.t=fit_summary[6], skew.t=fit_summary[7]),
                 n.quant=100,
                 low.lim=c(-100,-100),upp.lim=c(100,100))
sum(is.na(ph.roc$specificity));sum(is.na(ph.roc$sensitivity))

# parametric
auc.ph <- ph.roc %>% group_by(time) %>% summarise(AUC=AUC(1-specificity, sensitivity)) %>% data.frame()
auc.ph$mod <- 'PH'

# nonparametric
auc.km <- data.frame()
for (i in tt){
  res <- survivalROC(Stime=mayo$time, status=mayo$censor, marker=mayo$mayoscore5, predict.time=i, method='KM')$AUC
  auc.km <- rbind(auc.km, c(i, res))
}
colnames(auc.km) <- c('time','AUC')
auc.km$mod <- 'KM'

# semiparametric
auc.smooth <- data.frame()
df$status <- as.numeric(df$status)

for (i in tt){
  res <- mean(as.numeric(stRoc(data=df, t=i, tcr="C")$auc))
  auc.smooth <- rbind(auc.smooth, c(i, res))
}
colnames(auc.smooth) <- c('time','AUC')
auc.smooth$mod <- 'smoothTime'

aucplot <- rbind(auc.ph,auc.km,auc.smooth)
ggplot(aucplot, aes(x=time, y=AUC, linetype = mod)) +
  geom_line(linewidth = 0.7) +
  labs(title = "AUC trend for time-dependent ROC", x = 'Year', linetype='Model') +
  ylim(0.5,1) +
  theme_classic()

# devtools::install_github('FaizAzhar/MasterDeg') # install package from github repository
pkgs <- c('survivalROC','cubature','DescTools','dplyr','ggplot2','plotly','TH.data',
          'gridExtra','grid','devtools','parTimeROC', 'sn', 'rstan', 'smoothROCtime',
          'moments')
lapply(pkgs, require, character.only = TRUE)

load_all(".")

# # Setting used to check affect of sample size (Sim No.1) (Comment if not run in this round)
# muX = 0; sigmaX = 1; alphaX = 2;
# muT = 0; sigmaT = 1; alphaT = 2
# beta = 0.5
# # num_sample <- c(35, 60, 90)
# num_sample <- 100
# # ================================================
# # Sample size impact towards time-dependent ROC
# # ================================================
#
# params <- list(mu.x=muX, sigma.x=sigmaX, skew.x=alphaX,
#                mu.t=muT, sigma.t=sigmaT, skew.t=alphaT,
#                beta=beta, c.rate=0.1)
#
# bayes.result <- list()
#
# for(sam.size in num_sample){
#   set.seed(42)
#   simdata <- generate_dataPH(mod="snorm_snorm", params=params, n.sample=sam.size)
#
#   # estimate parameters using Bayesian inference
#   MCMCestimate_PH("snorm_snorm",simdata, f.path="C:/Users/farea/Documents/Rworkspace",
#                   n.chains=2, n.iter=2000, n.burnin=500)
#   sn_ph <- readRDS("snorm_snorm.rds")
#   bayes.result <- append(bayes.result,sn_ph)
#   fit_summary <- summary(sn_ph)$summary
#
#   tt <- quantile(simdata$t, probs=c(0.5))
#
#  # parametric time-dependent ROC curve
#  ph.roc <- roc_PH(max(simdata$x), min.x=min(simdata$x), time.t=tt, mod="snorm_snorm",
#                   params=list(mu.x=fit_summary[1], sigma.x=fit_summary[2],skew.x=fit_summary[3],
#                               beta=fit_summary[4],mu.t=fit_summary[5], sigma.t=fit_summary[6], skew.t=fit_summary[7]),
#                   n.quant=100,
#                   low.lim=c(-5,-5),upp.lim=c(5,5))
#
#  # nonparametric
#  km.roc <- survivalROC(simdata$t, simdata$status, simdata$x, predict.time=tt, method='KM')
#
#  # semiparametric
#  simdata$status <- as.numeric(simdata$status)
#  simdata$x <- as.numeric(simdata$x)
#  simdata$t <- as.numeric(simdata$t)
#  st.roc <- stRoc(data=simdata[,c(2,3,1)], t=tt, tcr="C")
#
#  # plotting
#  jpeg(file=paste0("sim2_size",sam.size,".jpeg"))
#  plot(0:1, 0:1, type='n', main=paste0('50th percentile of Time-to-event (Sample = ', sam.size,')'),
#       xlab='1-Specificity', ylab='Sensitivity')
#  lines(1-ph.roc$specificity, ph.roc$sensitivity, lty=1)
#  lines(km.roc$FP, km.roc$TP, lty=2)
#  lines(st.roc$FP, st.roc$TP, lty=3)
#  legend('bottomright', inset=0.05, c(paste0("PH (AUC = ", round(AUC(1-ph.roc$specificity,ph.roc$sensitivity),4),")"),
#                                      paste0("KM (AUC = ", round(km.roc$AUC,4),")"),
#                                      paste0("smoothTime (AUC = ",round(as.numeric(st.roc$auc[1]),4),")")), lty=c(1,2,3))
#  dev.off()
# }
#
# # saveRDS(bayes.result,"bayes_result.rds")
# bayes <- readRDS("bayes_result.rds")
# # bayes.result <- append(bayes,bayes.result)
# # saveRDS(bayes.result,"bayes_result.rds")
#
# for(i in 1:4){
#   traceplot(bayes[[i]])
#   plot(bayes[[i]])
#   print(bayes[[i]])
# }

# # ================================================
# # Bootstrapping random data generation procedure
# # ================================================
# library(moments)
#
# ## Setting used to check affect of varying beta (Sim No.2)
# bootstrap.skew <- list()
#
# muX = 0; sigmaX = 1; alphaX = -1;
# muT = 0; sigmaT = 1; alphaT = -1
# beta <- c(-0.9, -0.5, 0.5, 0.9)
# num_sample <- 60
# iter <- 500
#
# # Start bootstrapping
# sim2 <- data.frame()
#
# for(varbeta in beta){
#   params <- list(mu.x=muX, sigma.x=sigmaX, skew.x=alphaX,
#                  mu.t=muT, sigma.t=sigmaT, skew.t=alphaT,
#                  beta=varbeta, c.rate=0.1)
#   for(boot.s in 1:iter){
#     tryCatch({
#       simdata <- generate_dataPH(mod="snorm_snorm", params=params, n.sample=num_sample)
#       sim2 <- rbind(sim2,list(x=skewness(simdata$x), t=skewness(simdata$t), b=varbeta))
#       print(paste0("iter = ",boot.s," ,uncensor = ",skewness(simdata$t)))
#     }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
#
#   }
# }
#
# skewww <- sim2 %>% na.omit() %>% group_by(b) %>% summarise(skewness.x=mean(x), skewness.t = mean(t))
# hist(simdata$t, freq = F)
# bootstrap.skew <- append(bootstrap.skew,skewww)
# saveRDS(bootstrap.skew,"bootstrap_skew.rds")

# # ================================================
# # Varying beta impact towards AUC
# # ================================================
# muX = 0; sigmaX = 1; alphaX = 1;
# muT = 0; sigmaT = 1; alphaT = 1
# beta <- c(0.9, 0.5, -0.5, -0.9)
# num_sample <- 60
#
# # Start simulation
# sim2 <- data.frame()
#
# for(varbeta in beta){
#   params <- list(mu.x=muX, sigma.x=sigmaX, skew.x=alphaX,
#                mu.t=muT, sigma.t=sigmaT, skew.t=alphaT,
#                beta=varbeta, c.rate=0.1)
#   set.seed(42)
#   simdata <- generate_dataPH(mod="snorm_snorm", params=params, n.sample=num_sample)
#   print(paste0("X skewness = ", skewness(simdata$x), " min.x = ",min(simdata$x), " max.x = ",max(simdata$x)))
#   print(paste0("T skewness = ", skewness(simdata$x), " min.t = ",min(simdata$t), " max.t = ",max(simdata$t)))
#
#   tt <- quantile(df$t, probs=c(0.2,0.4,0.5,0.6,0.8,0.9))
#
# # parametric time-dependent ROC curve
# ph.roc <- roc_PH(max.x=max(simdata$x), min.x=min(simdata$x), time.t=tt, mod="snorm_snorm",
#                  params=list(mu.x=muX, sigma.x=sigmaX ,skew.x=alphaX,
#                  beta=varbeta ,mu.t=muT, sigma.t=sigmaT, skew.t=alphaT),
#                  n.quant=30,
#                  low.lim=c(-4,-4),upp.lim=c(4,4))
#
#   ph.roc$beta <- varbeta
#   sim2 <- rbind(sim2, ph.roc)
#
# }
#
# auc.roc <- sim2 %>% group_by(beta,time) %>% summarise(.,auc=AUC(1-specificity,sensitivity)) %>% data.frame()
# auc.roc$time <- rep(c(1:6),length(beta))
# ggplot(auc.roc, aes(x = time, y = auc, linetype = as.factor(beta))) +
#   geom_line(linewidth = 0.7) +
#   scale_linetype_manual(values = c("solid","dotted","dotdash","longdash"),breaks = c('0.9','0.5','-0.5','-0.9')) +
#   theme_classic() +
#   labs(x='time (quantile)',title = 'Varying association level (Skew X = 1, Skew T = 1)', lty="beta") +
#   scale_x_discrete(limits=c("0.1","0.25","0.4","0.6","0.75","0.9"))

# # # ================================================
# # # Real data application (Comment if not used in this round)
# # # ================================================
#
# # Clean dataset
# data(wpbc)
# df <- wpbc[,c(2,1,6)]
# df <- na.omit(df)
# colnames(df) <- c("t","status","x")
# df$status <- ifelse(df$status == 'R', 1, 0)
# XX <- max(df$t)
# TT <- max(df$x)
# df$t <- df$t/max(df$t)
# df$x <- df$x/max(df$x)
# print(paste0('min x = ', min(df$x),
#              ', max x = ', max(df$x),
#              ', min t = ', min(df$t),
#              ', max t = ', max(df$t)))
#
# s <- subplot(
#   plot_ly(x = df$x, type = "histogram"),
#   plotly_empty(),
#   plot_ly(x = df$x, y = df$t, type = "histogram2dcontour"),
#   plot_ly(y = df$t, type = "histogram"),
#   nrows = 2, heights = c(0.2, 0.8), widths = c(0.8, 0.2), margin = 0,
#   shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE
# )
# fig <- layout(s, showlegend = FALSE)
#
# fig
#
# # Parameter estimation
# f.path = "C://Users/farea/Documents/Rworkspace"
# MCMCestimate_PH('snorm_snorm', df, f.path=f.path, n.chains=2,n.iter=2000, n.burnin=500)
# sn_ph <- readRDS("snorm_snorm.rds")
# print(sn_ph)
# traceplot(sn_ph)
# plot(sn_ph)
# fit_summary <- summary(sn_ph)$summary
#
#
# # Producing time-dependent ROC curve and AUC
# tt <- quantile(df$t, probs=c(0.2,0.4,0.5,0.6,0.8,0.9))
# ph.roc <- roc_PH(max(df$x), min.x = min(df$x), time.t=tt, mod="snorm_snorm",
#                  params=list(mu.x=fit_summary[1], sigma.x=fit_summary[2],skew.x=fit_summary[3], beta=fit_summary[4],mu.t=fit_summary[5], sigma.t=fit_summary[6], skew.t=fit_summary[7]),
#                  n.quant=100,
#                  low.lim=c(-3,-3),upp.lim=c(3,3))
# sum(is.na(ph.roc$specificity));sum(is.na(ph.roc$sensitivity))
#
# # parametric
# auc.ph <- ph.roc %>% group_by(time) %>% summarise(AUC=AUC(1-specificity, sensitivity)) %>% data.frame()
# auc.ph$mod <- 'PH'
#
# # nonparametric
# auc.km <- data.frame()
# km.roc <- data.frame()
# for (i in tt){
#   store <- survivalROC(Stime=df$t, status=df$status, marker=df$x, predict.time=i, method='KM')
#   km.roc <- rbind(km.roc, data.frame(FP=store$FP,TP=store$TP,time=rep(i,length(store$FP))))
#   res <- store$AUC
#   auc.km <- rbind(auc.km, c(i, res))
# }
# colnames(auc.km) <- c('time','AUC')
# auc.km$mod <- 'KM'
#
# # semiparametric
# auc.smooth <- data.frame()
# df$status <- as.numeric(df$status)
# smooth.roc <- data.frame()
#
# for (i in tt){
#   store <- stRoc(data=df, t=i, tcr="C")
#   smooth.roc <- rbind(smooth.roc, data.frame(FP=store$FP,TP=store$TP,time=rep(i,length(store$FP))))
#   res <- mean(as.numeric(store$auc))
#   auc.smooth <- rbind(auc.smooth, c(i, res))
# }
# colnames(auc.smooth) <- c('time','AUC')
# auc.smooth$mod <- 'smoothTime'
#
# aucplot <- rbind(auc.ph,auc.km,auc.smooth)
# ggplot(aucplot, aes(x=time, y=AUC, linetype = mod)) +
#   geom_line(linewidth = 0.7) +
#   scale_linetype_manual(values=c("dotted","solid","dotdash")) +
#   labs(title = "AUC trend for time-dependent ROC", x = 'Normalized Month', linetype='Model') +
#   ylim(0.3,1) +
#   theme_classic()

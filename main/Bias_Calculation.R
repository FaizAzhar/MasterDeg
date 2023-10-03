PHpar1 <- c(1.4, 0.2, 1.2)
PHpar2 <- c(5, 0.2, 1.2,0.8)
PHpar3 <- c(5, 1.4, 1.2, 1.2, 0.8)
COPpar1 <- c(1.4, 0.2, 2*0.7/0.3)
COPpar2 <- c(5, 0.2, 2*0.7/0.3, 0.8)
# careful.. mu.x and shape.t switch in previous run
COPpar3 <- c(1.2, 5, 2*0.7/0.3, 0.8, 1.4)

PHbias1 <- sweep(store.param1,2,PHpar1)
# PHbias1$Average <- (PHbias1[,1] + PHbias1[,2] + PHbias1[,3])/3
PHbias1$Size <- c(30, 60, 120, 200)
PHbias1$Model <- 'Exponential-Exponential'
PHbias1$Method <- 'PH'
PHbias2 <- sweep(store.param2,2,PHpar2)
# PHbias2$Average <- (PHbias2[,1] + PHbias2[,2] + PHbias2[,3] + PHbias2[,4])/4
PHbias2$Size <- c(30, 60, 120, 200)
PHbias2$Model <- 'Normal-Exponential'
PHbias2$Method <- 'PH'
PHbias3 <- sweep(store.param3,2,PHpar3)
# PHbias3$Average <- (PHbias3[,1] + PHbias3[,2] + PHbias3[,3] + PHbias3[,4] + PHbias3[,5])/5
PHbias3$Size <- c(30, 60, 120, 200)
PHbias3$Model <- 'Normal-Weibull'
PHbias3$Method <- 'PH'
PHbias <- rbind(PHbias1[,c('beta','Size','Model','Method')], PHbias2[,c('beta','Size','Model','Method')],
                PHbias3[,c('beta','Size','Model','Method')])

COPbias1 <- sweep(store.paramCOP1,2,COPpar1)
# COPbias1$Average <- (COPbias1[,1] + COPbias1[,2] + COPbias1[,3])/3
COPbias1$Size <- c(30, 60, 120, 200)
COPbias1$Model <- 'Exponential-Exponential'
COPbias1$Method <- 'COP'
COPbias2 <- sweep(store.paramCOP2,2,COPpar2)
# COPbias2$Average <- (COPbias2[,1] + COPbias2[,2] + COPbias2[,3] + COPbias2[,4])/4
COPbias2$Size <- c(30, 60, 120, 200)
COPbias2$Model <- 'Normal-Exponential'
COPbias2$Method <- 'COP'
COPbias3 <- sweep(store.paramCOP3,2,COPpar3)
# COPbias3$Average <- (COPbias3[,1] + COPbias3[,2] + COPbias3[,3] + COPbias3[,4] + COPbias3[,5])/5
COPbias3$Size <- c(30, 60, 120, 200)
COPbias3$Model <- 'Normal-Weibull'
COPbias3$Method <- 'COP'
COPbias <- rbind(COPbias1[,c('theta','Size','Model', 'Method')], COPbias2[,c('theta','Size','Model', 'Method')],
                 COPbias3[,c('theta','Size','Model', 'Method')])

# overallbias <- rbind(PHbias, COPbias)

p1 <- ggplot(PHbias, aes(x=Size, y=beta)) +
  geom_line(linewidth = 0.8, aes(linetype = Model)) +
  scale_linetype_manual(values=c('dotted','dotdash','solid')) +
  ylim(-1,1) +
  labs(x=element_blank(), y=element_blank(), title='Proportional Hazard Model (Bias of Beta)') +
  theme(axis.text.x=element_blank(),legend.position=c(.75,0.8), legend.direction='horizontal')

p2 <- ggplot(COPbias, aes(x=Size, y=theta)) +
  geom_line(linewidth = 0.8, aes(linetype = Model), show.legend=FALSE) +
  scale_linetype_manual(values=c('dotted','dotdash','solid')) +
  ylim(-4.7,0.5) +
  labs(x=element_blank(), y=element_blank(), title='Copula function (Bias of Alpha)', lty='Model')

grid.arrange(p1,p2,nrow=2,ncol=1,
             bottom=textGrob("Sample Size",gp=gpar(fontsize=20,font=3)),
             left=textGrob("Bias",gp=gpar(fontsize=20,font=3), rot=90))

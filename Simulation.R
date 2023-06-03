library(survivalROC) # for PBC data
library(cubature) # to perform integration
library(VineCopula) # to calculate copula CDF
library(MASS) # for simulation of Bivariate Normal
library(DescTools) # to use trapezoid rule for AUC calculation

#=================================================
#  Function to compute sensitivity and specificity 
#  for method 1 (hazard function)
#=================================================

#   c=marker measurement, t=time, 
#   mod=function of joint model

haz.sens <- function(c,t,mod,par){
  prob.t <- cuhre(f=mod,
                  pars=par,
                  lower=rep(0,2),
                  upper=c(Inf,t))$integral
  prob.xt <- cuhre(f=mod,
                  pars=par,
                  lower=c(c,0),
                  upper=c(Inf,t))$integral
  
  return(prob.xt/prob.t)
}

haz.spec <- function(c,t,mod,par){
  prob.t <- cuhre(f=mod,
                  pars=par,
                  lower=c(0,t),
                  upper=c(Inf,Inf))$integral
  prob.xt <- cuhre(f=mod,
                   pars=par,
                   lower=c(0,t),
                   upper=c(c,Inf))$integral
  
  return(prob.xt/prob.t)
}

# =================================================
#  Function to compute sensitivity and specificity 
#  for method 2 (copula function)
# =================================================

#   c=marker measurement, t=time, family=refer VineCopula to choose copula
#   pars=parameter for selected copula

cop.sens <- function(c,t, family, pars){
  prob.xt <- BiCopCDF(c,t,family=family,
                      par = pars[1], par2 = ifelse(is.na[pars[2]],pars[2],0))
  return((t - prob.xt)/t)
}

cop.spec <- function(c,t, family, pars){
  prob.xt <- BiCopCDF(c,t,family=family,
                      par = pars[1], par2 = ifelse(is.na[pars[2]],pars[2],0))
  return((c - prob.xt)/(1-t))
}

# =================================================
#  Joint Model function
# =================================================

#     z[1]=marker measurement, z[2]=time

# expo-expo
joint.mod1 <- function(z,pars) return(z[1]*pars[1]*exp(-(z[1]*pars[1]*z[2]))*pars[2]*exp(-(pars[2]*z[1])))


# =================================================
#  log-likelihood function
# =================================================

# xi=marker measurement, ti=time, 
# ui=quantile of marker, vi=quantile of time

cop.mod.ll <- function(par,Ui,Vi,ui,vi){
  # par[1] = theta
  n <- length(Ui)
  logll <- n*log(1+par[1])-(2+1/par[1])*(sum(log((1-Ui)**(-par[1])+Vi**(-par[1])-1)))-(par[1]+1)*(sum(log((1-Ui)*Vi)))-sum(log(ui))+sum(log(vi))
  return(-logll)
}

#    expo-expo
haz.mod1.ll <- function(xi,ti,par){
  # par[1] = beta, par[2] = mu
  n <- length(xi)
  logll <- n*log(par[1]) - par[1]*(sum(xi*ti)) + n*log(par[2]) - par[2]*sum(xi)
  return(-logll)
}

marg.mod1T.ll <- function(ti,par){
  # par[1] = lambda
  n <- length(ti)
  logll <- n*log(par[1]) - par[1]*(sum(ti))
  return(-logll)
}

marg.mod1X.ll <- function(xi,par){
  # par[1] = lambda
  n <- length(xi)
  logll <- n*log(par[1]) - par[1]*(sum(xi))
  return(-logll)
}

#    normal-expo
haz.mod2.ll <- function(par,xi,ti){
  # par[1] = beta, par[2] = mu, par[3] = sigma
  n <- length(xi)
  logll <- n*log(par[1])+sum(xi)-par[1]*sum(xi*ti)-(n/2)*log(par[3]**2)-(1/(2*par[3]**2))*sum((xi-par[2])**2)
  return(-logll)
}

marg.mod2T.ll <- function(ti,par){
  # par[1] = lambda
  n <- length(ti)
  logll <- n*log(par[1]) - par[1]*(sum(ti))
  return(-logll)
}

marg.mod2X.ll <- function(xi,par){
  # par[1] = mu, par[2] = sigma
  n <- length(xi)
  logll <- (-n/2)*log((2*pi*par[2]**2))-sum((xi-par[1])**2)/(2*par[2]**2)
  return(-logll)
}

#    normal-weibull
haz.mod3.ll <- function(par,xi,ti){
  # par[1] = beta, par[2] = gamma, par[3] = mu, par[4] = sigma
  n <- length(xi)
  logll <- n*log(par[2])+(par[2]-1)*sum(log(xi*ti))+n*par[2]*log(par[1])-par[1]**par[2]*sum((xi*ti)**par[2])-(n/2)*log(par[4]**2)-(1/(2*par[4]**2))*sum((xi-par[3])**2)
  return(-logll)
}

marg.mod3T.ll <- function(ti,par){
  # par[1] = lambda, par[2] = gamma
  n <- length(ti)
  logll <- n*log(par[1]*par[2])+(par[2]-1)*sum(log(ti*par[1]))-sum((ti*par[1])**par[2])
  return(-logll)
}

marg.mod3X.ll <- function(xi,par){
  # par[1] = mu, par[2] = sigma
  n <- length(xi)
  logll <- (-n/2)*log((2*pi*par[2]**2))-sum((xi-par[1])**2)/(2*par[2]**2)
  return(-logll)
}

# =======================================
#  Function for simulation
# =======================================

sim.run <- function(nrun, sample_size=100, mu, sigma){
 
  set.seed(24)
  sim.data <- as.data.frame(mvrnorm(n=sample_size, mu=mu, Sigma=sigma))
  mod1.par.df <- data.frame()
  mod2.par.df <- data.frame()
  mod3.par.df <- data.frame()
  
  for(run in 1:nrun){
    boot.data <- sim.data[sample(nrow(sim.data), sample_size, replace=T), ]

    # Fitting for expo-expo
    #     hazard function
    haz.mod1.fit.par <- optim(c(0.1,0.1),haz.mod1.ll, xi=boot.data[,1], ti=boot.data[,2], 
                              hessian=T,method='L-BFGS-B',
                              lower=c(0.1,0.1), upper=c(999,999))
    haz.mod1.par <- haz.mod1.fit.par$par
    
    #     copula function
    marg.mod1X.fit.par <- optim(0.1, marg.mod1X.ll, xi=boot.data[,1],
                                method='Brent',hessian=T,
                                lower=0, upper=999) #stage 1 ~ MLE of marginals
    marg.mod1T.fit.par <- optim(0.1, marg.mod1T.ll, ti=boot.data[,2],
                                method='Brent',hessian=T,
                                lower=0, upper=999)
    marg.mod1.par <- c(marg.mod1X.fit.par$par,marg.mod1T.fit.par$par)
    U.i <- pexp(boot.data[,1], rate=marg.mod1.par[1]) #stage 2 ~ use estimated marginal for MLE of copula
    V.i <- pexp(boot.data[,2], rate=marg.mod1.par[2])
    u.i <- dexp(boot.data[,1], rate=marg.mod1.par[1])
    v.i <- dexp(boot.data[,2], rate=marg.mod1.par[2])
    cop.mod1.fit.par <- optim(2, cop.mod.ll, Ui=U.i, Vi=V.i, ui=u.i, vi=v.i,
                         hessian=T)
    
    #    storing all parameters for each run
    mod1.par.df <- rbind(mod1.par.df,c(haz.mod1.par[1],haz.mod1.par[2],
                                       marg.mod1.par[1],marg.mod1.par[2],
                                       cop.mod1.fit.par[1]))
    
    # Fitting for normal-expo
    #     hazard function
    haz.mod2.fit.par <- optim(rep(0.1,3),haz.mod2.ll, xi=boot.data[,1], ti=boot.data[,2],
                              hessian=T,method='L-BFGS-B',
                              lower=c(0.1,-999,0.1), upper=rep(999,3))
    haz.mod2.par <- haz.mod2.fit.par$par
    
    #     copula function
    marg.mod2X.fit.par <- optim(rep(0.1,2), marg.mod2X.ll, xi=boot.data[,1],
                                method='L-BFGS-B',hessian=T,
                                lower=c(-999,0.1), upper=rep(999,2)) #stage 1 ~ MLE of marginals
    marg.mod2T.fit.par <- optim(0.1, marg.mod2T.ll, ti=boot.data[,2],
                                method='Brent',hessian=T,
                                lower=0, upper=999)
    marg.mod2.par <- c(marg.mod2X.fit.par$par,marg.mod2T.fit.par$par)
    U.i <- pnorm(boot.data[,1], mean=marg.mod2.par[1], sd=marg.mod2.par[2]) #stage 2 ~ use estimated marginal for MLE of copula
    V.i <- pexp(boot.data[,2], rate=marg.mod2.par[3])
    u.i <- dnorm(boot.data[,1], mean=marg.mod2.par[1], sd=marg.mod2.par[2])
    v.i <- dexp(boot.data[,2], rate=marg.mod2.par[3])
    cop.mod2.fit.par <- optim(2, cop.mod.ll, Ui=U.i, Vi=V.i, ui=u.i, vi=v.i,
                              hessian=T)
    #    storing all parameters for each run
    mod2.par.df <- rbind(mod2.par.df,c(haz.mod2.par[1],haz.mod2.par[2],haz.mod2.par[3],
                                       marg.mod2.par[1],marg.mod2.par[2],marg.mod2.par[3],
                                       cop.mod2.fit.par[1]))

    # Fitting for normal-weibull
    #     hazard function
    haz.mod3.fit.par <- optim(rep(0.1,4),haz.mod3.ll, xi=boot.data[,1], ti=boot.data[,2],
                              hessian=T)
    haz.mod3.par <- haz.mod3.fit.par$par
    
    #     copula function
    marg.mod3X.fit.par <- optim(rep(0.1,2), marg.mod3X.ll, xi=boot.data[,1],
                                method='L-BFGS-B',hessian=T,
                                lower=c(-999,0.1), upper=rep(999,2)) #stage 1 ~ MLE of marginals
    marg.mod3T.fit.par <- optim(rep(0.1,2), marg.mod3T.ll, ti=boot.data[,2],
                                hessian=T)
    marg.mod3.par <- c(marg.mod3X.fit.par$par,marg.mod3T.fit.par$par)
    U.i <- pnorm(boot.data[,1], mean=marg.mod3.par[1], sd=marg.mod3.par[2]) #stage 2 ~ use estimated marginal for MLE of copula
    V.i <- pweibull(boot.data[,2], scale=marg.mod3.par[3], shape=marg.mod3.par[4])
    u.i <- dnorm(boot.data[,1], mean=marg.mod3.par[1], sd=marg.mod3.par[2])
    v.i <- dweibull(boot.data[,2], scale=marg.mod3.par[3], shape=marg.mod3.par[4])
    cop.mod3.fit.par <- optimize(cop.mod.ll, c(0.1,50),
                                 Ui=U.i, Vi=V.i, ui=u.i, vi=v.i)
    
    #    storing all parameters for each run
    mod3.par.df <- rbind(mod3.par.df,c(haz.mod3.par[1],haz.mod3.par[2],haz.mod3.par[3],haz.mod3.par[4],
                                       marg.mod3.par[1],marg.mod3.par[2],marg.mod3.par[3],marg.mod3.par[4],
                                       cop.mod3.fit.par[1]))
  }
  
  colnames(mod1.par.df) <- c('haz.beta','haz.mu','marg.lambda1','marg.lambda2','cop.rho')
  colnames(mod2.par.df) <- c('haz.beta','haz.mu','haz.sigma','marg.mu','marg.sigma','marg.lambda','cop.rho')
  colnames(mod3.par.df) <- c('haz.beta','haz.gamma','haz.mu','haz.sigma','marg.mu','marg.sigma','marg.lambda','marg.gamma','cop.rho')
  
  # calculating bootstrap S.E and CI
  mod1 <- list('Fitted.par'=apply(mod1.par.df,2,mean),
               'CI.lower.upper'=apply(mod1.par.df,2,quantile, probs=c(0.05,0.95)),
               'S.Err'=sqrt(apply(mod1.par.df,2,var)))
  mod2 <- list('Fitted.par'=apply(mod2.par.df,2,mean),
               'CI.lower.upper'=apply(mod2.par.df,2,quantile, probs=c(0.05,0.95)),
               'S.Err'=sqrt(apply(mod2.par.df,2,var)))
  mod3 <- list('Fitted.par'=apply(mod3.par.df,2,mean),
               'CI.lower.upper'=apply(mod3.par.df,2,quantile, probs=c(0.05,0.95)),
               'S.Err'=sqrt(apply(mod3.par.df,2,var)))
  return(list('expo-expo'=mod1,'normal-expo'=mod2,'normal-weibull'=mod3))
}

sim.run(1,sample_size=100,mu=c(5,7),sigma=matrix(c(1,-.9,-.9,1),ncol=2))

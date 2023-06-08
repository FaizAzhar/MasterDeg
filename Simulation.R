library(survivalROC) # for PBC data
library(cubature) # to perform integration
library(copula) # to calculate copula CDF
library(MASS) # for simulation of Bivariate Normal
library(DescTools) # to use trapezoid rule for AUC calculation
library(psych)
library(dplyr)
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

#    expo-expo
haz.mod1.ll <- function(xi,ti,par){
  n <- length(xi)
  beta <- par[1]
  lambda <- par[2]
  logll <- n*log(beta) - beta*(sum(xi*ti)) + n*log(lambda) - lambda*sum(xi)
  return(-logll)
}

cop.mod1.ll <- function(pars,xi,ti){
  theta <- pars[1]
  mu    <- exp(-pars[2])
  lambda <- exp(-pars[3])
  rot.cop <- rotCopula(claytonCopula(theta),flip=c(F,T)) # rotated 90 degrees
  f.X <- dexp(xi, rate=mu) #first stage ~ estimate marginal
  f.T <- dexp(ti, rate = lambda)
  F.X <- pexp(xi, rate=mu)
  F.T <- pexp(ti, rate = lambda)
  f.XT <- dCopula(copula=rot.cop, cbind(F.X, F.T))
#  f.XT <- f.X*f.T*((1+theta)*((1-F.X)**(-theta)+F.T**(-theta)-1)**(-1/theta-2)*(F.T**(-theta-1))*((1-F.X)**(-theta-1))) #second stage ~ estimate copula 
  logll <- sum(log(f.XT))
  return(-logll)
}

#    normal-expo
haz.mod2.ll <- function(par,xi,ti){
  n <- length(xi)
  beta <- par[1]
  mu <- par[2]
  sigma <- par[3]
  logll <- n*log(beta)+sum(xi)-beta*sum(xi*ti)-(n/2)*log(sigma**2)-(1/(2*sigma**2))*sum((xi-mu)**2)
  return(-logll)
}

cop.mod2.ll <- function(pars,xi,ti){
  theta   <- pars[1]
  mu    <- pars[2]
  sigma <- pars[3]
  lambda <- exp(-pars[4])
  rot.cop <- rotCopula(claytonCopula(theta),flip=c(F,T)) #rotated 90 degrees
  f.X <- dnorm(xi, mean = mu, sd=sigma) #first stage ~ estimate marginal
  f.T <- dexp(ti, rate = lambda)
  F.X <- pnorm(xi, mean = mu, sd=sigma)
  F.T <- pexp(ti, rate = lambda)
  f.XT <- f.X*f.T*dCopula(copula=rot.cop, cbind(F.X,F.T)) #second stage ~ estimate copula 
  # f.XT <- f.X*f.T*((1+theta)*((1-F.X)**(-theta)+F.T**(-theta)-1)**(-1/theta-2)*(F.T**(-theta-1))*((1-F.X)**(-theta-1))) #second stage ~ estimate copula 
  logll <- sum(log(f.XT))
  return(-logll)
}

#    normal-weibull
haz.mod3.ll <- function(par,xi,ti){
  n <- length(xi)
  beta <- par[1]
  mu <- par[2]
  sigma <- par[3]
  gamma <- exp(par[4])
  logll <- n*log(gamma)+(gamma-1)*sum(log(xi*ti))+n*gamma*log(beta)-beta**gamma*sum((xi*ti)**gamma)-(n/2)*log(sigma**2)-(1/(2*sigma**2))*sum((xi-mu)**2)
  return(-logll)
}

cop.mod3.ll <- function(pars,xi,ti){
  theta   <- pars[1]
  mu    <- pars[2]
  sigma <- pars[3]
  shape <- exp(pars[4])
  scale <- exp(pars[5])
  rot.cop <- rotCopula(claytonCopula(theta),flip=c(F,T)) #rotated 90 degrees
  f.X <- dnorm(xi, mean = mu, sd=sigma) #first stage ~ estimate marginal
  f.T <- dweibull(ti, shape = shape, scale = scale)
  F.X <- pnorm(xi, mean = mu, sd=sigma)
  F.T <- pweibull(ti, shape = shape, scale = scale)
  f.XT <- f.X*f.T*dCopula(copula=rot.cop,cbind(F.X,F.T)) #second stage ~ estimate copula 
  # f.XT <- f.X*f.T*((1+theta)*((1-F.X)**(-theta)+F.T**(-theta)-1)**(-1/theta-2)*(F.T**(-theta-1))*((1-F.X)**(-theta-1))) #second stage ~ estimate copula 
  logll <- sum(log(f.XT))
  return(-logll)
}

#    normal-normal
haz.mod4.ll <- function(par,xi,ti){
  n <- length(xi)
  mu.X <- par[1]
  sigma.X <- par[2]
  beta <- par[3]
  sigma.T <- par[4]
  logll <- -n*(log(sigma.T)+log(sigma.X))-.5*(sum((ti-xi*beta)**2))/(sigma.T**2)-.5*(sum((xi-xi*mu.X)**2))/(sigma.X**2)
  return(-logll)
}

cop.mod4.ll <- function(pars,xi,ti){
  theta   <- exp(pars[1])
  mu.X    <- pars[2]
  sigma.X <- pars[3]
  mu.T <- pars[4]
  sigma.T <- pars[5]
  rot.cop <- rotCopula(claytonCopula(theta),flip=c(F,T)) #rotated 90 degrees
  f.X <- dnorm(xi, mean = mu.X, sd=sigma.X) #first stage ~ estimate marginal
  f.T <- dnorm(ti, mean = mu.T, sd=sigma.T)
  F.X <- pnorm(xi, mean = mu.X, sd=sigma.X)
  F.T <- pnorm(ti, mean = mu.T, sd=sigma.T)
  f.XT <- f.X*f.T*dCopula(copula=rot.cop,cbind(F.X,F.T)) #second stage ~ estimate copula 
  # f.XT <- f.X*f.T*((1+theta)*((1-F.X)**(-theta)+F.T**(-theta)-1)**(-1/theta-2)*(F.T**(-theta-1))*((1-F.X)**(-theta-1))) #second stage ~ estimate copula 
  logll <- sum(log(f.XT))
  return(-logll)
}
# =================================================
#  Function to run bootstrapping for each scenario
# =================================================

# expo-expo
sim.mod1 <- function(nrun,sample_size,theta,lambda1,lambda2){
  clay.cop <- rotCopula(archmCopula(family='clayton',dim=2, param=theta), flip=c(F,T))
  mymvd <- mvdc(copula=clay.cop, margins=c('exp','exp'), 
                paramMargins=list(rate=lambda1,rate=lambda2))
  set.seed(123456)
  sim.data <- as.data.frame(rMvdc(500, mymvd))
  mod.par.df <- data.frame()
  for(run in 1:nrun){
    if(run%%10 == 0) cat("Runs =",run,"\n")
    boot.data <- sim.data[sample(nrow(sim.data), sample_size, replace=T), ]
    
    # Fitting for expo-expo
    #     hazard function
    haz.mod.fit.par <- optim(c((sample_size/sum(boot.data[,1]*boot.data[,2])),(sample_size/sum(boot.data[,1]))),haz.mod1.ll,
                              xi=boot.data[,1], ti=boot.data[,2])
    haz.mod.par <- haz.mod.fit.par$par
    
    #     copula function
    # cop.mod.fit.par <- optim(c(theta,(sample_size/sum(boot.data[,1])),(sample_size/sum(boot.data[,2]))), cop.mod1.ll,
    #                           xi=boot.data[,1], ti=boot.data[,2]) 
    cop.mod.fit.par <- optim(c(theta,-log(lambda1),-log(lambda2)), cop.mod1.ll,
                             xi=boot.data[,1], ti=boot.data[,2]) 
    #    storing all parameters for each run
    mod.par.df <- rbind(mod.par.df,c(haz.mod.par[1],haz.mod.par[2],
                                       cop.mod.fit.par$par[1],exp(-cop.mod.fit.par$par[2]),
                                       exp(-cop.mod.fit.par$par[3])))
  }
  
  colnames(mod.par.df) <- c('haz.beta','haz.lambda','cop.theta','marg.lambda1','marg.lambda2')
  # calculating bootstrap S.E and CI
  mod <- data.frame('Estimated'=apply(mod.par.df,2,mean),
                    'CI.lower'=apply(mod.par.df,2,quantile, probs=0.05),
                    'CI.upper'=apply(mod.par.df,2,quantile, probs=0.95),
                    'S.Err'=sqrt(apply(mod.par.df,2,var)))
  write.csv(mod, "C:\\Users\\farea\\Documents\\RWorkspace\\expo_expo_sim.csv")
  
  return(mod)
}

#normal-expo
sim.mod2 <- function(nrun,sample_size,theta,mu,sigma,lambda){
  clay.cop <- rotCopula(archmCopula(family='clayton',dim=2, param=theta), flip=c(F,T))
  mymvd <- mvdc(copula=clay.cop, margins=c('norm','exp'), 
                paramMargins=list(list(mean=mu,sd=sigma),rate=lambda))
  set.seed(10)
  sim.data <- as.data.frame(rMvdc(500, mymvd))
  mod.par.df <- data.frame()
  for(run in 1:nrun){
    if(run%%10 == 0) cat("Runs =",run,"\n")
    boot.data <- sim.data[sample(nrow(sim.data), sample_size, replace=T), ]
    
    # Fitting for normal-expo
    #     hazard function
    haz.mod.fit.par <- optim(c((sample_size/sum(boot.data[,1]*boot.data[,2])),mean(boot.data[,1]),var(boot.data[,1])),haz.mod2.ll,
                             xi=boot.data[,1], ti=boot.data[,2])
    haz.mod.par <- haz.mod.fit.par$par
    
    #     copula function
    cop.mod.fit.par <- optim(c(theta, mu,sigma,-log(lambda)), cop.mod2.ll,
                             xi=boot.data[,1], ti=boot.data[,2]) 
    
    #    storing all parameters for each run
    mod.par.df <- rbind(mod.par.df,c(haz.mod.par[1],haz.mod.par[2],haz.mod.par[3],
                                     cop.mod.fit.par$par[1],cop.mod.fit.par$par[2],
                                     cop.mod.fit.par$par[3],exp(-cop.mod.fit.par$par[4])))
  }
  
  colnames(mod.par.df) <- c('haz.beta','haz.mu','haz.sigma','cop.theta','marg.mu','marg.sigma','marg.lambda')
  # calculating bootstrap S.E and CI
  mod <- data.frame('Estimated'=apply(mod.par.df,2,mean),
                    'CI.lower'=apply(mod.par.df,2,quantile, probs=0.05),
                    'CI.upper'=apply(mod.par.df,2,quantile, probs=0.95),
                    'S.Err'=sqrt(apply(mod.par.df,2,var)))
  write.csv(mod, "C:\\Users\\farea\\Documents\\RWorkspace\\normal_expo_sim.csv")
  
  return(mod)
}

#normal-weibull
sim.mod3 <- function(nrun,sample_size,theta,mu,sigma,shape,scale){
  clay.cop <- rotCopula(archmCopula(family='clayton',dim=2, param=theta), flip=c(F,T))
  mymvd <- mvdc(copula=clay.cop, margins=c('norm','weibull'), 
                paramMargins=list(list(mean=mu,sd=sigma),list(shape=shape,scale=scale)))
  set.seed(7)
  sim.data <- as.data.frame(rMvdc(500, mymvd))
  mod.par.df <- data.frame()
  for(run in 1:nrun){
    if(run%%10 == 0) cat("Runs =",run,"\n")
    boot.data <- sim.data[sample(nrow(sim.data), sample_size, replace=T), ]
    
    # Fitting for normal-weibull
    #     hazard function
    haz.mod.fit.par <- optim(c(.1,mu,sigma,shape),haz.mod3.ll,
                             xi=boot.data[,1], ti=boot.data[,2])
    haz.mod.par <- haz.mod.fit.par$par
    
    #     copula function
    cop.mod.fit.par <- optim(c(log(theta),mu,sigma,log(shape),log(scale)), cop.mod3.ll,
                             xi=boot.data[,1], ti=boot.data[,2]) 
    
    #    storing all parameters for each run
    mod.par.df <- rbind(mod.par.df,c(haz.mod.par[1],haz.mod.par[2],haz.mod.par[3],haz.mod.par[4],
                                     cop.mod.fit.par$par[1],cop.mod.fit.par$par[2],cop.mod.fit.par$par[3],
                                     exp(cop.mod.fit.par$par[4]),exp(cop.mod.fit.par$par[5])))
  }
  
  colnames(mod.par.df) <- c('haz.beta','haz.mu','haz.sigma','haz.shape','cop.theta','marg.mu','marg.sigma','marg.shape','marg.scale')
  # calculating bootstrap S.E and CI
  mod <- data.frame('Estimated'=apply(mod.par.df,2,mean),
                    'CI.lower'=apply(mod.par.df,2,quantile, probs=0.05),
                    'CI.upper'=apply(mod.par.df,2,quantile, probs=0.95),
                    'S.Err'=sqrt(apply(mod.par.df,2,var)))
  write.csv(mod, "C:\\Users\\farea\\Documents\\RWorkspace\\normal_weibull_sim.csv")
  
  return(mod)
}

#normal-normal
sim.mod4 <- function(nrun,sample_size,theta,mu.X,mu.T,sigma.X,sigma.T){
  clay.cop <- rotCopula(archmCopula(family='clayton',dim=2, param=theta), flip=c(F,T))
  mymvd <- mvdc(copula=clay.cop, margins=c('norm','norm'), 
                paramMargins=list(list(mean=mu.X,sd=sigma.X),list(mean=mu.T, sd=sigma.T)))
  set.seed(7)
  sim.data <- as.data.frame(rMvdc(500, mymvd))
  mod.par.df <- data.frame()
  for(run in 1:nrun){
    if(run%%10 == 0) cat("Runs =",run,"\n")
    boot.data <- sim.data[sample(nrow(sim.data), sample_size, replace=T), ]
    
    # Fitting for normal-normal
    #     hazard function
    haz.mod.fit.par <- optim(c(mu.X,sigma.X,.1,sigma.T),haz.mod4.ll,
                             xi=boot.data[,1], ti=boot.data[,2])
    haz.mod.par <- haz.mod.fit.par$par
    
    #     copula function
    cop.mod.fit.par <- optim(c(log(theta),mu.X,sigma.X,mu.T,sigma.T), cop.mod4.ll,
                             xi=boot.data[,1], ti=boot.data[,2]) 
    
    #    storing all parameters for each run
    mod.par.df <- rbind(mod.par.df,c(haz.mod.par[1],haz.mod.par[2],haz.mod.par[3],haz.mod.par[4],
                                     exp(cop.mod.fit.par$par[1]),cop.mod.fit.par$par[2],cop.mod.fit.par$par[3],
                                     cop.mod.fit.par$par[4],cop.mod.fit.par$par[5]))
  }
  
  colnames(mod.par.df) <- c('haz.muX','haz.sigmaX','haz.beta','haz.sigmaT','cop.theta','marg.muX','marg.sigmaX','marg.muT','marg.sigmaT')
  # calculating bootstrap S.E and CI
  mod <- data.frame('Estimated'=apply(mod.par.df,2,mean),
                    'CI.lower'=apply(mod.par.df,2,quantile, probs=0.05),
                    'CI.upper'=apply(mod.par.df,2,quantile, probs=0.95),
                    'S.Err'=sqrt(apply(mod.par.df,2,var)))
  write.csv(mod, "C:\\Users\\farea\\Documents\\RWorkspace\\normal_normal_sim.csv")
  
  return(mod)
}

# =======================================
#  Function for simulation
# =======================================

sim.run <- function(nrun, sample_size=100, sim,...){
  
  stime <- proc.time()
  
  theta <- list(...)$theta
  lambda <- list(...)$lambda
  lambda1 <- list(...)$lambda1
  lambda2 <- list(...)$lambda2
  mu <- list(...)$mu
  sigma <- list(...)$sigma
  shape <- list(...)$shape
  scale <- list(...)$scale
  mu.X <- list(...)$mu.X
  mu.T <- list(...)$mu.T
  sigma.X <- list(...)$sigma.X
  sigma.T <- list(...)$sigma.T
  
  if(sim == 1) result <- sim.mod1(nrun=nrun,sample_size=sample_size,theta=theta,lambda1=lambda1,lambda2=lambda2)
  else if(sim == 2) result <- sim.mod2(nrun=nrun,sample_size=sample_size,theta=theta,mu=mu,sigma=sigma,lambda=lambda)
  else if(sim == 3) result <- sim.mod3(nrun=nrun,sample_size=sample_size,theta=theta,mu=mu,sigma=sigma,shape=shape,scale=scale)
  else if(sim == 4) result <- sim.mod4(nrun=nrun,sample_size=sample_size,theta=theta,mu.X=mu.X,sigma.X=sigma.X,mu.T=mu.T,sigma.T=sigma.T)

  print(proc.time() - stime)
  return(result)
}

sim.run(nrun=100,sample_size=100, sim=1,theta=3.5, lambda1=.3, lambda2=.2) # expo-expo
sim.run(nrun=100,sample_size=100, sim=2,theta=3.5,mu=5,sigma=2,lambda=.2) # normal-expo
sim.run(nrun=100,sample_size=100, sim=3,theta=3.5,mu=5,sigma=2,shape=2,scale=7) # normal-weibull
sim.run(nrun=100,sample_size=100, sim=4,theta=3.5,mu.X=5,sigma.X=2,mu.T=7,sigma.T=2) # normal-normal

# clay.cop <- rotCopula(archmCopula(family='clayton',dim=2, param=4.693), flip=c(F,T))
# mymvd <- mvdc(copula=clay.cop, margins=c('norm','weibull'), 
#               paramMargins=list(list(mean=6.82,sd=7.66),list(shape=0.512,scale=2.613)))
# X <- rMvdc(200, mymvd) # generate X
# contour(mymvd,dMvdc, xlim=c(-1,25),ylim=c(-3,10))
# plot(X)
# pairs.panels(X)

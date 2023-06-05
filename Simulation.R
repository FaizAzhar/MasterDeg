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

#    expo-expo
haz.mod1.ll <- function(xi,ti,par){
  n <- length(xi)
  beta <- par[1]
  lambda <- par[2]
  logll <- n*log(beta) - beta*(sum(xi*ti)) + n*log(lambda) - lambda*sum(xi)
  return(-logll)
}

cop.mod1.ll <- function(pars,xi,ti){
  theta   <- exp(pars[1])
  mu    <- pars[2]
  lambda <- pars[3]
  f.X <- dexp(xi, rate=mu) #first stage ~ estimate marginal
  f.T <- dexp(ti, rate = lambda)
  F.X <- pexp(xi, rate=mu)
  F.T <- pexp(ti, rate = lambda)
  f.XT <- f.X*f.T*((1+theta)*((1-F.X)**(-theta)+F.T**(-theta)-1)**(-1/theta-2)*(F.T**(-theta-1))*((1-F.X)**(-theta-1))) #second stage ~ estimate copula 
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
  theta   <- exp(pars[1])
  mu    <- pars[2]
  sigma <- pars[3]
  lambda <- pars[4]
  f.X <- dnorm(xi, mean = mu, sd=sigma) #first stage ~ estimate marginal
  f.T <- dexp(ti, rate = lambda)
  F.X <- pnorm(xi, mean = mu, sd=sigma)
  F.T <- pexp(ti, rate = lambda)
  f.XT <- f.X*f.T*((1+theta)*((1-F.X)**(-theta)+F.T**(-theta)-1)**(-1/theta-2)*(F.T**(-theta-1))*((1-F.X)**(-theta-1))) #second stage ~ estimate copula 
  logll <- sum(log(f.XT))
  return(-logll)
}

#    normal-weibull
haz.mod3.ll <- function(par,xi,ti){
  n <- length(xi)
  beta <- par[1]
  gamma <- exp(par[2])
  mu <- par[3]
  sigma <- par[4]
  logll <- n*log(gamma)+(gamma-1)*sum(log(xi*ti))+n*gamma*log(beta)-beta**gamma*sum((xi*ti)**gamma)-(n/2)*log(sigma**2)-(1/(2*sigma**2))*sum((xi-mu)**2)
  return(-logll)
}

cop.mod3.ll <- function(pars,xi,ti){
  theta   <- exp(pars[1])
  mu    <- pars[2]
  sigma <- pars[3]
  shape <- exp(pars[4])
  scale <- exp(pars[5])
  f.X <- dnorm(xi, mean = mu, sd=sigma) #first stage ~ estimate marginal
  f.T <- dweibull(ti, shape = shape, scale = scale)
  F.X <- pnorm(xi, mean = mu, sd=sigma)
  F.T <- pweibull(ti, shape = shape, scale = scale)
  f.XT <- f.X*f.T*((1+theta)*((1-F.X)**(-theta)+F.T**(-theta)-1)**(-1/theta-2)*(F.T**(-theta-1))*((1-F.X)**(-theta-1))) #second stage ~ estimate copula 
  logll <- sum(log(f.XT))
  return(-logll)
}

# =======================================
#  Function for simulation
# =======================================

sim.run <- function(nrun, sample_size=100, mu, sigma){
  
  stime <- proc.time()
  set.seed(24)
  sim.data <- as.data.frame(mvrnorm(n=sample_size, mu=mu, Sigma=sigma))
  mod1.par.df <- data.frame()
  mod2.par.df <- data.frame()
  mod3.par.df <- data.frame()
  
  for(run in 1:nrun){
    boot.data <- sim.data[sample(nrow(sim.data), sample_size, replace=T), ]

    # Fitting for expo-expo
    #     hazard function
    haz.mod1.fit.par <- optim(c(0.1,0.1),haz.mod1.ll,
                              xi=boot.data[,1], ti=boot.data[,2],
                              hessian=T)
    haz.mod1.par <- haz.mod1.fit.par$par

    #     copula function
    cop.mod1.fit.par <- optim(c(1, 1, 1), cop.mod1.ll,
                              xi=boot.data[,1], ti=boot.data[,2],
                              hessian=T) 

    #    storing all parameters for each run
    mod1.par.df <- rbind(mod1.par.df,c(haz.mod1.par[1],haz.mod1.par[2],
                                       exp(cop.mod1.fit.par$par[1]),cop.mod1.fit.par$par[2],
                                       cop.mod1.fit.par$par[3]))

    # Fitting for normal-expo
    #     hazard function
    haz.mod2.fit.par <- optim(rep(0.1,3),haz.mod2.ll,
                              xi=boot.data[,1], ti=boot.data[,2],
                              hessian=T)
    haz.mod2.par <- haz.mod2.fit.par$par

    #     copula function
    cop.mod2.fit.par <- optim(c(0, 6, 6,0.1), cop.mod2.ll,
                              xi=boot.data[,1], ti=boot.data[,2],
                              hessian=T) 
    
    #    storing all parameters for each run
    mod2.par.df <- rbind(mod2.par.df,c(haz.mod2.par[1],haz.mod2.par[2],haz.mod2.par[3],
                                       exp(cop.mod2.fit.par$par[1]),cop.mod2.fit.par$par[2],
                                       cop.mod2.fit.par$par[3],cop.mod2.fit.par$par[4]))

    # Fitting for normal-weibull
    #     hazard function
    haz.mod3.fit.par <- optim(rep(0.1,4),haz.mod3.ll,
                              xi=boot.data[,1], ti=boot.data[,2],
                              hessian=T)
    haz.mod3.par <- haz.mod3.fit.par$par
    
    #     copula function
    cop.mod3.fit.par <- optim(c(0, 6, 6, 0, 8), cop.mod3.ll,
                              xi=boot.data[,1], ti=boot.data[,2])   

    #    storing all parameters for each run
    mod3.par.df <- rbind(mod3.par.df,c(haz.mod3.par[1],exp(haz.mod3.par[2]),haz.mod3.par[3],haz.mod3.par[4],
                                       exp(cop.mod3.fit.par$par[1]),cop.mod3.fit.par$par[2],
                                       cop.mod3.fit.par$par[3],exp(cop.mod3.fit.par$par[4]),
                                       exp(cop.mod3.fit.par$par[5])))
  }
  
  colnames(mod1.par.df) <- c('haz.beta','haz.lambda','cop.rho','marg.lambda1','marg.lambda2')
  colnames(mod2.par.df) <- c('haz.beta','haz.mu','haz.sigma','cop.rho','marg.mu','marg.sigma','marg.lambda')
  colnames(mod3.par.df) <- c('haz.beta','haz.gamma','haz.mu','haz.sigma','cop.rho','marg.mu','marg.sigma','marg.lambda','marg.gamma')

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
  print(proc.time() - stime)
  
  return(list('expo-expo'=mod1,'normal-expo'=mod2,'normal-weibull'=mod3))
}

sim.run(100,sample_size=100,mu=c(5,7),sigma=matrix(c(1,-.9,-.9,1),ncol=2))

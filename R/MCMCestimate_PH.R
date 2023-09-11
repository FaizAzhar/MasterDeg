MCMCestimate_PH <- function(model, df, f.path, p.jags, i.jags, n.chains=3,
                             n.iter=50000, n.thin=1, n.burnin=2000){
  n <- nrow(df)
  t <- df$t # time column
  x <- df$x # biomarker column
  d <- df$status # censor or not (0 or 1)
  d.jags <- list(n=n,Time=t,X=x,delta=d,zeros=rep(0,n)) # data for JAGS

  if(model == 'exp_exp'){
    mod = function(){
      for(i in 1:n){X[i]~dexp(mu.X)}
      for(i in 1:n){
        ph.const[i] <- exp(beta*X[i])
        logHaz[i] <- log(lambda.T * ph.const[i])
        logSurv[i] <- log(pow(exp(-lambda.T*Time[i]),ph.const[i]))
        logfX[i] <- log(mu.X*exp(-mu.X*X[i]))

        # using zeros-trick instead of standard log-likelihood
        L[i] <- delta[i] * (logHaz[i] + logfX[i]) + logSurv[i]
        phi[i] <- 10^5 - L[i]
        zeros[i]~dpois(phi[i])
      }
      # prior distribution
      Deviance <- -2*sum(L[1:n])
      beta~dnorm(0,0.001)
      lambda.T~dgamma(0.001,0.001)
      mu.X~dgamma(0.001,0.001)
    }
    model.file = file.path(f.path,"model.txt")
    write.model(mod,model.file)
  }
  else if(model == 'norm_exp'){
    mod = function(){
      for(i in 1:n){X[i]~dnorm(mu.X,tau)}
      for(i in 1:n){
        ph.const[i] <- exp(beta*X[i])
        logHaz[i] <- log(lambda.T) + log(ph.const[i])
        logSurv[i] <- -lambda.T*Time[i]*ph.const[i]
        logfX[i] <- -pow(X[i]-mu.X,2)/(2*pow(sigma.X,2)) - log(sigma.X*pow(2*(22/7),1/2))

        # using zeros-trick instead of standard log-likelihood
        L[i] <- delta[i] * (logHaz[i] + logfX[i]) + logSurv[i]
        phi[i] <- 10^5 - L[i]
        zeros[i]~dpois(phi[i])
      }
      # prior distribution
      Deviance <- -2*sum(L[1:n])
      beta~dnorm(0,0.001)
      lambda.T~dgamma(0.001,0.001)
      mu.X~dnorm(0,10)
      tau ~ dnorm(0,0.000001)
      sigma.X <- 1/sqrt(tau)
    }
    model.file = file.path(f.path,"model.txt")
    write.model(mod,model.file)
  }
  else if(model == 'norm_weib'){
    mod = function(){
      for(i in 1:n){X[i]~dnorm(mu.X,tau)}
      for(i in 1:n){
        ph.const[i] <- exp(beta*X[i])
        logHaz[i] <- log(shape.T) + log(scale.T) + (shape.T-1)*log(scale.T*Time[i]) + log(ph.const[i])
        logSurv[i] <- -pow(scale.T*Time[i],shape.T)*ph.const[i]
        logfX[i] <- -pow(X[i]-mu.X,2)/(2*pow(sigma.X,2)) - log(sigma.X*pow(2*(22/7),1/2))

        # using zeros-trick instead of standard log-likelihood
        L[i] <- delta[i] * (logHaz[i] + logfX[i]) + logSurv[i]
        phi[i] <- 10^5 - L[i]
        zeros[i]~dpois(phi[i])
      }
      # prior distribution
      Deviance <- -2*sum(L[1:n])
      beta~dnorm(0,0.001)
      shape.T~dgamma(0.001,0.001)
      scale.T~dgamma(0.001,0.001)
      mu.X~dnorm(0,0.001)
      tau ~ dnorm(0,0.000001)
      sigma.X <- 1/sqrt(tau)
    }
    model.file = file.path(f.path,"model.txt")
    write.model(mod,model.file)
  }

  jmod <- jags(model=model.file, data=d.jags, n.chains=n.chains, inits=i.jags,
               parameters.to.save=p.jags, n.iter=n.iter, n.thin=n.thin, n.burnin=n.burnin)
  print(jmod,intervals=c(0.025,0.975),digits=4)
  jmod.mcmc <- as.mcmc(jmod)
  print(gelman.diag(jmod.mcmc))

  if(model == 'exp_exp'){
    write.table(MCMCsummary(jmod.mcmc,round=4),file=file.path(f.path,"Simulation_result","MCMC_ph_expo_expo_result.txt"))
    MCMCtrace(jmod.mcmc,params=p.jags,ISB=F,exact=T,pdf=T, filename='MCMC_ph_expo_expo_traceplot.pdf',
              open_pdf=F)
  }
  else if(model == 'norm_exp'){
    write.table(MCMCsummary(jmod.mcmc,round=4),file.path(f.path,'Simulation_result','MCMC_ph_norm_expo_result.txt'))
    MCMCtrace(jmod.mcmc,params=p.jags,ISB=F,exact=T,pdf=T, filename='MCMC_ph_norm_expo_traceplot.pdf',
              open_pdf=F)
  }
  else if(model == 'norm_weib'){
    write.table(MCMCsummary(jmod.mcmc,round=4),file.path(f.path,'Simulation_result','MCMC_ph_norm_weib_result.txt'))
    MCMCtrace(jmod.mcmc,params=p.jags,ISB=F,exact=T,pdf=T, filename='MCMC_ph_norm_weib_traceplot.pdf',
              open_pdf=F)
  }
}

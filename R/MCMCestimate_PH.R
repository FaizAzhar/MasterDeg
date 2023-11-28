#' MCMCestimate_PH
#'
#' Function to run estimation process of joint model discussed in this package using MCMC Gibbs Sampling.
#'
#' @param model A string which specify the model of interest.\cr
#'              model = c('exp_exp','norm_exp','norm_weib').
#' @param df A dataframe with 3 columns ('x' = biomarker, 't' = time-to-event, 'status' = censoring status (0 or 1))
#' @param f.path A string that specify file path to store value of estimated parameters.
#' @param p.jags A vector of parameters that need to be monitored.\cr
#'               'exp_exp':  c('lambda.T','beta','mu.X','Deviance')\cr
#'               'norm_exp': c('lambda.T','beta','mu.X','sigma.X','Deviance')\cr
#'               'norm_weib': c('scale.T','shape.T','beta','mu.X','sigma.X','Deviance')\cr
#'               Deviance: Average of '-2*loglikelihood' (Dbar) - Will be used for manual DIC calculation.
#' @param i.jags A function that generate the initial value of parameter used in JAGS.\cr
#'               'exp_exp': function(){list(lambda.T=lambda.t,beta=beta, mu.X=mu.x)}\cr
#'               'norm_exp': function(){list(lambda.T=lambda.t,beta=beta, mu.X=mu.x, tau=runif(1))}\cr
#'               'norm_weib': function(){list(shape.T=shape.t,scale.T=scale.t,beta=beta, mu.X=mu.x, tau=runif(1))}\cr
#'               beta: A value of association measurement\cr
#'               lambda.T: A scale parameter of T\cr
#'               scale.T: A scale parameter of T (for weibull)\cr
#'               shape.T: A shape parameter of T (for weibull)\cr
#'               mu.x: A scale parameter for X\cr
#'               tau: A standard deviation for X (JAGS suggest using tau)
#' @param n.chains A number of chains to run MCMC process.
#' @param n.iter A number of iteration run for each chain.
#' @param n.thin A number of iteration thinning.
#' @param n.burnin A number that specify how many iteration will be discarded.
#' @export
#' @returns Print MCMC summary & automatically save the summary to the specified file path.
#' @examples
#' mu.x=5; sigma.x=0.8; shape.t=1.2; scale.t=1/1.4; tau=0.7; theta=2*tau/(1-tau)
#' f.path='C:/Users/farea/Documents/RWorkspace'
#' p.jags <- c('scale.T','shape.T','theta','mu.X','sigma.X','Deviance')
#' i.jags <- function(){list(scale.T=scale.t,shape.T=shape.t,theta=theta, mu.X=mu.x, tau=runif(1))}
#' MCMCestimate_COP(model='norm_weib', df=df.COP, f.path=f.path, p.jags=p.jags, i.jags=i.jags)



MCMCestimate_PH <- function(model, df, f.path, p.jags=NULL, i.jags=NULL, n.chains=3,
                             n.iter=20000, n.thin=5, n.burnin=2000){
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
  else if(model == 'snorm_snorm'){
    mod.stan =
      "
      data {
        int<lower=0> n; // no. of obs
        real X[n]; // marker value
        real Time[n]; // time-to-event
        real<lower=0> delta[n]; // status (0=censor, 1=not censor)
      }

      parameters {
        real muX;
        real<lower=0> sigmaX;
        real skewX;
        real beta;
        real muT;
        real<lower=0> sigmaT;
        real skewT;
      }

      model {
        beta~normal(0,3);
        muX~normal(0,3);
        skewX ~ normal(0,3);
        sigmaX ~ gamma(0.5,0.5);
        muT~normal(0,3);
        skewT ~ normal(0,3);
        sigmaT ~ gamma(0.5,0.5);
        //X ~ skew_normal(muX, sigmaX, skewX);
        //Time ~ skew_normal(muT, sigmaT, skewT);
        for(i in 1:n)
         target += delta[i]*((beta*X[i]) + (exp(beta*X[i])-1)*skew_normal_lccdf(Time[i]| muT, sigmaT, skewT) + skew_normal_lpdf(Time[i] | muT, sigmaT, skewT) + skew_normal_lpdf(X[i] | muX, sigmaX, skewX)) + (1-delta[i])*(exp(beta*X[i])*skew_normal_lccdf(Time[i]| muT, sigmaT, skewT));
      }
      "
    model.file = file.path(f.path, "model.stan")
    write(mod.stan, model.file)
  }

  if(model != "snorm_snorm"){
   jmod <- jags(model=model.file, data=d.jags, n.chains=n.chains, inits=i.jags,
               parameters.to.save=p.jags, n.iter=n.iter, n.thin=n.thin, n.burnin=n.burnin)
  print(jmod,intervals=c(0.025,0.975),digits=4)
  jmod.mcmc <- as.mcmc(jmod)
  print(gelman.diag(jmod.mcmc))
  }else{
    fit <- stan(file = model.file, data = d.jags, chains=n.chains, iter=n.iter, warmup=n.burnin)
    saveRDS(fit, "snorm_snorm.rds")
    print(fit)
    }


  if(model == 'exp_exp'){
    write.table(MCMCsummary(jmod.mcmc,round=4),file=file.path(f.path,"MCMC_ph_expo_expo_result.txt"))
    MCMCtrace(jmod.mcmc,params=p.jags,ISB=F,exact=T,pdf=T, filename='MCMC_ph_expo_expo_traceplot.pdf',
              open_pdf=F)
  }
  else if(model == 'norm_exp'){
    write.table(MCMCsummary(jmod.mcmc,round=4),file.path(f.path,'MCMC_ph_norm_expo_result.txt'))
    MCMCtrace(jmod.mcmc,params=p.jags,ISB=F,exact=T,pdf=T, filename='MCMC_ph_norm_expo_traceplot.pdf',
              open_pdf=F)
  }
  else if(model == 'norm_weib'){
    write.table(MCMCsummary(jmod.mcmc,round=4),file.path(f.path,'MCMC_ph_norm_weib_result.txt'))
    MCMCtrace(jmod.mcmc,params=p.jags,ISB=F,exact=T,pdf=T, filename='MCMC_ph_norm_weib_traceplot.pdf',
              open_pdf=F)
  }
}

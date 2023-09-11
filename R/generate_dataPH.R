generate_dataPH <- function(mod,params,n.sample){
  mu.x <- params$mu.x; lambda.t <- params$lambda.t; beta <- params$beta; c.rate <- params$c.rate;
  sigma.x <- params$sigma.x; shape.t <- params$shape.t; scale.t <- params$scale.t
  c <- rep(NA,n.sample)

  if(mod == 'exp_exp'){
    x <- rexp(n.sample,mu.x)
    t.true <- -log(runif(n.sample))/(lambda.t*exp(x*beta))

    #censoring data
    for(i in 1:n.sample){c[i] <- runif(1,0,1/(exp(x[i]*beta)*lambda.t*c.rate))}

    # setting observable time-to-event
    t <- ifelse(t.true < c, t.true,c)
    status <- ifelse(t.true < c, 1,0)
    simdata <- data.frame(x=x,t=t,status=status)
  }
  else if(mod == 'norm_exp'){
    x <- rnorm(n.sample,mu.x,sigma.x)
    t.true <- -log(runif(n.sample))/(lambda.t*exp(x*beta))

    #censoring data
    for(i in 1:n.sample){c[i] <- runif(1,0,1/(exp(x[i]*beta)*lambda.t*c.rate))}

    # setting observable time-to-event
    t <- ifelse(t.true < c, t.true,c)
    status <- ifelse(t.true < c, 1,0)
    simdata <- data.frame(x=x,t=t,status=status)
  }
  else if(mod == 'norm_weib'){
    x <- rnorm(n.sample,mu.x,sigma.x)
    t.true <- (-log(runif(n.sample))/(scale.t*exp(x*beta)))^(1/shape.t)

    #censoring data
    for(i in 1:n.sample){c[i] <- runif(1,0,gamma(1+1/shape.t)/((exp(x[i]*beta)^(1/shape.t))*scale.t*c.rate))}

    # setting observable time-to-event
    t <- ifelse(t.true < c, t.true,c)
    status <- ifelse(t.true < c, 1,0)
    simdata <- data.frame(x=x,t=t,status=status)
  }
  return(simdata)
}

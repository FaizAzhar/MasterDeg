generate_dataCOP <- function(mod, params, n.sample){
  mu.x <- params$mu.x; lambda.t <- params$lambda.t; theta <- params$theta; c.rate <- params$c.rate;
  sigma.x <- params$sigma.x;
  shape.t <- params$shape.t; scale.t <- params$scale.t;

  if(mod == 'exp_exp'){
    clay.cop <- rotCopula(archmCopula(family='clayton',dim=2,param=theta), flip=c(T,F)) # rotate 90 clayton
    clay.dist <- mvdc(copula=clay.cop, margins=c('exp','exp'),
                      paramMargins=list(rate=mu.x,rate=lambda.t))

    simdata <- as.data.frame(rMvdc(n.sample, clay.dist))
    colnames(simdata) <- c('x','t')

    # censoring data
    c <- runif(n.sample,0,1/(c.rate*lambda.t))

    # setting observable time-to-event
    t.true <- simdata$t
    t <- ifelse(t.true < c, t.true,c)
    status <- ifelse(t.true < c, 1,0)

    simdata <- data.frame(x=simdata$x,t=t,status=status)
  }
  else if(mod == 'norm_exp'){
    clay.cop <- rotCopula(archmCopula(family='clayton',dim=2,param=theta), flip=c(T,F)) # rotate 90 clayton
    clay.dist <- mvdc(copula=clay.cop, margins=c('norm','exp'),
                      paramMargins=list(list(mean=mu.x,sd=sigma.x),rate=lambda.t))

    simdata <- as.data.frame(rMvdc(n.sample, clay.dist))
    colnames(simdata) <- c('x','t')

    # censoring data
    c <- runif(n.sample,0,1/(c.rate*lambda.t))

    # setting observable time-to-event
    t.true <- simdata$t
    t <- ifelse(t.true < c, t.true,c)
    status <- ifelse(t.true < c, 1,0)

    simdata <- data.frame(x=simdata$x,t=t,status=status)
  }
  else if(mod == 'norm_weib'){
    clay.cop <- rotCopula(archmCopula(family='clayton',dim=2,param=theta), flip=c(T,F)) # rotate 90 clayton
    clay.dist <- mvdc(copula=clay.cop, margins=c('norm','weibull'),
                      paramMargins=list(list(mean=mu.x,sd=sigma.x),list(shape=shape.t,scale=scale.t)))

    simdata <- as.data.frame(rMvdc(n.sample, clay.dist))
    colnames(simdata) <- c('x','t')

    # 0.2 censoring rate
    c <- runif(n.sample,0,gamma(1/shape.t+1)/(c.rate/scale.t))

    # setting observable time-to-event
    t.true <- simdata$t
    t <- ifelse(t.true < c, t.true,c)
    status <- ifelse(t.true < c, 1,0)

    simdata <- data.frame(x=simdata$x,t=t,status=status)
  }
  return(simdata)
}

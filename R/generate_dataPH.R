#' generate_dataPH
#'
#' Function to generate bivariate data using Proportional Hazard Model using inverse sampling technique
#'
#' @param mod A string that specify which model to be implemented to join the bivariate data (join biomarker X and time-to-event T).\cr
#'            mod = c('exp_exp', 'norm_exp', 'norm_weib')
#' @param params A list of parameters for selected model.\cr
#'               'exp_exp': params = list(c.rate, beta, lambda.t, mu.x).\cr
#'               'norm_exp': params = list(c.rate, beta, lambda.t, mu.x, sigma.x).\cr
#'               'norm_weib': params = list(c.rate, beta, scale.t, shape.t, mu.x, sigma.x).\cr
#'               c.rate: censoring rate (0 to 1)\cr
#'               beta: association measurement\cr
#'               lambda.t: scale parameter of T\cr
#'               scale.t: scale parameter of T (for weibull)\cr
#'               shape.t: shape parameter of T (for weibull)\cr
#'               mu.x: scale parameter of X\cr
#'               sigma.x: standard deviation of X
#' @param n.sample A number of sample size to be generates.
#' @returns A dataframe with 3 columns (x = biomarker value, t = observable time-to-event, status = censored/not censor (0 or 1))
#' @examples
#' data <- generate_dataPH('exp_exp', params=list(c.rate=0.9, beta=1.2, lambda.t=0.5, mu.x=1), n.sample=200)
#' data <- generate_dataPH('norm_exp', params=list(c.rate=0.9, beta=0.8, lambda.t=0.5, mu.x=5, sigma.x=0.8), n.sample=200)
#' data <- generate_dataPH('norm_weib', params=list(c.rate=0.9, beta=0.2, scale.t=1/1.5, shape.t=1.2, mu.x=5, sigma.x=0.8), n.sample=200)

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

#' generate_dataCOP
#'
#' Function to generate bivariate data using Rotated 90 Clayton copula
#'
#' @param mod A string that specify which model to be implemented to join the bivariate data (join biomarker X and time-to-event T).\cr
#'            mod = c('exp_exp', 'norm_exp', 'norm_weib')
#' @param params A list of parameters for selected model.\cr
#'               'exp_exp': params = list(c.rate, theta, lambda.t, mu.x).\cr
#'               'norm_exp': params = list(c.rate, theta, lambda.t, mu.x, sigma.x).\cr
#'               'norm_weib': params = list(c.rate, theta, scale.t, shape.t, mu.x, sigma.x).\cr
#'               c.rate: censoring rate (0 to 1)\cr
#'               theta: dependence measurement (tau (0 to 1);theta=2*tau/(1-tau))\cr
#'               lambda.t: scale parameter of T\cr
#'               scale.t: scale parameter of T (for weibull)\cr
#'               shape.t: shape parameter of T (for weibull)\cr
#'               mu.x: scale parameter of X\cr
#'               sigma.x: standard deviation of X
#' @param n.sample A number of sample size to be generates.
#' @returns A dataframe with 3 columns (x = biomarker value, t = observable time-to-event, status = censored/not censor (0 or 1))
#' @examples
#' data <- generate_dataCOP('exp_exp', params=list(c.rate=0.9, theta=2*(0.8)/(1-0.8), lambda.t=0.5, mu.x=1), n.sample=200)
#' data <- generate_dataCOP('norm_exp', params=list(c.rate=0.9, theta=2*(0.8)/(1-0.8), lambda.t=0.5, mu.x=5, sigma.x=0.8), n.sample=200)
#' data <- generate_dataCOP('norm_weib', params=list(c.rate=0.9, theta=2*(0.8)/(1-0.8), scale.t=1/1.5, shape.t=1.2, mu.x=5, sigma.x=0.8), n.sample=200)

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

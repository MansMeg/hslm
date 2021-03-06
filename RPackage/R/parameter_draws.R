#' Draw beta parameters
#'
#' @details 
#' See the file derivations for details on how to draw the betas.
#'
#' @param mu_n The posterior beta expectation
#' @param Lambda_n_inv The inverse of the posterior prescision matrix.
#' @param sigma The sigma parameter
#'
draw_beta <- function(mu_n, Lambda_n_inv, sigma){
  mvtnorm::rmvnorm(n=1, mean=mu_n, sigma = sigma^2 * Lambda_n_inv)
}



#' Draw sigma
#'
#' @details 
#' See the file derivations for details on how to draw sigma.
#'
#' @param a_n posterior shape
#' @param b_0 prior scale
#' @param yty y^t * y
#' @inheritParams draw_beta
#'
draw_sigma <- function(a_n, b_0, yty, mu_n, Lambda_n){
  b_n <- b_0 + 0.5 * (yty - t(mu_n) %*% Lambda_n %*% mu_n)
  sqrt(1/rgamma(n = 1, shape = a_n, rate = b_n))
}


#' Draw lambda
#'
#' @details 
#' See the file derivations for details on how to draw lambda.
#'
#' @param lambda Previous lambda draw
#' @param beta Beta parameters
#' @param sigma Sigma parameter
#' @param tau Tau parameter
#'
#'
draw_lambda <- function(lambda, beta, sigma, tau){
  gamma_l <- 1 / lambda^2
  u1 <- runif(length(lambda), 0, 1 / (1 + gamma_l))
  trunc_limit <- (1 - u1) / u1
  mu2_j <- (beta / (sigma * tau))^2
  rate_lambda <- (mu2_j / 2)
  ub_lambda <- pexp(trunc_limit, rate_lambda)
  u2 <- runif(length(ub_lambda), 0, ub_lambda)
  gamma_l <- qexp(u2, rate_lambda)
  return(1 / sqrt(gamma_l))
}



#' Draw tau
#'
#' @details 
#' See the file derivations for details on how to draw tau.
#'
#' @inheritParams draw_lambda
#' @param shape_tau Tau shape parameter
#'
#'
draw_tau <- function(lambda, beta, sigma, tau){
  shape_tau <- 0.5 * (length(lambda) + 1)
  gamma_t <- 1 / tau^2
  u1 <- runif(1, 0, 1 / (1 + gamma_t))
  trunc_limit_tau <- (1 - u1) / u1
  mu2_tau <- sum((beta / (sigma * lambda))^2)
  rate_tau <- (mu2_tau / 2)
  ub_tau <- pgamma(trunc_limit_tau, shape=shape_tau, rate=rate_tau)
  u2 <- runif(1, 0, ub_tau)
  gamma_t <- qgamma(p = u2, shape = shape_tau, rate = rate_tau)
  return(1 / sqrt(gamma_t))
}


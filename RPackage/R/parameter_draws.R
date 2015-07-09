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




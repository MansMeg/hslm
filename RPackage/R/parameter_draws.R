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




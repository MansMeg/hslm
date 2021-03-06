#' Bayeisan linear regression using the horseshoe prior
#'
#' @description 
#' Function to estimate a linear regression with the horseshoe prior. 
#' Model specification:
#' \deqn{y ~ N(X\beta^T,\sigma ^2)}
#' \deqn{\beta | \lambda_i, \tau ~ N(0, \tau^2  \lambda_i^2 \sigma^2)}
#' \deqn{\sigma^2 ~ IG(a, b)}
#' \deqn{\lambda_i ~ C_+(0, 1)}
#' \deqn{\tau ~ C_+(0, 1)}
#'
#' @param y dependent variabel as a column matrix. 
#' @param x independent variables in a matrix. 
#' @param intercept Should an intercept be included in the model. Default is FALSE.
#' @param ab prior for sigma in model. If NULL, sigma is assumed to be known and 1. Default is c(1, 1).
#' @param iter The number of MCMC iterations. Default is 2000.
#'
#' @return 
#' List with
#' \describe{
#'   \item{\code{beta}}{Matrix of size \code{iter} * (\code{ncol(x) + intercept})}
#'   \item{\code{lambda}}{Matrix of size \code{iter} * (\code{ncol(x) + intercept})}
#'   \item{\code{tau}}{Vector of length \code{iter}}
#'   \item{\code{sigma}}{Vector of length \code{iter}, if estimated.}
#' }
#' 
#' @examples
#' data(diabetes_x_train)
#' data(diabetes_y_train)
#' hs_res <- hslm(diabetes_y_train, diabetes_x_train)
#'
#' @export
hslm <- function(y, x, iter=2000, intercept=FALSE, ab=c(1,1)){
  # Assertions
  stopifnot(is.matrix(y) | is.data.frame(y),
            is.matrix(x) | is.data.frame(x),
            is.numeric(iter) & length(iter) == 1,
            is.logical(intercept) & length(intercept) == 1,
            is.numeric(ab) & length(ab) == 2 & all(ab > 0))
  
  # Priors
  if(!is.null(ab)){
    stopifnot(length(ab)==2)
    a_0 <- ab[1]
    b_0 <- ab[2]
  }
  
  # Param storage
  beta_bayes_hs <- matrix(as.numeric(NA),ncol=ncol(x), nrow=iter)
  sigma <- rep(as.numeric(NA), iter)
  lambda <- gamma_l <- matrix(as.numeric(NA), ncol=length(x), nrow=iter)
  tau <- gamma_t <- rep(as.numeric(NA), iter) 
  Lambda0 <- diag(lambda[1,])

  # Inits
  beta_bayes_hs[1,] <- 0
  lambda[1,] <- 1
  sigma[1] <- tau[1] <- 1
  
  # Convert y and x to matrices
  y <- as.matrix(y)
  if(intercept){
    ff <- y~.
  } else {
    ff <- y~. -1
  }
  m <- model.frame(ff, x)
  x <- model.matrix(ff, m)
  rm(m)

  # Precalculations
  XtX <- t(x)%*%x
  Xty <- t(x)%*%y
  yty <- t(y)%*%y
  if(!is.null(ab)){
    a_n <- a_0 + length(y)/2
  } else {
    sigma <- rep(1, iter)
  }
  
  # Sampling 
  for(it in 2:iter){
    # Calculate Lambda0, Lambda_n, Lambda_n_inv, mu_n
    diag(Lambda0) <- 1 / (lambda[it - 1,]^2 * tau[it - 1]^2)
    Lambda_n <- XtX + Lambda0
    Lambda_n_inv <- solve(Lambda_n)
    mu_n <- Lambda_n_inv %*% Xty
    
    # Sampling beta
    beta_bayes_hs[it,] <- draw_beta(mu_n, Lambda_n_inv, sigma = sigma[it-1])
    
    # Sampling sigma
    if(!is.null(ab)){
      sigma[it] <- draw_sigma(a_n, b_0, yty, mu_n, Lambda_n)
    }
    
    # Sampling lambda
    lambda[it,] <- draw_lambda(lambda = lambda[it-1,], 
                               beta = beta_bayes_hs[it, ], 
                               sigma = sigma[it], 
                               tau = tau[it-1])
    
    # Sampling tau
    tau[it] <- draw_tau(lambda = lambda[it,], 
                        beta = beta_bayes_hs[it, ], 
                        sigma = sigma[it], 
                        tau = tau[it-1])
  }
  # Return results
  res <- list(beta=beta_bayes_hs,
              lambda=lambda,
              tau=tau)
  if(!is.null(ab)){
    res <- c(res, list(sigma=sigma))
  }
  return(res)
}

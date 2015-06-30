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
  
  # Convert to matrices
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
  shape_tau <- 0.5 * (ncol(lambda) + 1)
  if(!is.null(ab)){
    a_n <- a_0 + length(y)/2
  } else {
    sigma <- rep(1, iter)
  }
  
  # Sampling 
  for(it in 2:iter){
    # Create Lambda0 matrix
    diag(Lambda0) <- 1 / (lambda[it - 1,]^2 * tau[it - 1]^2)
    
    # Sampling beta
    Lambda_n <- XtX + Lambda0
    Lambda_n_inv <- solve(Lambda_n)
    mu_n <- Lambda_n_inv %*% Xty
    beta_bayes_hs[it,] <- mvtnorm::rmvnorm(n=1, mean=mu_n, sigma = sigma[it-1]^2 * Lambda_n_inv)
    
    # Sampling sigma
    if(!is.null(ab)){
      b_n <- b_0 + 0.5 * (yty - t(mu_n) %*% Lambda_n %*% mu_n)
      sigma[it] <- sqrt(MCMCpack::rinvgamma(n=1, a_n, b_n))
    }
    
    # Sampling lambda
    gamma_l[it,] <- 1 / lambda[it-1,]^2
    u1 <- runif(ncol(lambda), 0, 1 / (1 + gamma_l[it,]))
    trunc_limit <- (1 - u1) / u1
    mu2_j <- (beta_bayes_hs[it, ] / (sigma[it] * tau[it-1]))^2
    rate_lambda <- (mu2_j / 2)
    ub_lambda <- pexp(trunc_limit, rate_lambda)
    u2 <- runif(length(ub_lambda), 0, ub_lambda)
    gamma_l[it,] <- qexp(u2, rate_lambda)
    lambda[it,] <- 1 / sqrt(gamma_l[it,])
    
    # Sampling tau
    gamma_t[it] <- 1 / tau[it-1]^2
    u1 <- runif(1, 0, 1 / (1 + gamma_t[it]))
    trunc_limit_tau <- (1 - u1) / u1
    mu2_tau <- sum((beta_bayes_hs[it, ] / (sigma[it] * lambda[it,]))^2)
    rate_tau <- (mu2_tau / 2)
    ub_tau <- pgamma(trunc_limit_tau, shape=shape_tau, rate=rate_tau)
    u2 <- runif(1, 0, ub_tau)
    gamma_t[it] <- qgamma(p = u2, shape = shape_tau, rate = rate_tau)
    tau[it] <- 1 / sqrt(gamma_t[it])
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

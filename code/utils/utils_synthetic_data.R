library(confsurv)

## Define mixture distribution
sample_mixture <- function(generator, generator_shift , num_samples, prop_shift) {
    n_shift <- round(prop_shift * num_samples)
    n_noshift <- num_samples - n_shift
    parts <- list()
    if (n_shift > 0) parts[[length(parts) + 1]] <- generator_shift$sample(n_shift)
    if (n_noshift > 0) parts[[length(parts) + 1]] <- generator$sample(n_noshift)
    data_all <- do.call(rbind, parts)
    data_all <- data_all[sample(nrow(data_all)), , drop = FALSE]
    return(data_all)
}

generate_correlated_features <- function(num_samples, num_features, rho = 0.3) {
  # Create a covariance matrix with constant correlation
  Sigma <- matrix(rho, nrow = num_features, ncol = num_features)
  diag(Sigma) <- 1  # Variance of 1 on the diagonal
  
  # Sample from multivariate normal
  raw_X <- MASS::mvrnorm(n = num_samples, mu = rep(0, num_features), Sigma = Sigma)
  
  # Optionally transform to [0,1] scale (similar to runif)
  X <- pnorm(raw_X)  # Transforms to uniform via CDF of standard normal
  
    return(X)
}

init_data_generator <- function(setting) {
    if(setting==4) {
        ## New setting (more interesting)
        num_features <- 100
        covariate_generator <- function(num_samples) {
            matrix(runif(num_samples * num_features, -1, 1), nrow = num_samples)
        }
        surv_mu_fun <- function(X) {
            0*X[,1] + 0.2 * (1 + X[,1])*X[,2] + ifelse(X[,3] > 0, log(2), log(10))
        }
        surv_sigma_fun <- function(X) 0*X[,1] + 0.25
        survival_generator <- LogNormalDistribution$new(mu_fun = surv_mu_fun, sigma_fun = surv_sigma_fun)
        
        cens_mu_fun <- function(X) 0*X[,1] + 2 + 0.5 * X[,1]
        cens_sigma_fun <- function(X) 0*X[,1] + 0.1
        censoring_generator <- LogNormalDistribution$new(mu_fun = cens_mu_fun, sigma_fun = cens_sigma_fun)
        ## Mixture distribution (covariate shift)
        surv_mu_fun_shift <- function(X) {
            0*X[,1] + 0.2 * (1 + X[,1])*X[,2] + log(10)
        }
        survival_generator_shift <- LogNormalDistribution$new(mu_fun = surv_mu_fun_shift, sigma_fun = surv_sigma_fun)        
        
    } else if(setting==1) {
        ## Setting from Sesia and Svetnik (difficult)
        ## Initialize the covariate model
        num_features <- 100
        covariate_generator <- function(num_samples) {
            matrix(runif(num_samples * num_features, 0, 1), nrow = num_samples)
        }
        ## Initialize the survival time distribution
        surv_mu_fun <- function(X) 0*X[,1] + (X[,2]>0.5) + (X[,3]<0.5) + (1-X[,1])^0.25
        surv_sigma_fun <- function(X) 0*X[,1] + 0.1 * (1-X[,1])
        survival_generator <- LogNormalDistribution$new(mu_fun = surv_mu_fun, sigma_fun = surv_sigma_fun)
        ## Initialize the censoring time distribution
        cens_mu_fun <- function(X) 0*X[,1] + (X[,2]>0.5) + (X[,3]<0.5) + (1-X[,1])^4 + 0.4 #+ 10
        cens_sigma_fun <- function(X) 0*X[,1] + 0.1 * X[,2]
        censoring_generator <- LogNormalDistribution$new(mu_fun = cens_mu_fun, sigma_fun = cens_sigma_fun)
        survival_generator_shift <- survival_generator
        
    } else if(setting==2) {
        ## Setting from Sesia and Svetnik (medium)
        ## Initialize the covariate model
        num_features <- 100
        covariate_generator <- function(num_samples) {
            matrix(runif(num_samples * num_features, 0, 1), nrow = num_samples)
        }
        ## Initialize the survival time distribution
        surv_mu_fun <- function(X) 0*X[,1] + X[,1]^0.25
        surv_sigma_fun <- function(X) 0*X[,1] + 0.1
        survival_generator <- LogNormalDistribution$new(mu_fun = surv_mu_fun, sigma_fun = surv_sigma_fun)
        ## Initialize the censoring time distribution
        cens_mu_fun <- function(X) 0*X[,1] + X[,1]^4 + 0.4
        cens_sigma_fun <- function(X) 0*X[,1] + 0.1
        censoring_generator <- LogNormalDistribution$new(mu_fun = cens_mu_fun, sigma_fun = cens_sigma_fun)
        survival_generator_shift <- survival_generator

    } else if(setting==3) {
        ## Setting from Sesia and Svetnik (easy)
        ## Initialize the covariate model
        num_features <- 100
        covariate_generator <- function(num_samples) {
            matrix(runif(num_samples * num_features, -1, 1), nrow = num_samples)
        }
        ## Initialize the survival time distribution
        surv_mu_fun <- function(X) 0*X[,1] + log(2) + 1 + 0.55*(X[,1]^2-X[,3]*X[,5])
        surv_sigma_fun <- function(X) 0*X[,1] + abs(X[,10]) + 1
        survival_generator <- LogNormalDistribution$new(mu_fun = surv_mu_fun, sigma_fun = surv_sigma_fun)
        ## Initialize the censoring time distribution
        cens_rate_fun <- function(X) 0*X[,1] + 0.4
        censoring_generator <- ExponentialDistribution$new(rate_fun = cens_rate_fun)
        survival_generator_shift <- survival_generator

    }

    return(list(num_features=num_features,
                covariate_generator=covariate_generator,
                survival_generator=survival_generator,
                censoring_generator=censoring_generator,
                survival_generator_shift=survival_generator_shift))

}

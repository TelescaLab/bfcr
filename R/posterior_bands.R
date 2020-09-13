#' Get means and upper/lower credible bands for each subject
#' 
#' Posterior processing function for computing means and credible bands. Refer 
#' to our manuscript for details on credible band calculation
#' 
#' @param mcmc_output Object returned from run_mcmc()
#' @param alpha_level Type I error
#' @export get_posterior_subject_bands
#' @return a tibble with columns subject ID, response, mean, lower limit,
#' upper limit, and time
#' @examples subject_bands <- get_posterior_subject_bands(mcmc_output)
#' subj <- 80
#' p1 <- subject_bands %>% 
#'   filter(id == subj) %>%
#'  ggplot() +
#'  geom_point(aes(x = time, y = response)) +
#'  geom_ribbon(aes(x = time, ymin = lower, ymax = upper), alpha = 0.3) +
#'  theme_bw()
#' p1
get_posterior_subject_bands <- function(mcmc_output, alpha_level = .05) {
  subject_summaries <- 
    get_posterior_subject_bands_cpp(mcmc_output, alpha_level)
  num_subjects <- dim(mcmc_output$data$response)[1]
  num_time_pts <- dim(mcmc_output$data$response)[2]
  response <- c(t(mcmc_output$data$response))
  mean <- c(subject_summaries$subject_means)
  lower <- c(subject_summaries$subject_lower)
  upper <- c(subject_summaries$subject_upper)
  id <- rep(1:num_subjects, each = num_time_pts)
  time <- rep(mcmc_output$data$time, times = num_subjects)
  
  posterior_subject_bands <- tibble(id = id,
                                    response,
                                    mean = mean,
                                    lower = lower,
                                    upper = upper,
                                    time = time)
}

#' Get posterior mean function
#' 
#' Process samples from an object returned from run_mcmc to obtain a posterior
#' mean with upper/lower bounds. The mean can be evaluated at a specific 
#' value of covariate.
#' 
#' @param mcmc_output object returned from run_mcmc
#' @param xi specified covariate value to obtain mean function
#' @param alpha_level Type I error 
#' @export get_posterior_means
#' @return A tibble with columns for lower bound, mean, upper bound, 
#' and user specified time grid
#' @examples # Intercept only model
#' xi <- c(1)
#' mean_bands <- get_posterior_means(mcmc_output, xi)
get_posterior_means <- function(mcmc_output, xi, alpha_level = .05) {
  mean_matrix <- get_posterior_means_cpp_correct(mcmc_output, xi, alpha_level)
  mean_tibble <- tibble(lower = mean_matrix[,1],
                        mean = mean_matrix[,2],
                        upper = mean_matrix[,3],
                        time  = mcmc_output$data$time)
}

#' Posterior inference for covariate-adjusted covariance function
#' 
#' @param mcmc_results mcmc object
#' @param eigenvals Number of eigenvalues to keep
#' @param zi Covariate vector of interest
#' @param alpha_level Type I error rate
#' @details Generates posterior inference for covariate adjusted
#'  eigenfunctions, surfaces, and magnitudes
#' @return
#' An R list containing the following elements 
#' 
#' \code{lower} A matrix containing the lower bound of the simultaneous 
#' credible band of eigenfunctions  
#' 
#' \code{mean} A matrix containing the posterior mean of eigenfunctions  
#' 
#' \code{upper} A matrix containing the upper bound of the simultaneous 
#' credible band of eigenfunctions  
#' 
#' \code{eigenval_intervals} A matrix containing a 1-alpha credible interval
#' for eigenvalues  
#' 
#' \code{eigenval_pve_intervals} A matrix containing a 1-alpha credible 
#' interval for relative eigenvalues  
#' 
#' \code{surface} Posterior covariance surface  
#' 
#' \code{magnitude} Total variance of the fitted covariance surface across, one
#' for each sample
#' 
#' \code{raw_magnitude} 1-alpha credible interval for total variance
#' @export get_posterior_eigen
get_posterior_eigen <- function(mcmc_results, eigenvals, zi, alpha_level=0.05) {
  post_eigen <- get_posterior_eigen_cpp_correct(mcmc_results, eigenvals, zi, alpha_level)
  eigenfunctions <- tibble(mean = c(post_eigen$mean_eigen),
                           lower = c(post_eigen$lower_eigen),
                           upper = c(post_eigen$upper_eigen),
                           number = rep(factor(1:eigenvals),
                                        each = length(post_eigen$time)),
                           time = rep(post_eigen$time, eigenvals))
  return(list(eigenfunctions = eigenfunctions,
         eigenvalues = post_eigen$eigenval_intervals,
         prop_var_explained = post_eigen$eigenval_pve_intervals,
         surface = post_eigen$surface,
         magnitude = post_eigen$magnitude,
         raw_magnitude = post_eigen$raw_magnitude,
         time = post_eigen$time))
}

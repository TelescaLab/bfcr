#' Get means and upper/lower credible bands for each subject
#' 
#' Posterior processing function for computing means and credible bands. Refer 
#' to our manuscript for details on credible band calculation
#' 
#' @param mcmc_output Object returned from run_mcmc()
#' @param alpha_level Specifies Type I error of credible band calculation
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
#' 
get_posterior_subject_bands <- function(mcmc_output, alpha_level = .05) {
  subject_summaries <- 
    get_posterior_predictive_bands69(mcmc_output, alpha_level)
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


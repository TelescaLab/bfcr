get_posterior_subject_bands <- function(mcmc_output, alpha_level = .05) {
  subject_summaries <- 
    get_posterior_predictive_bands69(mcmc_output, alpha_level)
  num_subjects <- dim(mcmc_output$data$response)[1]
  num_time_pts <- dim(mcmc_output$data$response)[2]
  id <- rep(rep(1:num_subjects, each = num_time_pts), 4)
  time <- rep(rep(mcmc_output$data$time, times = num_subjects), 4)
  label <- rep(c("actual", "subject_mean", "lower", "upper"),
               each = num_subjects * num_time_pts)
  value <- c(c(t(mcmc_output$data$response)),
             c(subject_summaries$subject_means),
             c(subject_summaries$subject_lower),
             c(subject_summaries$subject_upper))
  
  posterior_subject_bands <- tibble(id = id,
                                    label = label,
                                    value = value,
                                    time = time)
}
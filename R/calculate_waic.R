calculate_waic <- function(mcmc_output) {
  full_ind <- 0:(length(mcmc_output$data$time) -1)
  num_sub <- dim(mcmc_output$data$response)[1]
  observed_time <- vector(mode = "list", length = num_sub)
  num_post <- mcmc_output$control$iterations -
    mcmc_output$control$burnin
  waic_return <- matrix(0, nrow = num_post, 
                        ncol = dim(mcmc_output$data$response)[1])
  calc_observed_time <- function(i) {
    setdiff(full_ind, mcmc_output$data$missing_time[
      which(mcmc_output$data$missing_subjects == i-1)])
  }
  observed_time <- 1:num_sub %>% purrr::map(calc_observed_time)
  calculate_waic_cpp(mcmc_output, observed_time, waic_return)
  rel_n_eff <- loo::relative_eff(exp(waic_return), chain_id = rep(1, num_post))
  loo_result <- loo::loo(waic_return, r_eff = rel_n_eff)
  waic_result <- loo::waic(waic_return)
  list(loo_result = loo_result, waic_result = waic_result)
}
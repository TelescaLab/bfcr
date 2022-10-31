library(tidyverse)
library(mgcv)
library(future.apply)
library(bfcr)
library(pracma)

# Uncomment and change if working directory needs to be changed
#setwd("")
dir.create("Simulations")
dir.create("Simulations/Metrics")

compute_metrics <- function(this_seed) {
  file_name <- paste0(getwd(),
                      "/Simulation_parameters.RData")
  
  load(file_name)
  n <- 100
  k <- simparam$k
  times <- simparam$times
  beta <- simparam$beta
  lambda <- simparam$lambda
  time_basis <- simparam$time_basis
  x <- seq(from = 0, to = 1, length.out = n)
  df_x <- simparam$df_x
  x_basis <- smoothCon(s(x, k = df_x, m = 2, bs = "ps"), data.frame(x),
                       absorb.cons = TRUE)
  df_times <- simparam$df_times
  x_basis_fit <- smoothCon(s(x, k = df_x + 2, m = 2, bs = "ps"),
                           absorb.cons = TRUE, data.frame(x))
  time_basis_fit <- smoothCon(s(times, k = df_times, m = 2, bs = "ps"),
                              data.frame(times))
  design_mean <- cbind(1, x_basis_fit[[1]]$X)
  indices_mean <- c(1, 2, 2)
  time_penalty <- time_basis_fit[[1]]$S.scale * time_basis[[1]]$S[[1]]
  x_penalty <- x_basis_fit[[1]]$S.scale * x_basis_fit[[1]]$S[[1]]
  S_list <- list(x_penalty, time_penalty)
  S <- tensor.prod.penalties(S_list)
  penalties_mean <- list(time_penalty, S[[1]], S[[2]])
  k_fit <- k + 2
  evaluate_basis <- function(x_basis, x) {
    c(1, PredictMat(x_basis, data.frame(x)))
  }
  design_mean_truth <- cbind(1, x_basis[[1]]$X)
  set.seed(this_seed)
  j <- 1
  for (j in 1:4) {
    if (j == 1) {
      print(j)
      lambda_truth <- lambda
      design_var_truth <- cbind(1, x_basis[[1]]$X)
      design_var <- cbind(1, x_basis_fit[[1]]$X)
      evaluate_basis_var <- function(x_basis, x) {
        c(1, PredictMat(x_basis, data.frame(x)))
      }
      evaluate_basis_var_truth <- function(x_basis, x) {
        c(1, PredictMat(x_basis, data.frame(x)))
      }
      penalties_var <- list(time_penalty, S[[1]], S[[2]])
      indices_var <- c(1, 2, 2)
      str_file <- "ca"
    }
    else if (j == 2) {
      {
        print(j)
        lambda_truth <- lambda
        design_var_truth <- cbind(1, x_basis[[1]]$X)
        design_var <- cbind(rep(1, n))
        evaluate_basis_var <- function(x_basis, x) {
          c(1)
        }
        evaluate_basis_var_truth <- function(x_basis, x) {
          c(1, PredictMat(x_basis, data.frame(x)))
        }
        penalties_var <- list(time_penalty)
        indices_var <- c(1)
        str_file <- "base"
      }
    } else if (j == 3) {
      {
        print(j)
        lambda_truth <- array(lambda[,1,], c(df_times, 1, k))
        design_var_truth <- cbind(rep(1, n))
        design_var <- cbind(1, x_basis_fit[[1]]$X)
        evaluate_basis_var <- function(x_basis, x) {
          c(1, PredictMat(x_basis, data.frame(x)))
        }
        evaluate_basis_var_truth <- function(x_basis, x) {
          c(1)
        }
        penalties_var <- list(time_penalty, S[[1]], S[[2]])
        indices_var <- c(1, 2, 2)
        str_file <- "nocov"
      }
    } else if (j == 4) {
      print(j)
      lambda_truth <- array(lambda[,1,], c(df_times, 1, k))
      design_var_truth <- cbind(rep(1, n))
      design_var <- cbind(rep(1, n))
      evaluate_basis_var <- function(x_basis, x) {
        c(1)
      }
      evaluate_basis_var_truth <- function(x_basis, x) {
        c(1)
      }
      penalties_var <- list(time_penalty)
      indices_var <- c(1)
      str_file <- "nocovbase"
    }
    print(str_file)
    eta <- rnorm(n * k) %>% array(dim = c(n, k))
    Y <- time_basis[[1]]$X %*% beta %*% t(design_mean_truth) %>% t()
    for (kp in 1:k) {
      Y <- Y + time_basis[[1]]$X %*% lambda_truth[,,kp] %*%
        t(design_var_truth) %*% diag(eta[,kp]) %>% t()
    }
    Y_errors <- Y + rnorm(n * length(times), sd = .20) %>% matrix(n, length(times))
    mcmc_results <- run_mcmc(Y_errors, design_mean, design_var,
                             time_basis_fit[[1]]$X, times, penalties_mean,
                             penalties_var, indices_mean, indices_var,
                             k_fit, 50000, 25000, 5)
    subject_bands <- get_posterior_subject_bands(mcmc_results)
    subject_bands <- subject_bands %>% mutate(truth = c(t(Y))) %>%
      mutate(in_bounds = (lower < truth) & (upper & truth))
    subj_coverage <- subject_bands %>%
      summarise(sum(in_bounds) / n()) %>% pull()
    subj_width <- subject_bands %>%
      summarise(mean(upper-lower)) %>% pull()
    subj_error <- subject_bands %>%
      group_by(id) %>%
      summarise(rel_error =
                  trapz(times, 100 * (mean - truth)^2 / trapz(times, truth^2))) %>%
      summarise(mean(rel_error)) %>% pull()
    
    evaluated_points <- seq(from = .1, to = .9, length.out = 10)
    mean_tibble <- tibble(est_mean = numeric(10 * length(times)),
                          lower = numeric(10 * length(times)),
                          upper = numeric(10 * length(times)),
                          truth = numeric(10 * length(times)),
                          x = numeric(10 * length(times)))
    start <- 1
    
    for (i in 1:10) {
      end <- i * length(times)
      xi <- evaluate_basis(x_basis_fit[[1]], evaluated_points[i])
      mean_bands <- get_posterior_means(mcmc_results, xi)
      mean_tibble[start:end, ]$est_mean <- mean_bands$mean
      mean_tibble[start:end, ]$lower <- mean_bands$lower
      mean_tibble[start:end, ]$upper <- mean_bands$upper
      mean_tibble[start:end, ]$truth <- time_basis[[1]]$X %*%
        beta %*% evaluate_basis(x_basis[[1]], evaluated_points[i]) %>% c()
      mean_tibble[start:end, ]$x <- rep(evaluated_points[i], length(times))
      start <- end + 1
    }
    mean_coverage <- mean_tibble %>%
      mutate(in_bounds = (lower < truth) & (upper & truth)) %>%
      summarise(sum(in_bounds) / n()) %>%
      pull()
    mean_width <- mean_tibble %>%
      summarise(mean(upper-lower)) %>% pull()
    mean_error <- mean_tibble %>%
      group_by(x) %>%
      summarise(rel_error =
                  trapz(times, 100 * (est_mean - truth)^2 / trapz(times, truth^2))) %>%
      summarise(mean(rel_error)) %>% pull()
    
    covariance_tibble <- tibble(mean = numeric(10 * length(times) * length(times)),
                                lower = numeric(10 * length(times) * length(times)),
                                upper = numeric(10 * length(times) * length(times)),
                                truth = numeric(10 * length(times) * length(times)),
                                row_num = numeric(10 * length(times) * length(times)),
                                col_num = numeric(10 * length(times) * length(times)),
                                x = numeric(10 * length(times) * length(times)))
    start <- 1
    
    for (i in 1:10) {
      end <- i * length(times) * length(times)
      zi <- evaluate_basis_var(x_basis_fit[[1]], evaluated_points[i])
      cov_summ <- get_posterior_covariance(mcmc_results, zi, .05)
      covariance_tibble$mean[start:end] <- cov_summ$mean
      covariance_tibble$lower[start:end] <- cov_summ$lower
      covariance_tibble$upper[start:end] <- cov_summ$upper
      covariance_tibble$row_num[start:end] <- rep(1:length(times), each = length(times))
      covariance_tibble$col_num[start:end] <- rep(1:length(times), length(times))
      temp_truth <- numeric(length(times) * length(times))
      for (kp in 1:k) {
        temp_truth <- temp_truth + c(tcrossprod(time_basis[[1]]$X %*% lambda_truth[,,kp] %*% evaluate_basis_var_truth(x_basis[[1]], evaluated_points[i])))
      }
      covariance_tibble$truth[start:end] <- temp_truth
      covariance_tibble$x[start:end] <- rep(evaluated_points[i],
                                            length(times) * length(times))
      start <- end + 1
    }
    covariance_coverage <- covariance_tibble %>%
      filter(row_num >= col_num) %>%
      mutate(in_bounds = (lower < truth) & (upper > truth)) %>%
      summarise(mean(in_bounds)) %>% pull()
    covariance_width <- covariance_tibble %>%
      filter(row_num >= col_num) %>%
      summarise(mean(upper - lower)) %>% pull()
    covariance_error <- covariance_tibble %>%
      mutate(groupx = (1:n() - 1) %/% length(times)) %>%
      group_by(x, groupx) %>%
      summarize(error = trapz(times, (mean - truth)^2), normalize = trapz(times, truth^2)) %>%
      ungroup(groupx) %>%
      summarize(normalizing = trapz(times, normalize), error = trapz(times, error),
                rise = 100 * error / normalizing) %>%
      summarize(rise = mean(rise)) %>% pull()
    metrics <- list(subj_coverage = subj_coverage, subj_error = subj_error,
                    subj_width = subj_width, mean_coverage = mean_coverage,
                    mean_error = mean_error, mean_width = mean_width,
                    covariance_error = covariance_error,
                    covariance_width = covariance_width,
                    covariance_coverage = covariance_coverage)
    file_name <- paste0(getwd(),
                        "/Metrics/n", n, "_", str_file, "_seed", this_seed,
                        ".RData")
    save(metrics, file = file_name)
  }
}

run_100 <- function(){
  ### Limit the number of CPU cores used (can change depending on computer)
  ncpu <- min(4, availableCores())
  # 
  plan(multisession, workers = ncpu)
  already_ran <- dir(paste0(getwd(), "/Metrics"))
  to_run <- which(!paste0("n100_nocovbase_seed", 1:300, ".RData") %in% already_ran)
  seeds <- to_run
  future_lapply(seeds[1:8], function(this_seed) compute_metrics(this_seed))
  more <- F
  if(length(seeds) > 8){
    more <- T
  }
  return(cat(as.numeric(more)))
}

run_100()
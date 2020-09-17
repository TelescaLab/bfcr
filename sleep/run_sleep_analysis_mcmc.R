# Running the MCMC routine for all subjects in sleep study.
# Using principal components by conditional expectation (PACE), 8 components
# are needed to explain 99% of variability.
# This model will have 12 components (a conservative estimate)  
library(BayesianConditionalFPCA)
library(tidyverse)
library(mgcv)
library(spam)
library(plotly)
library(dlnm)
library(fdapace)

load(
  paste0("/Users/johnshamshoian/Documents/R_projects/",
         "bfcr/sleep/data/",
         "relative_psd.RData"))

sleep_tabulated <- read.csv(
  paste0("/Users/johnshamshoian/Documents/R_projects/",
         "bfcr/sleep/tabulated_data/", 
         "shhs1-dataset-0.15.0.csv"),
  stringsAsFactors = FALSE)

num_epochs <- 120
id_range <- 200001:206000

sleep_tabulated_filtered <- sleep_tabulated %>%
  filter(EEG1qual == 4) %>%
  select(nsrrid, age_s1) %>%
  drop_na()

sleep_data <- inner_join(sleep_tabulated_filtered, relative_psd)

sleep_data_filtered <- sleep_data %>%
  group_by(nsrrid) %>%
  filter(n() >= num_epochs, epoch <= num_epochs,
         nsrrid %in% id_range) %>%
  ungroup()

num_subjects <- as.numeric(summarise(
  sleep_data_filtered, n_distinct(nsrrid)))

epoch_grid <- 1:num_epochs
epoch_df <- ceiling(20 / 100 * num_epochs)
age_grid <- sleep_data_filtered %>%
  group_by(nsrrid) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(age_s1)
age_df <- ceiling(10 / 100 * (max(age_grid) - min(age_grid)))

epoch_basis_object <- smoothCon(
  s(epoch_grid, bs = "ps", k = epoch_df, m = 2),
  data = data.frame(epoch_grid),
  absorb.cons = TRUE
)

epoch_basis <- cbind(1, epoch_basis_object[[1]]$X)    
epoch_penalty <- matrix(0, epoch_df, epoch_df)
epoch_penalty[2:epoch_df, 2:epoch_df] <- epoch_basis_object[[1]]$S[[1]]

evaluate_spline <- function(smoothCon_list, age) {
  c(1, PredictMat(smoothCon_list[[1]], data = data.frame(Age = age)))
}

# Age adjusted, using mgcv
{
  Age <- age_grid$age_s1
  age_basis <- smoothCon(s(Age, bs = "ps", k = 6, m = 2),
                         data = data.frame(Age),
                         absorb.cons = TRUE)
  design_mean <- cbind(1, age_basis[[1]]$X)
  design_var <- cbind(1, age_basis[[1]]$X)
  model_penalties <- tensor.prod.penalties(list(age_basis[[1]]$S[[1]], epoch_penalty))
  mean_indices <- c(1, 2, 2)
  var_indices <- c(1, 2, 2)
  mean_penalty <- list(epoch_penalty, model_penalties[[1]], model_penalties[[2]])
  var_penalty <- list(epoch_penalty, model_penalties[[1]], model_penalties[[2]])
  new_ages <- 40:79
}

response <- t(matrix(sleep_data_filtered$psd,
                     nrow = num_epochs,
                     ncol = num_subjects))

k <- 12
iter <- 20000
burnin <- 5000
thin <- 10

mcmc_results <- run_mcmc(response, design_mean,
                         design_var, epoch_basis,
                         epoch_grid,
                         mean_penalty, var_penalty,
                         mean_indices, var_indices,
                         k, iter, burnin, thin = thin,
                         var = "unequal")
saveRDS(
  mcmc_results, 
  file = paste0("/Users/johnshamshoian/Documents/R_projects/", 
                "bfcr/sleep/mcmc_output/mcmc_results.rds")
)

### NOT RUN

# Not age adjusted
if (FALSE)
{
  design_mean <- cbind(rep(1, num_subjects))
  design_var <- cbind(rep(1, num_subjects))
  mean_indices <- c(1)
  var_indices <- c(1)
  mean_penalty <- list(epoch_penalty)
  var_penalty <- list(epoch_penalty)
}

# Age adjusted, custom splines
if (FALSE)
{
  age_basis <- ps(age_grid$age_s1, df = age_df, intercept = FALSE)
  age_penalty <- attr(age_basis, "S")
  if (!(attr(age_basis, "intercept"))) {
    model_penalties <- tensor.prod.penalties(list(age_penalty, epoch_penalty))
    
    mean_indices <- c(1, 2, 2)
    var_indices <- c(1, 2, 2)
    var_penalty <- list(epoch_penalty, model_penalties[[1]], model_penalties[[2]])
    # var_penalty <- list(epoch_penalty)
    # var_penalty <- list(epoch_penalty, diag(180), diag(180))
    mean_penalty <- list(epoch_penalty, model_penalties[[1]], model_penalties[[2]])
    # mean_penalty <- list(epoch_penalty, diag(180), diag(180))
    design_mean <- cbind(1, age_basis)
    design_var <- cbind(1, age_basis)
    # design_var <- cbind(rep(1, num_subjects))
  } else {
    model_penalties <- tensor.prod.penalties(list(age_penalty, epoch_penalty))
    
    mean_indices <- c(1, 1)
    var_indices <- c(1, 1)
    mean_penalty <- model_penalties
    var_penalty <- model_penalties
    design_mean <- age_basis
    design_var <- age_basis
  }
}
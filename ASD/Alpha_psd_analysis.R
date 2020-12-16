library(BayesianConditionalFPCA)
library(mgcv)
library(tidyverse)
library(fdapace)
library(plotly)
source("/Users/johnshamshoian/Documents/R_projects/bfcr/ASD/Peak_Alpha_Data_Transfer.R")
num_subjects <- pa.dat %>%
  summarise(n_distinct(ID)) %>%
  pull()
freq_grid <- pa.dat %>%
  distinct(func) %>%
  pull()
age_grid <- pa.dat %>% 
  group_by(ID) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(Age) %>%
  pull()
group_grid <- pa.dat %>%
  group_by(ID) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  mutate(group = group - 1) %>%
  select(group) %>%
  pull()

num_times <- pa.dat %>%
  filter(ID == 2, reg == 1) %>%
  count() %>%
  pull()
response <- pa.dat %>%
  filter(reg == 15) %>%
  select(y) %>%
  pull() %>%
  matrix(nrow = num_times, ncol = num_subjects) %>%
  t()

freq_grid <- seq(from = 6, to = 14, by = .25)
freq_basis_obj <- smoothCon(s(freq_grid, k = 12, bs = "ps", m = 2),
                        data.frame(freq_grid),
                        absorb.cons = FALSE)
freq_basis <- freq_basis_obj[[1]]$X
freq_penalty <- freq_basis_obj[[1]]$S[[1]] * freq_basis_obj[[1]]$S.scale

{
  age_basis_obj <- smoothCon(s(age_grid, k = 6, bs = "ps", by = as.factor(group_grid)),
                             data = data.frame(age_grid, group_grid),
                             absorb.cons = TRUE)
  age_penalty <- age_basis_obj[[1]]$S[[1]] * age_basis_obj[[1]]$S.scale
  model_penalties <- tensor.prod.penalties(list(age_penalty, freq_penalty))
  design_mean <- cbind(1, group_grid,
                       age_basis_obj[[1]]$X + age_basis_obj[[2]]$X,
                       age_basis_obj[[2]]$X)
  design_var <- cbind(1, group_grid,
                       age_basis_obj[[1]]$X + age_basis_obj[[2]]$X,
                       age_basis_obj[[2]]$X)
  mean_penalty <- list(freq_penalty, freq_penalty,
                       model_penalties[[1]], model_penalties[[2]],
                       model_penalties[[1]], model_penalties[[2]])
  var_penalty <- list(freq_penalty, freq_penalty,
                      model_penalties[[1]], model_penalties[[2]],
                      model_penalties[[1]], model_penalties[[2]])
  var_indices <- c(1, 2, 3, 3, 4, 4)
  mean_indices <- c(1, 2, 3, 3, 4, 4)
  evaluate_basis <- function(obj, age, group) {
    spline_part1 <- 
      PredictMat(obj[[1]], data.frame(age_grid = age, group_grid = group)) +
      PredictMat(obj[[2]], data.frame(age_grid = age, group_grid = group))
      
    spline_part2 <- 
      PredictMat(obj[[2]], data.frame(age_grid = age, group_grid = group))
    c(1, group, spline_part1, spline_part2)
  }
}

# {
#   age_basis_obj <- smoothCon(s(age_grid, k = 6, bs = "ps"),
#                              data = data.frame(age_grid, group_grid),
#                              absorb.cons = TRUE)
#   age_penalty <- age_basis_obj[[1]]$S[[1]]
#   model_penalties <- tensor.prod.penalties(list(age_penalty, freq_penalty))
#   design_mean <- cbind(1, group_grid, age_basis_obj[[1]]$X)
#   design_var <- cbind(1, group_grid, age_basis_obj[[1]]$X)
#   mean_penalty <- list(freq_penalty, freq_penalty,
#                        model_penalties[[1]], model_penalties[[2]])
#   var_penalty <- list(freq_penalty, freq_penalty,
#                       model_penalties[[1]], model_penalties[[2]])
#   var_indices <- c(1, 2, 3, 3)
#   mean_indices <- c(1, 2, 3, 3)
#   evaluate_basis <- function(age_basis_obj, age, group) {
#     spline_part <- PredictMat(age_basis_obj[[1]],
#                               data.frame(age_grid = age, group_grid = group))
#     c(1, group, spline_part)
#   }
# }

# {
#   design_mean <- cbind(rep(1, num_subjects))
#   design_var <- cbind(rep(1, num_subjects))
#   mean_indices <- c(1)
#   var_indices <- c(1)
#   mean_penalty <- list(freq_penalty)
#   var_penalty <- list(freq_penalty)
# }

k <- 8
iter <- 10000
burnin <- 5000
thin <- 20

mcmc_results <- run_mcmc(response, design_mean,
                         design_var, freq_basis,
                         freq_grid,
                         mean_penalty, var_penalty,
                         mean_indices, var_indices,
                         k, iter, burnin, thin = thin,
                         var = "pooled")
file_name <- paste0("/Users/johnshamshoian/Documents/R_projects/bfcr/", 
                    "Peak Alpha Data/mcmc_results.rds")
saveRDS(mcmc_results, file = file_name)
mcmc_results <- readRDS(file_name)
freq_grid <- seq(from = 6, to = 14, by = .05)
freq_basis_obj <- smoothCon(s(freq_grid, k = 12, bs = "ps", m = 2),
                            data.frame(freq_grid),
                            absorb.cons = FALSE)
mcmc_results$data$time <- freq_grid
mcmc_results$data$basis <- freq_basis_obj[[1]]$X
num_times <- length(freq_grid)
new_age_grid <- c(50, 70, 90, 110)
mean_tibble <- tibble(mean = numeric(num_times * length(new_age_grid) * 2),
                      lower = numeric(num_times * length(new_age_grid) * 2),
                      upper = numeric(num_times * length(new_age_grid) * 2),
                      age = rep(new_age_grid, each = num_times * 2),
                      group = factor(rep(rep(c(0, 1), each = num_times),
                                  length(new_age_grid))),
                      freq = numeric(2 * num_times * length(new_age_grid)))
first_index <- 1
for (i in seq_along(new_age_grid)) {
  for (j in c(0, 1)) {
    xi <- evaluate_basis(age_basis_obj, new_age_grid[i], j)
    mean_bands <- get_posterior_means(mcmc_results, xi)
    last_index <- first_index + num_times - 1
    mean_tibble$mean[first_index:last_index] <- mean_bands$mean
    mean_tibble$lower[first_index:last_index] <- mean_bands$lower
    mean_tibble$upper[first_index:last_index] <- mean_bands$upper
    mean_tibble$freq[first_index:last_index] <- freq_grid
    first_index <- last_index + 1
  }
}
mean_tibble %>%
  ggplot(aes(x = freq)) +
  geom_line(aes(y = mean, lty = group), size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = .3, show.legend = FALSE) +
  facet_wrap(~ age, nrow = 1) +
  labs(x = expression(omega), y = expression(paste(mu,"(", omega, ", ", x, ")"))) +
  scale_linetype_manual(values = c(1,4), name = "Group", labels = c("TD", "ASD")) +
  scale_fill_manual(values = c("black", "black")) + 
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 3))    +
  theme(axis.title.y = element_text(margin = margin(r = 10, b = 0, l = 0)))


N <- 20
lower_age <- quantile(age_grid, .1)
upper_age <- quantile(age_grid, .9)
new_age_grid <- seq(from = lower_age, to = upper_age, length.out = N)

eigen_surface_td <- matrix(0, N, num_times)
magnitude_td <- matrix(0, N, 3)
for (i in 1:N) {
  print(i)
  zi <- evaluate_basis(age_basis_obj, new_age_grid[i], 0)
  eigen_bands <- get_posterior_eigen(mcmc_results, 1, zi)
  eigen_surface_td[i,] <- eigen_bands$eigenfunctions$mean
  if (eigen_surface_td[i, 1] < 0) {
    eigen_surface_td[i,] <- -1 * eigen_surface_td[i,]
  }
  magnitude_td[i, ] <- eigen_bands$magnitude
}


eigen_surface_asd <- matrix(0, N, num_times)
magnitude_asd <- matrix(0, N, 3)
for (i in 1:N) {
  zi <- evaluate_basis(age_basis_obj, new_age_grid[i], 1)
  eigen_bands <- get_posterior_eigen(mcmc_results, 1, zi)
  eigen_surface_asd[i,] <- eigen_bands$eigenfunctions$mean
  if (eigen_surface_asd[i, 1] < 0) {
    eigen_surface_asd[i,] <- -1 * eigen_surface_asd[i,]
  }
  magnitude_asd[i, ] <- eigen_bands$magnitude
}

x <- freq_grid
y <- new_age_grid
data <- expand.grid(X=x, Y=y, G = factor(c(0, 1)))
data$Z <- c(c(t(eigen_surface_td)), c(t(eigen_surface_asd)))
group.labs <- c("TD", "ASD")
names(group.labs) <- c(0, 1)
ggplot(data, aes(X, Y, fill = Z)) +
  geom_tile() +
  labs(x = expression(paste(omega)), y = "Age (months)") +
  facet_wrap(~ G, labeller = labeller(G = group.labs)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines")) +
  scale_fill_gradient(name = "Eigenfunction 1", low = "black", high = "white")


freq_basis <- freq_basis_obj[[1]]$X
intercept_effect <- numeric(num_times)
age_effect <- numeric(num_times)
group_effect <- numeric(num_times)
interaction_effect <- numeric(num_times)
for (i in seq_along(new_age_grid)) {
  age_zi <- evaluate_basis(age_basis_obj, new_age_grid[i], 0)
  age_zi[1] <- 0
  group_zi <- c(0,1, rep(0, 10))
  interaction_zi <- c(rep(0, 7), age_zi[3:7])
  for (iter in 5001:10000) {
    for (kp in 1:k) {
      intercept_effect <- intercept_effect + c((freq_basis %*% mcmc_results$samples$lambda[[iter]][,1,kp])^2)
      age_effect <- age_effect + c((freq_basis %*% mcmc_results$samples$lambda[[iter]][,,kp] %*% age_zi)^2)
      group_effect <- group_effect + c((freq_basis %*% mcmc_results$samples$lambda[[iter]][,,kp] %*% group_zi)^2)
      interaction_effect <- interaction_effect + c((freq_basis %*% mcmc_results$samples$lambda[[iter]][,,kp] %*% interaction_zi)^2)
    }
  }
}
intercept_effect <- intercept_effect / (5000 * length(new_age_grid))
age_effect <- age_effect / (5000 * length(new_age_grid))
group_effect <- group_effect / (5000 * length(new_age_grid))
interaction_effect <- interaction_effect / (5000 * length(new_age_grid))
var_effects <- data.frame(effect = c(intercept_effect,
                                     age_effect,
                                     group_effect,
                                     interaction_effect),
                          freq = rep(freq_grid, 4),
                          type = factor(rep(1:4, each = length(freq_grid))))
type.labs <- c("Baseline", "Age effect", "Group effect", "Interaction effect")
names(type.labs) <- c(1, 2, 3, 4)
var_effects %>%
  ggplot(aes(x = freq)) +
  geom_line(aes(y = effect, linetype = type)) +
  labs(x = expression(paste(omega)), y = substitute(paste(italic('g(t, x)'))),
       linetype = "Covariate effect") +
  scale_linetype_manual(labels = c("Baseline", "Age effect", "Group effect", "Interaction effect"), values = c(1,3,5,6)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

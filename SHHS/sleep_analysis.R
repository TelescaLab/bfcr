library(BayesianConditionalFPCA)
library(tidyverse)
library(mgcv)
library(spam)
library(plotly)
library(dlnm)
library(fdapace)

load(paste0("/Users/johnshamshoian/Documents/R_projects/",
            "bfcr/sleep/data/",
            "relative_psd.RData"))
sleep_tabulated <- read.csv(
  paste0("/Users/johnshamshoian/Documents/R_projects/",
         "bfcr/sleep/tabulated_data/",
         "shhs1-dataset-0.15.0.csv"), stringsAsFactors = FALSE)

num_epochs <- 240
id_range <- 200001:206000

# Filter sleep delta power spectral density waveforms
sleep_data <- sleep_tabulated %>%
  filter(EEG1qual >= 3) %>%
  select(nsrrid, age_s1, HTNDerv_s1) %>%
  drop_na()
sleep_data <- inner_join(sleep_data, relative_psd)
sleep_data <- sleep_data %>%
  group_by(nsrrid) %>%
  filter(n() >= num_epochs, epoch <= num_epochs,
         nsrrid %in% id_range) %>%
  ungroup()


num_subjects <- as.numeric(summarise(
  sleep_data, n_distinct(nsrrid)))
epoch_grid <- 1:num_epochs

# Use one basis function every ten epochs
epoch_df <- ceiling(10 / 100 * num_epochs)

covariate_grid <- sleep_data %>%
  group_by(nsrrid) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(-nsrrid, -epoch, -psd) %>%
  as.data.frame()

epoch_smooth <- smoothCon(s(epoch_grid, k = epoch_df, bs = "ps", m = 2),
                          data.frame(epoch_grid))
epoch_basis <- epoch_smooth[[1]]$X
epoch_penalty <- epoch_smooth[[1]]$S[[1]] * epoch_smooth[[1]]$S.scale

{
  # Set up smoothing p-splines over age, by hypertension diagnostic status
  age_basis_obj <- smoothCon(s(age_s1, k = 8, bs = "ps", m = 2,
                               by = as.factor(HTNDerv_s1)),
                             data = covariate_grid,
                             absorb.cons = TRUE)
  age_penalty <- age_basis_obj[[1]]$S[[1]] * age_basis_obj[[1]]$S.scale
  model_penalties <- tensor.prod.penalties(list(age_penalty, epoch_penalty))
  X <- cbind(1, covariate_grid$HTNDerv_s1,
                       age_basis_obj[[1]]$X + age_basis_obj[[2]]$X,
                       age_basis_obj[[2]]$X)
  penalties <- list(epoch_penalty, epoch_penalty,
                       model_penalties[[1]], model_penalties[[2]],
                       model_penalties[[1]], model_penalties[[2]])
  
  # Set up smoothing parameters over epochs and age
  indices <- c(1, 2, 3, 3, 4, 4)
  evaluate_basis <- function(obj, age, group) {
    spline_part1 <- 
      PredictMat(obj[[1]], data.frame(age_s1 = age, HTNDerv_s1 = group)) +
      PredictMat(obj[[2]], data.frame(age_s1 = age, HTNDerv_s1 = group))
    
    spline_part2 <- 
      PredictMat(obj[[2]], data.frame(age_s1 = age, HTNDerv_s1 = group))
    c(1, group, spline_part1, spline_part2)
  }
}

response <- t(matrix(sleep_data$psd,
                     nrow = num_epochs,
                     ncol = num_subjects))

### Need 9 components to explain 99% variability, not adjusting for age
k <- 12
iter <- 10000
burnin <- 5000
thin <- 20
rm(sleep_data)
rm(relative_psd)
rm(sleep_tabulated)
mcmc_results <- run_mcmc(response, X,
                         X, epoch_basis,
                         epoch_grid,
                         penalties, penalties,
                         indices, indices,
                         k, iter, burnin, thin,
                         var = "pooled")

file_name <- paste0("/Users/johnshamshoian/Documents/R_projects/bfcr/sleep/",
                    "mcmc_output/mcmc_results.RData")
save(mcmc_results, file = file_name)
load(file_name)
new_age_grid <- c(50, 60, 70, 80)
mean_tibble <- tibble(mean = numeric(num_epochs * length(new_age_grid) * 2),
                      lower = numeric(num_epochs * length(new_age_grid) * 2),
                      upper = numeric(num_epochs * length(new_age_grid) * 2),
                      age = rep(new_age_grid, each = num_epochs * 2),
                      group = factor(rep(rep(c(0, 1), each = num_epochs),
                                         length(new_age_grid))),
                      epoch = numeric(2 * num_epochs * length(new_age_grid)))
first_index <- 1
for (i in seq_along(new_age_grid)) {
  for (j in c(0, 1)) {
    xi <- evaluate_basis(age_basis_obj, new_age_grid[i], j)
    mean_bands <- get_posterior_means(mcmc_results, xi)
    last_index <- first_index + num_epochs - 1
    mean_tibble$mean[first_index:last_index] <- mean_bands$mean
    mean_tibble$lower[first_index:last_index] <- mean_bands$lower
    mean_tibble$upper[first_index:last_index] <- mean_bands$upper
    mean_tibble$epoch[first_index:last_index] <- epoch_grid
    first_index <- last_index + 1
  }
}

mean_tibble %>%
  ggplot(aes(x = epoch)) +
  geom_line(aes(y = mean, lty = group), size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = group),
              alpha = .3, show.legend = FALSE) +
  facet_wrap(~ age, nrow = 1) +
  labs(x = "30 second epoch", y = expression(paste(mu,"(t, x)"))) +
  scale_linetype_manual(values = c(1,4), name = "Group", labels = c("No hypertension", "Hypertension")) +
  scale_fill_manual(values = c("black", "black")) +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 3))    +
  theme(axis.title.y = element_text(margin = margin(r = 10, b = 0, l = 0)))

# Can look at heatmaps of delta spectral power. Not included in paper.
N <- 20
lower_age <- quantile(covariate_grid$age_s1, .1)
upper_age <- quantile(covariate_grid$age_s1, .9)
new_age_grid <- seq(from = lower_age, to = upper_age, length.out = N)
mean_surface_no <- matrix(0, N, num_epochs)
mean_surface_yes <- matrix(0, N, num_epochs)
for (i in 1:N) {
  xi <- evaluate_basis(age_basis_obj, new_age_grid[i], 0)
  mean_bands <- get_posterior_means(mcmc_results, xi)
  mean_surface_no[i,] <- mean_bands$mean
  xi <- evaluate_basis(age_basis_obj, new_age_grid[i], 1)
  mean_bands <- get_posterior_means(mcmc_results, xi)
  mean_surface_yes[i,] <- mean_bands$mean
}
x <- epoch_grid
y <- new_age_grid
data <- expand.grid(X=x, Y=y, G = factor(c(0, 1)))
data$Z <- c(c(t(mean_surface_no)), c(t(mean_surface_yes)))
group.labs <- c("No hypertension", "Hypertension")
names(group.labs) <- c(0, 1)
ggplot(data, aes(X, Y, fill = Z)) +
  geom_tile() +
  labs(x = "30 second epoch", y = "Age (years)") +
  facet_wrap(~ G, labeller = labeller(G = group.labs)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines")) +
  scale_fill_gradient(name = expression(paste(mu,"(t, x)")), low = "black", high = "white")

eigen_surface_no <- matrix(0, N, num_epochs)
eigen_surface_yes <- matrix(0, N, num_epochs)
for (i in 1:N) {
  print(i)
  zi <- evaluate_basis(age_basis_obj, new_age_grid[i], 0)
  eigen_bands <- get_posterior_eigen(mcmc_results, 1, zi)
  eigen_surface_no[i,] <- eigen_bands$eigenfunctions  %>% filter(number == 1) %>% select(mean) %>% pull()
  zi <- evaluate_basis(age_basis_obj, new_age_grid[i], 1)
  eigen_bands <- get_posterior_eigen(mcmc_results, 1, zi)
  eigen_surface_yes[i,] <- eigen_bands$eigenfunctions %>% filter(number == 1) %>% select(mean) %>% pull()
  if (eigen_surface_no[i, 150] < 0) {
    eigen_surface_no[i,] <- -1 * eigen_surface_no[i,]
  }
  if (eigen_surface_yes[i, 150] < 0) {
    eigen_surface_yes[i,] <- -1 * eigen_surface_yes[i,]
  }
}

# for (i in 2:N) {
#   if (sum(eigen_surface_yes[i,] + eigen_surface_yes[1,])^2 <
#       sum(eigen_surface_yes[i,] - eigen_surface_yes[1,])^2) {
#     print("flipped")
#     eigen_surface_yes[i,] <- -1 * eigen_surface_yes[i,]
#   }
# }
# 
data <- expand.grid(X=x, Y=y, G = factor(c(0, 1)))
data$Z <- c(c(t(eigen_surface_no)), c(t(eigen_surface_yes)))
group.labs <- c("No hypertension", "Hypertension")
names(group.labs) <- c(0, 1)
ggplot(data, aes(X, Y, fill = Z)) +
  geom_tile() +
  labs(x = "30 second epoch", y = "Age (years)") +
  facet_wrap(~ G, labeller = labeller(G = group.labs)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines")) +
  scale_fill_gradient(name = "Eigenfunction 1", low = "black", high = "white")
# 
plot(eigen_bands$eigenfunctions$mean, type = "l")
lines(eigen_bands$eigenfunctions$lower)
lines(eigen_bands$eigenfunctions$upper)
matplot(t(eigen_surface_yes), type = "l")
intercept_effect <- numeric(num_epochs)
age_effect <- numeric(num_epochs)
group_effect <- numeric(num_epochs)
interaction_effect <- numeric(num_epochs)
group_zi <- c(0,1, rep(0, 14))
for (iter in 5001:10000) {
  for (kp in 1:k) {
    intercept_effect <- intercept_effect + (epoch_basis%*%mcmc_results$samples$lambda[[iter]][,1,kp])^2
    group_effect <- group_effect + c((epoch_basis %*% mcmc_results$samples$lambda[[iter]][,,kp] %*% group_zi)^2)
  }
}
for (i in seq_along(new_age_grid)) {
  age_zi <- evaluate_basis(age_basis_obj, new_age_grid[i], 0)
  age_zi[1] <- 0
  interaction_zi <- c(rep(0, 9), age_zi[3:9])
  for (iter in 5001:10000) {
    for (kp in 1:k) {
      age_effect <- age_effect + c((epoch_basis %*% mcmc_results$samples$lambda[[iter]][,,kp] %*% age_zi)^2)
      interaction_effect <- interaction_effect + c((epoch_basis %*% mcmc_results$samples$lambda[[iter]][,,kp] %*% interaction_zi)^2)
    }
  }
}

intercept_effect <- intercept_effect / 5000
age_effect <- age_effect / (5000 * length(new_age_grid))
group_effect <- group_effect / 5000
interaction_effect <- interaction_effect / (5000 * length(new_age_grid))
var_effects <- data.frame(effect = c(intercept_effect,
                                     age_effect,
                                     group_effect,
                                     interaction_effect),
                          epoch = rep(epoch_grid, 4),
                          type = factor(rep(1:4, each = length(epoch_grid))))
type.labs <- c("Baseline", "Age effect", "Group effect", "Interaction effect")
names(type.labs) <- c(1, 2, 3, 4)
var_effects %>%
  ggplot(aes(x = epoch)) +
  geom_line(aes(y = effect, linetype = type)) +
  labs(x = "30 second epoch" , y = substitute(paste(italic('g(t, x)'))),
       linetype = "Covariate effect") +
  scale_linetype_manual(labels = c("Baseline", "Age effect", "Hypertension effect", "Interaction effect"), values = c(1,3,5,6)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

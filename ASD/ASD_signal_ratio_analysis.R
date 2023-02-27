#########################################
### ASD case study with signal ratio ####
#########################################
library(latex2exp)
library(gridExtra)
library(bfcr)
library(mgcv)
library(tidyverse)
library(fdapace)
library(plotly)
library(rlang)
setwd()
source("./Peak_Alpha_Data_Transfer.R")

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
signal_ratio_grid <- pa.dat %>%
  group_by(ID) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(Signal_ratio) %>%
  pull()
group_grid <- pa.dat %>%
  group_by(ID) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  mutate(group = group - 1) %>%
  select(group) %>%
  pull()
num_times <- length(freq_grid)
response <- pa.dat %>%
  filter(reg == 15) %>%
  select(y) %>%
  pull() %>%
  matrix(nrow = num_times, ncol = num_subjects) %>%
  t()

# Set up p-spline for frequency dimension, using 12 basis functions
freq_basis_obj <- smoothCon(s(freq_grid, k = 12, bs = "ps", m = 2),
                            data.frame(freq_grid),
                            absorb.cons = FALSE)
freq_basis <- freq_basis_obj[[1]]$X
freq_penalty <- freq_basis_obj[[1]]$S[[1]] * freq_basis_obj[[1]]$S.scale

{
  # Set up p-splines for age dimension by diagnostic group, using 5
  # basis functions (6 including an intercept)
  age_basis_obj <- smoothCon(s(age_grid, k = 6, bs = "ps",
                               by = as.factor(group_grid)), data =
                               data.frame(age_grid, group_grid),
                             absorb.cons = TRUE)
  age_penalty <- age_basis_obj[[1]]$S[[1]] * age_basis_obj[[1]]$S.scale
  signal_ratio_basis_obj <- smoothCon(s(signal_ratio_grid, k = 6, bs = "ps"),
                                      data = data.frame(signal_ratio_grid),
                                      absorb.cons = TRUE)
  signal_ratio_penatly <- signal_ratio_basis_obj[[1]]$S[[1]] * signal_ratio_basis_obj[[1]]$S.scale
  model_penalties <- tensor.prod.penalties(list(age_penalty, freq_penalty))
  model_penalties1 <- tensor.prod.penalties(list(signal_ratio_penatly, freq_penalty))
  design_mean <- cbind(1, group_grid,
                       age_basis_obj[[1]]$X + age_basis_obj[[2]]$X,
                       age_basis_obj[[2]]$X, signal_ratio_basis_obj[[1]]$X)
  design_var <- cbind(1, group_grid,
                      age_basis_obj[[1]]$X + age_basis_obj[[2]]$X,
                      age_basis_obj[[2]]$X, signal_ratio_basis_obj[[1]]$X)
  mean_penalty <- list(freq_penalty, freq_penalty,
                       model_penalties[[1]], model_penalties[[2]],
                       model_penalties[[1]], model_penalties[[2]],
                       model_penalties1[[1]], model_penalties1[[2]])
  var_penalty <- list(freq_penalty, freq_penalty,
                      model_penalties[[1]], model_penalties[[2]],
                      model_penalties[[1]], model_penalties[[2]],
                      model_penalties1[[1]], model_penalties1[[2]])
  
  # Set up smoothing parameters
  var_indices <- c(1, 2, 3, 3, 4, 4, 5, 5)
  mean_indices <- c(1, 2, 3, 3, 4, 4, 5, 5)
  evaluate_basis <- function(obj, obj2, age, signal_ratio, group) {
    spline_part1 <-
      PredictMat(obj[[1]], data.frame(age_grid = age, group_grid = group)) +
      PredictMat(obj[[2]], data.frame(age_grid = age, group_grid = group))
    
    spline_part2 <-
      PredictMat(obj[[2]], data.frame(age_grid = age, group_grid = group))
    spline_part3 <- 
      PredictMat(obj2[[1]], data.frame(signal_ratio_grid = signal_ratio))
    c(1, group, spline_part1, spline_part2, spline_part3)
  }
}

k <- 8
iter <- 25000
burnin <- 12500
thin <- 20

mcmc_results <- run_mcmc(response, design_mean,
                         design_var, freq_basis,
                         freq_grid,
                         mean_penalty, var_penalty,
                         mean_indices, var_indices,
                         k, iter, burnin, thin = thin,
                         var = "pooled")
file_name <- paste0("./signal_ratio_group_mcmc_results_final.rds")
saveRDS(mcmc_results, file = file_name)
mcmc_results <- readRDS(file_name)
freq_grid <- seq(from = 6, to = 14, by = .05)
num_times <- length(freq_grid)
freq_basis_obj <- smoothCon(s(freq_grid, k = 12, bs = "ps", m = 2),
                            data.frame(freq_grid),
                            absorb.cons = FALSE)
mcmc_results$data$time <- freq_grid
mcmc_results$data$basis <- freq_basis_obj[[1]]$X

# Plot age/group adjusted alpha power spectral density means
new_age_grid <- c(40, 120)
new_signal_ratio_grid <- c(0, 0.5, 1)
mean_tibble <- tibble(mean = numeric(num_times * 6 * 2),
                      lower = numeric(num_times * 6 * 2),
                      upper = numeric(num_times * 6 * 2),
                      age = numeric(num_times * 6 * 2),
                      signal_ratio = numeric(num_times * 6 * 2),
                      group = factor(rep(rep(c(0, 1), each = num_times),
                                         6)),
                      freq = numeric(2 * num_times * 6))
first_index <- 1
for(j in c(0, 1)) {
  for(i in seq_along(new_age_grid)) {
    for(l in 1:length(new_signal_ratio_grid)){
      xi <- evaluate_basis(age_basis_obj, signal_ratio_basis_obj ,new_age_grid[i], new_signal_ratio_grid[l], j)
      mean_bands <- get_posterior_means(mcmc_results, xi)
      last_index <- first_index + num_times - 1
      mean_tibble$age[first_index:last_index] <- new_age_grid[i]
      mean_tibble$signal_ratio[first_index:last_index] <- new_signal_ratio_grid[l]
      mean_tibble$mean[first_index:last_index] <- mean_bands$mean
      mean_tibble$group[first_index:last_index] <- j
      mean_tibble$lower[first_index:last_index] <- mean_bands$lower
      mean_tibble$upper[first_index:last_index] <- mean_bands$upper
      mean_tibble$freq[first_index:last_index] <- freq_grid
      first_index <- last_index + 1
    }
  }
}

mean_tibble %>%
  ggplot(aes(x = freq)) +
  geom_line(aes(y = mean, lty = group), size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = group),
              alpha = .3, show.legend = FALSE) +
  facet_wrap(~ signal_ratio + age, nrow = 2, labeller = labeller(
    signal_ratio = ~ paste("Signal Ratio: ", .),
    age = ~ paste("Age: ", ., "months"),
    .multi_line = T
  )) +
  labs(x = expression(omega), y = expression(paste(mu,"(", omega, ", ", x, ")"))) +
  scale_linetype_manual(values = c(1,4), name = "Group", labels = c("TD", "ASD")) +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 3))    +
  theme(axis.title.y = element_text(margin = margin(r = 10, b = 0, l = 0)))


mean_tibble$age <- as.factor(mean_tibble$age)
mean_tibble <- mean_tibble[mean_tibble$group == 1,]
mean_tibble %>%
  ggplot(aes(x = freq)) +
  geom_line(aes(y = mean, lty = age), size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = age),
              alpha = .3, show.legend = FALSE) +
  facet_wrap(~ signal_ratio, nrow = 1) +
  labs(x = expression(omega), y = expression(paste(mu,"(", omega, ", ", x, ")"))) +
  scale_linetype_manual(values = c(1,4), name = "Age", labels = c("40 months", "120 months")) +
  scale_fill_manual(values = c("black", "black")) +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 3))    +
  theme(axis.title.y = element_text(margin = margin(r = 10, b = 0, l = 0)))


# Age-group adjusted eigenfunction
N <- 50
lower_signal_ratio <- 0
upper_signal_ratio <- 1
new_signal_ratio_grid <- seq(from = lower_signal_ratio, to = upper_signal_ratio, length.out = N)
eigen_surface_td <- matrix(0, N, num_times)
covariance_td <- array(0, c(num_times, num_times, N))
magnitude_td <- matrix(0, N, 3)
for (i in 1:N) {
  zi <- evaluate_basis(age_basis_obj, signal_ratio_basis_obj , 120, new_signal_ratio_grid[i], 0)
  covariance_td[,,i] <- matrix(get_posterior_covariance(mcmc_results, zi)$mean, nrow = num_times, ncol = num_times)
  eigen_surface_td[i,] <- eigen(covariance_td[,,i])$vectors[,1]
  if (eigen_surface_td[i, 1] < 0) {
    eigen_surface_td[i,] <- -1 * eigen_surface_td[i,]
  }
  magnitude_td[i,] <- eigen(covariance_td[,,i])$values[1:3]
}
eigen_surface_asd <- matrix(0, N, num_times)
covariance_asd <- array(0, c(num_times, num_times, N))
magnitude_asd <- matrix(0, N, 3)
for (i in 1:N) {
  zi <- evaluate_basis(age_basis_obj, signal_ratio_basis_obj , 120, new_signal_ratio_grid[i], 1)
  covariance_asd[,,i] <- matrix(get_posterior_covariance(mcmc_results, zi)$mean, nrow = num_times, ncol = num_times)
  eigen_surface_asd[i,] <- eigen(covariance_asd[,,i])$vectors[,1]
  if (eigen_surface_asd[i, 1] < 0) {
    eigen_surface_asd[i,] <- -1 * eigen_surface_asd[i,]
  }
  magnitude_asd[i,] <- eigen(covariance_asd[,,i])$values[1:3]
}
x <- freq_grid
y <- new_signal_ratio_grid
data <- expand.grid(X=x, Y=y, G = factor(c(0, 1)))
data$Z <- c(c(t(eigen_surface_td)), c(t(eigen_surface_asd)))
group.labs <- c("TD", "ASD")
names(group.labs) <- c(0, 1)
ggplot(data, aes(X, Y, fill = Z)) +
  geom_tile() +
  labs(x = expression(paste(omega)), y = "Signal Ratio") +
  facet_wrap(~ G, labeller = labeller(G = group.labs)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines")) +
  scale_fill_gradient(name = "Eigenfunction 1", low = "blue", high = "red")

x <- freq_grid
y <- new_signal_ratio_grid
var_td <- eigen_surface_td
for(i in 1:nrow(var_td)){
  var_td[i,] <- diag(covariance_td[,,i])
}
var_asd <- eigen_surface_asd
for(i in 1:nrow(var_asd)){
  var_asd[i,] <- diag(covariance_asd[,,i])
}
data <- expand.grid(X=x, Y=y, G = factor(c(0, 1)))
data$Z <- c(c(t(var_td)), c(t(var_asd)))
group.labs <- c("TD", "ASD")
names(group.labs) <- c(0, 1)
ggplot(data, aes(X, Y, fill = Z)) +
  geom_tile() +
  labs(x = expression(paste(omega)), y = "Signal Ratio") +
  facet_wrap(~ G, labeller = labeller(G = group.labs)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines")) +
  scale_fill_gradient(name = "Variance", low = "blue", high = "red")



x <- freq_grid
y <- freq_grid
data <- expand.grid(X=x, Y=y, G = factor(c(0, 1, 2, 3, 4, 5)))
data$Z <- c(c(t(covariance_td[,,50])), c(t(covariance_td[,,25])), c(t(covariance_td[,,1])),
            c(t(covariance_asd[,,50])), c(t(covariance_asd[,,25])), c(t(covariance_asd[,,1])))
group.labs <- c("TD (Signal Ratio = 1)", "TD (Signal Ratio = 0.5)", "TD (Signal Ratio = 0)",
                "ASD (Signal Ratio = 1)", "ASD (Signal Ratio = 0.5)", "ASD (Signal Ratio = 0)")
names(group.labs) <- c(0, 1, 2, 3, 4, 5)
ggplot(data, aes(X, Y, fill = Z)) +
  geom_tile() +
  labs(x = expression(paste(omega)), y = expression(paste(omega))) +
  facet_wrap(~ G, labeller = labeller(G = group.labs), nrow = 2) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines")) +
  scale_fill_gradient(name = "Covariance", low = "blue", high = "red")

svd_ph <- eigen(covariance_td[,,50])
len <- length(svd_ph$vectors[,1])
if(svd_ph$vectors[1,1] < 0){
  svd_ph$vectors[,1] <- -1 * svd_ph$vectors[,1]
}
if(svd_ph$vectors[1,2] < 0){
  svd_ph$vectors[,2] <- -1 * svd_ph$vectors[,2]
}
if(svd_ph$vectors[1,3] < 0){
  svd_ph$vectors[,3] <- -1 * svd_ph$vectors[,3]
}
eigen_ph <- matrix(0, 3*len, 3)
eigen_ph[1:len,1] <- svd_ph$vectors[,1]
eigen_ph[1:len,2] <- "eigenfunction 1"
eigen_ph[1:len,3] <- freq_grid
eigen_ph[(len + 1):(2*len),1] <- svd_ph$vectors[,2]
eigen_ph[(len + 1):(2*len),2] <- "eigenfunction 2"
eigen_ph[(len + 1):(2*len),3] <- freq_grid
eigen_ph[(2*len + 1):(3*len),1] <- svd_ph$vectors[,3]
eigen_ph[(2*len + 1):(3*len),2] <- "eigenfunction 3"
eigen_ph[(2*len + 1):(3*len),3] <- freq_grid
eigen_ph <- as.data.frame(eigen_ph)
colnames(eigen_ph) <- c("y", "eigenfunction", "frequency")
eigen_ph$eigenfunction <- as.factor(eigen_ph$eigenfunction)
eigen_ph$y <- as.numeric(eigen_ph$y)
eigen_ph$frequency <- as.numeric(eigen_ph$frequency)
p1 <- ggplot(eigen_ph, aes(x= `frequency`, y=y, group = `eigenfunction`, colour= `eigenfunction`)) +  ggtitle("TD (Signal Ratio = 1)") +
  geom_line() +  theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                           ,plot.title = element_text(hjust = 0.5))

svd_ph <- eigen(covariance_td[,,25])
len <- length(svd_ph$vectors[,1])
if(svd_ph$vectors[1,1] < 0){
  svd_ph$vectors[,1] <- -1 * svd_ph$vectors[,1]
}
if(svd_ph$vectors[1,2] < 0){
  svd_ph$vectors[,2] <- -1 * svd_ph$vectors[,2]
}
if(svd_ph$vectors[1,3] < 0){
  svd_ph$vectors[,3] <- -1 * svd_ph$vectors[,3]
}
eigen_ph <- matrix(0, 3*len, 3)
eigen_ph[1:len,1] <- svd_ph$vectors[,1]
eigen_ph[1:len,2] <- "eigenfunction 1"
eigen_ph[1:len,3] <- freq_grid
eigen_ph[(len + 1):(2*len),1] <- svd_ph$vectors[,2]
eigen_ph[(len + 1):(2*len),2] <- "eigenfunction 2"
eigen_ph[(len + 1):(2*len),3] <- freq_grid
eigen_ph[(2*len + 1):(3*len),1] <- svd_ph$vectors[,3]
eigen_ph[(2*len + 1):(3*len),2] <- "eigenfunction 3"
eigen_ph[(2*len + 1):(3*len),3] <- freq_grid
eigen_ph <- as.data.frame(eigen_ph)
colnames(eigen_ph) <- c("y", "eigenfunction", "frequency")
eigen_ph$eigenfunction <- as.factor(eigen_ph$eigenfunction)
eigen_ph$y <- as.numeric(eigen_ph$y)
eigen_ph$frequency <- as.numeric(eigen_ph$frequency)
p2 <- ggplot(eigen_ph, aes(x= `frequency`, y=y, group = `eigenfunction`, colour= `eigenfunction`)) +  ggtitle("TD (Signal Ratio = 0.5)") +
  geom_line() +  theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                         ,plot.title = element_text(hjust = 0.5))

svd_ph <- eigen(covariance_td[,,1])
len <- length(svd_ph$vectors[,1])
if(svd_ph$vectors[1,1] < 0){
  svd_ph$vectors[,1] <- -1 * svd_ph$vectors[,1]
}
if(svd_ph$vectors[1,2] < 0){
  svd_ph$vectors[,2] <- -1 * svd_ph$vectors[,2]
}
if(svd_ph$vectors[1,3] < 0){
  svd_ph$vectors[,3] <- -1 * svd_ph$vectors[,3]
}
eigen_ph <- matrix(0, 3*len, 3)
eigen_ph[1:len,1] <- svd_ph$vectors[,1]
eigen_ph[1:len,2] <- "eigenfunction 1"
eigen_ph[1:len,3] <- freq_grid
eigen_ph[(len + 1):(2*len),1] <- svd_ph$vectors[,2]
eigen_ph[(len + 1):(2*len),2] <- "eigenfunction 2"
eigen_ph[(len + 1):(2*len),3] <- freq_grid
eigen_ph[(2*len + 1):(3*len),1] <- svd_ph$vectors[,3]
eigen_ph[(2*len + 1):(3*len),2] <- "eigenfunction 3"
eigen_ph[(2*len + 1):(3*len),3] <- freq_grid
eigen_ph <- as.data.frame(eigen_ph)
colnames(eigen_ph) <- c("y", "eigenfunction", "frequency")
eigen_ph$eigenfunction <- as.factor(eigen_ph$eigenfunction)
eigen_ph$y <- as.numeric(eigen_ph$y)
eigen_ph$frequency <- as.numeric(eigen_ph$frequency)
p3 <- ggplot(eigen_ph, aes(x= `frequency`, y=y, group = `eigenfunction`, colour= `eigenfunction`)) +  ggtitle("TD (Signal Ratio = 0)") +
  geom_line() +  theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                         ,plot.title = element_text(hjust = 0.5))

svd_ph <- eigen(covariance_asd[,,50])
len <- length(svd_ph$vectors[,1])
if(svd_ph$vectors[1,1] < 0){
  svd_ph$vectors[,1] <- -1 * svd_ph$vectors[,1]
}
if(svd_ph$vectors[1,2] < 0){
  svd_ph$vectors[,2] <- -1 * svd_ph$vectors[,2]
}
if(svd_ph$vectors[1,3] < 0){
  svd_ph$vectors[,3] <- -1 * svd_ph$vectors[,3]
}
eigen_ph <- matrix(0, 3*len, 3)
eigen_ph[1:len,1] <- svd_ph$vectors[,1]
eigen_ph[1:len,2] <- "eigenfunction 1"
eigen_ph[1:len,3] <- freq_grid
eigen_ph[(len + 1):(2*len),1] <- svd_ph$vectors[,2]
eigen_ph[(len + 1):(2*len),2] <- "eigenfunction 2"
eigen_ph[(len + 1):(2*len),3] <- freq_grid
eigen_ph[(2*len + 1):(3*len),1] <- svd_ph$vectors[,3]
eigen_ph[(2*len + 1):(3*len),2] <- "eigenfunction 3"
eigen_ph[(2*len + 1):(3*len),3] <- freq_grid
eigen_ph <- as.data.frame(eigen_ph)
colnames(eigen_ph) <- c("y", "eigenfunction", "frequency")
eigen_ph$eigenfunction <- as.factor(eigen_ph$eigenfunction)
eigen_ph$y <- as.numeric(eigen_ph$y)
eigen_ph$frequency <- as.numeric(eigen_ph$frequency)
p4 <- ggplot(eigen_ph, aes(x= `frequency`, y=y, group = `eigenfunction`, colour= `eigenfunction`)) +  ggtitle("ASD (Signal Ratio = 1)") +
  geom_line() +  theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                         ,plot.title = element_text(hjust = 0.5))

svd_ph <- eigen(covariance_asd[,,25])
len <- length(svd_ph$vectors[,1])
if(svd_ph$vectors[1,1] < 0){
  svd_ph$vectors[,1] <- -1 * svd_ph$vectors[,1]
}
if(svd_ph$vectors[1,2] < 0){
  svd_ph$vectors[,2] <- -1 * svd_ph$vectors[,2]
}
if(svd_ph$vectors[1,3] < 0){
  svd_ph$vectors[,3] <- -1 * svd_ph$vectors[,3]
}
eigen_ph <- matrix(0, 3*len, 3)
eigen_ph[1:len,1] <- svd_ph$vectors[,1]
eigen_ph[1:len,2] <- "eigenfunction 1"
eigen_ph[1:len,3] <- freq_grid
eigen_ph[(len + 1):(2*len),1] <- svd_ph$vectors[,2]
eigen_ph[(len + 1):(2*len),2] <- "eigenfunction 2"
eigen_ph[(len + 1):(2*len),3] <- freq_grid
eigen_ph[(2*len + 1):(3*len),1] <- svd_ph$vectors[,3]
eigen_ph[(2*len + 1):(3*len),2] <- "eigenfunction 3"
eigen_ph[(2*len + 1):(3*len),3] <- freq_grid
eigen_ph <- as.data.frame(eigen_ph)
colnames(eigen_ph) <- c("y", "eigenfunction", "frequency")
eigen_ph$eigenfunction <- as.factor(eigen_ph$eigenfunction)
eigen_ph$y <- as.numeric(eigen_ph$y)
eigen_ph$frequency <- as.numeric(eigen_ph$frequency)
p5 <- ggplot(eigen_ph, aes(x= `frequency`, y=y, group = `eigenfunction`, colour= `eigenfunction`)) +  ggtitle("ASD (Signal Ratio = 0.5)") +
  geom_line() +  theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                         ,plot.title = element_text(hjust = 0.5))

svd_ph <- eigen(covariance_asd[,,1])
len <- length(svd_ph$vectors[,1])
if(svd_ph$vectors[1,1] < 0){
  svd_ph$vectors[,1] <- -1 * svd_ph$vectors[,1]
}
if(svd_ph$vectors[1,2] < 0){
  svd_ph$vectors[,2] <- -1 * svd_ph$vectors[,2]
}
if(svd_ph$vectors[1,3] < 0){
  svd_ph$vectors[,3] <- -1 * svd_ph$vectors[,3]
}
eigen_ph <- matrix(0, 3*len, 3)
eigen_ph[1:len,1] <- svd_ph$vectors[,1]
eigen_ph[1:len,2] <- "eigenfunction 1"
eigen_ph[1:len,3] <- freq_grid
eigen_ph[(len + 1):(2*len),1] <- svd_ph$vectors[,2]
eigen_ph[(len + 1):(2*len),2] <- "eigenfunction 2"
eigen_ph[(len + 1):(2*len),3] <- freq_grid
eigen_ph[(2*len + 1):(3*len),1] <- svd_ph$vectors[,3]
eigen_ph[(2*len + 1):(3*len),2] <- "eigenfunction 3"
eigen_ph[(2*len + 1):(3*len),3] <- freq_grid
eigen_ph <- as.data.frame(eigen_ph)
colnames(eigen_ph) <- c("y", "eigenfunction", "frequency")
eigen_ph$eigenfunction <- as.factor(eigen_ph$eigenfunction)
eigen_ph$y <- as.numeric(eigen_ph$y)
eigen_ph$frequency <- as.numeric(eigen_ph$frequency)
p6 <- ggplot(eigen_ph, aes(x= `frequency`, y=y, group = `eigenfunction`, colour= `eigenfunction`)) +  ggtitle("ASD (Signal Ratio = 0)") +
  geom_line() +  theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                         ,plot.title = element_text(hjust = 0.5))

grid.arrange(p1, p2, p3, p4, p5, p6, layout_matrix = rbind(c(1, 1, 2, 2, 3, 3),
                                                           c(4, 4, 5, 5, 6, 6)))


# Low dimensional g effects
new_age_grid <- c(seq(40, 120, by = 5))
new_signal_ratio_grid <- c(seq(0, 1, by = 0.05))

freq_basis <- freq_basis_obj[[1]]$X
intercept_effect <- numeric(num_times)
age_effect <- numeric(num_times)
group_effect <- numeric(num_times)
interaction_effect <- numeric(num_times)
signal_ratio_effect <- numeric(num_times)

for (iter in 12501:25000) {
  for (kp in 1:k) {
    intercept_effect <- intercept_effect + (freq_basis %*% mcmc_results$samples$lambda[[iter]][,1,kp])^2
    group_effect <- group_effect + c((freq_basis %*% mcmc_results$samples$lambda[[iter]][,2,kp])^2)
  }
}

for (i in seq_along(new_age_grid)) {
  age_zi <- evaluate_basis(age_basis_obj, signal_ratio_basis_obj, new_age_grid[i], 0, 1)
  interaction_zi <- age_zi[8:12]
  age_zi <- age_zi[3:7]
  for (iter in 12501:25000) {
    for (kp in 1:k) {
      age_effect <- age_effect +
        5 * c((freq_basis %*% mcmc_results$samples$lambda[[iter]][,3:7,kp] %*% age_zi)^2)
      interaction_effect <- interaction_effect +
        c((freq_basis %*% mcmc_results$samples$lambda[[iter]][,8:12,kp] %*% interaction_zi)^2)
    }
  }
}

for(i in seq_along(new_signal_ratio_grid)){
  signal_ratio_zi <- evaluate_basis(age_basis_obj, signal_ratio_basis_obj, 40, new_signal_ratio_grid, 0)
  signal_ratio_zi <- signal_ratio_zi[13:17]
  for (iter in 12501:25000) {
    for (kp in 1:k) {
      signal_ratio_effect <- signal_ratio_effect + 0.05*
        c((freq_basis %*% mcmc_results$samples$lambda[[iter]][,13:17,kp] %*% signal_ratio_zi)^2)
    }
  }
}
intercept_effect <- intercept_effect / 12500
age_effect <- age_effect / (12500 * 80)
group_effect <- group_effect / 12500
interaction_effect <- interaction_effect / (12500 * 80)
signal_ratio_effect <- signal_ratio_effect / (12500)

var_effects <- data.frame(effect = c(intercept_effect,
                                     age_effect,
                                     group_effect,
                                     interaction_effect,
                                     signal_ratio_effect),
                          freq = rep(freq_grid, 5),
                          type = factor(rep(1:5, each = length(freq_grid))))
type.labs <- c("Baseline", "Age effect", "Group effect", "Interaction effect", "Signal Ratio effect")
names(type.labs) <- c(1, 2, 3, 4, 5)
var_effects %>%
  ggplot(aes(x = freq),) +
  geom_line(aes(y = effect, linetype = type, col = type)) +
  labs(x = expression(paste(omega)), y = TeX("$g_r(t)$"),
       linetype = "Covariate effect") +
  scale_linetype_manual(
    labels = c("Baseline", "Age effect", "Group effect", "Interaction effect", "Signal Ratio effect"),
    values = c(1,3,5,6,2)) +
  scale_colour_discrete(name = "Covariate effect", labels = c("Baseline", "Age effect", "Group effect", "Interaction effect", "Signal Ratio effect")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))



library(MASS)
library(plot3D)
library(dlnm)
library(BayesianConditionalFPCA)
library(loo)
library(tidyverse)
setwd("E:\\PAF_visualization")
source('Peak Alpha Data Transfer.r')
pa.dat.tibble.reg %>%
  ggplot(aes(x = func, y = meanx, group = ID)) +
  facet_wrap(~group, labeller = labeller(group = group.labs)) +
  geom_line() +
  xlab("Frequency") +
  ylab("Spectral power") +
  theme_bw()

Basis <- ps(t_paf, df = 16, intercept = TRUE)

# X_paf columns
# Column 1: Intercept
# Column 2: Group. 0 = TD, 1 = ASD
# Column 3: Age
# Column 4: Interaction of group and age

### MCMC sanity check ###
K <- 3
mcmc_results <- run_mcmc_Morris(Y_paf, t_paf, X_paf, X_paf, Basis, K, iter = 5000, burnin = 10000, nchains = 1, thin = 5, loglik = 1)

### Visualization ###
sub <- 10
posterior_intervals <- get_posterior_predictive_bands2(mcmc_results, c(.025, .5, .975))
colnames(posterior_intervals) <- c("ID", "Frequency", "Y", "Lower_P", "Median_P", "Upper_P", "Lower_M", "Median_M", "Upper_M")
posterior_intervals <- as_tibble(posterior_intervals)
posterior_intervals %>%
  filter(ID == sub) %>%
  ggplot(aes(x = Frequency, y = Y)) +
  geom_point(na.rm = TRUE) +
  geom_ribbon(aes(ymin = Lower_P, ymax = Upper_P), alpha = 0.3) +
  theme_bw()
posterior_intervals %>%
  filter(ID == sub) %>%
  ggplot(aes(x = Frequency, y = Y)) +
  geom_point(na.rm = TRUE) +
  geom_line(aes(y = Median_M)) +
  theme_bw()
posterior_intervals %>%
  filter(ID == sub) %>%
  ggplot(aes(x = Frequency, y = Y)) +
  geom_point(na.rm = TRUE) +
  geom_ribbon(aes(ymin = Lower_M, ymax = Upper_M), alpha = 0.3) +
  theme_bw()
posterior_intervals %>%
  mutate(coverage_flag = ifelse(Y > Lower_P & Y < Upper_P, 1, 0)) %>%
  group_by(ID) %>%
  summarize(coverage_summary = sum(coverage_flag, na.rm = TRUE) / n())
posterior_intervals %>%
  mutate(coverage_flag = ifelse(Y > Lower_P & Y < Upper_P, 1, 0)) %>%
  summarize(coverage_summary = sum(coverage_flag, na.rm = TRUE) / n())

### Goodness of fit ###
mcmc_omnibus_fit <- get_omnibus_fit2(mcmc_results)
plot(mcmc_omnibus_fit$statistic_obs, mcmc_omnibus_fit$statistic_rep, xlab = "Chi2 observed", ylab = "Chi2 repitition")
abline(a = 0, b = 1)
sum(mcmc_omnibus_fit$statistic_rep > mcmc_omnibus_fit$statistic_obs) / 5000

### LOO ###
r_eff <- relative_eff(exp(mcmc_results$log_lik), cores = 4)
loo_1 <- loo(mcmc_results$log_lik, r_eff = r_eff, cores = 4)
print(loo_1)

### Posterior mean bands ###
alpha <- .05

# Younger TD
Group <- 0
Month <- 30
poly_coef <- predict(poly_age, Month)
xi <- c(1, Group, poly_coef, Group * poly_coef)
coef_bands <- get_posterior_means(mcmc_results, xi, alpha)
coef_bands <- cbind(t_paf, coef_bands)
colnames(coef_bands) <- c("Frequency", "Lower", "Mean", "Upper")
coef_bands <- as_tibble(coef_bands)
p1 <- coef_bands %>%
  ggplot(aes(x = Frequency)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, fill = "blue") +
  geom_line(aes(y = Mean)) +
  xlab("Frequency") +
  ylab("Power") +
  ggtitle("Young TD") + 
  theme_bw() 

# Older TD
Group <- 0
Month <- 120
poly_coef <- predict(poly_age, Month)
xi <- c(1, Group, poly_coef, Group * poly_coef)
coef_bands <- get_posterior_means(mcmc_results, xi, alpha)
coef_bands <- cbind(t_paf, coef_bands)
colnames(coef_bands) <- c("Frequency", "Lower", "Mean", "Upper")
coef_bands <- as_tibble(coef_bands)
p2 <- coef_bands %>%
  ggplot(aes(x = Frequency)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, fill = "blue") +
  geom_line(aes(y = Mean)) +
  xlab("Frequency") +
  ylab("Power") +
  ggtitle("Old TD") + 
  theme_bw() 
gridExtra::grid.arrange(p1, p2, nrow = 1)

# Younger ASD
Group <- 1
Month <- 30
poly_coef <- predict(poly_age, Month)
xi <- c(1, Group, poly_coef, Group * poly_coef)
coef_bands <- get_posterior_means(mcmc_results, xi, alpha)
coef_bands <- cbind(t_paf, coef_bands)
colnames(coef_bands) <- c("Frequency", "Lower", "Mean", "Upper")
coef_bands <- as_tibble(coef_bands)
p3 <- coef_bands %>%
  ggplot(aes(x = Frequency)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, fill = "blue") +
  geom_line(aes(y = Mean)) +
  xlab("Frequency") +
  ylab("Response") +
  ggtitle("Young ASD") + 
  theme_bw() 

# Older ASD
Group <- 1
Month <- 120
poly_coef <- predict(poly_age, Month)
xi <- c(1, Group, poly_coef, Group * poly_coef)
coef_bands <- get_posterior_means(mcmc_results, xi, alpha)
coef_bands <- cbind(t_paf, coef_bands)
colnames(coef_bands) <- c("Frequency", "Lower", "Mean", "Upper")
coef_bands <- as_tibble(coef_bands)
p4 <- coef_bands %>%
  ggplot(aes(x = Frequency)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, fill = "blue") +
  geom_line(aes(y = Mean)) +
  xlab("Frequency") +
  ylab("Response") +
  ggtitle("Old ASD") + 
  theme_bw() 
gridExtra::grid.arrange(p3, p4, nrow = 1)

### Some covariance visualization ###
Group <- 0
Month <- 60
poly_coef <- predict(poly_age, Month)
zi <- c(1, Group, poly_coef, Group * poly_coef)
evals <- 3
eigen_bands <- get_posterior_eigen2(mcmc_results, evals, zi, alpha)
eig_names <- c()
eig_labs <- c()
for(k in 1:evals){
  eig_names <- c(eig_names, paste("Eigenfunction", k))
  eig_labs <- c(eig_labs, paste("Eigenfunction", k, " ", round(eigen_bands$eigenval_pve_intervals[1,k], 2), "-", round(eigen_bands$eigenval_pve_intervals[3,k], 2)))
}
names(eig_labs) <- eig_names
eigen_bands_tibble <- tibble(Frequency = rep(t_paf, evals),
                             number = factor(rep(eig_names, each = length(t_paf))),
                             lower = c(eigen_bands$lower),
                             mean = c(eigen_bands$mean),
                             upper = c(eigen_bands$upper),
                             val_lower = rep(eigen_bands$eigenval_intervals[1,], each = length(t_paf)),
                             val_median = rep(eigen_bands$eigenval_intervals[2,], each = length(t_paf)),
                             val_upper = rep(eigen_bands$eigenval_intervals[3,], each = length(t_paf)))

eigen_bands_tibble %>%
  ggplot(aes(x = Frequency)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = number), alpha = 0.5) +
  geom_line(aes(y = mean)) +
  facet_wrap(number ~., labeller = labeller(number = eig_labs)) +
  ylab("") +
  theme_bw() +
  theme(legend.position="none")
eigen_bands$magnitude
aX <- list(title = "Frequency")
aY <- list(title = "Frequency")
aZ <- list(title = "Response")
plotly::plot_ly(x = t_paf, y = t_paf, z = eigen_bands$surface, type = "surface") %>%
  plotly::layout(scene = list(xaxis = aX, yaxis = aY, zaxis = aZ, dragmode = "turntable"))
# plot(eigen_bands$raw_magnitude, type = "l")
persp3D(t_paf, t_paf, eigen_bands$surface)


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

Basis <- ps(t_paf, df = 12, intercept = TRUE)
# X_paf <- cbind(X_paf, X_paf[,3]^2, X_paf[,4]^2, X_paf[,3]^3, X_paf[,4]^3)
# ps_age <- ps(X_paf[,3], df = 4, intercept = FALSE)
# X_paf <- cbind(X_paf[,1], X_paf[,2], ps_age, X_paf[,2] * ps_age)
# X_paf columns
# Column 1: Intercept
# Column 2: Group. 0 = TD, 1 = ASD
# Column 3: Age
# Column 4: Interaction of group and age

### MCMC ###
K <- 4
mcmc_results <- run_mcmc_Morris(Y_paf, t_paf, X_paf, X_paf[,1:4], Basis, K, iter = 10000, burnin = 2500, nchains = 1, thin = 5, loglik = 1)

### Visualization ###
sub <- 2
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
Month <- .4
# poly_coef <- predict(poly_age, Month)
# xi <- c(1, Group, poly_coef, Group * poly_coef)
xi <- c(1, Group, Month, Group * Month)
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
Month <- 1.5
# poly_coef <- predict(poly_age, Month)
# xi <- c(1, Group, poly_coef, Group * poly_coef)
xi <- c(1, Group, Month, Group * Month)
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
Month <- -2
# poly_coef <- predict(poly_age, Month)
# xi <- c(1, Group, poly_coef, Group * poly_coef)
xi <- c(1, Group, Month, Group * Month)
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
Month <- 1.5
# poly_coef <- predict(poly_age, Month)
# xi <- c(1, Group, poly_coef, Group * poly_coef)
xi <- c(1, Group, Month, Group * Month)
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


plot(Basis %*%mcmc_results$Beta[[1]][,,5000]%*%xi, type = "l", ylim = c(0,.6))
for(i in 1:5000){
  lines(Basis %*%mcmc_results$Beta[[1]][,,i]%*%xi)
}
coef_bands <- get_posterior_coefs(mcmc_results, .05)

coef_bands <- tibble(Frequency = rep(t_paf, times = dim(X_paf)[2]),
                     Covariate = rep(c("Intercept", "Group", "Age", "Group x Age"), each = dim(Basis)[1]),
                     Lower = c(coef_bands$lower),
                     Mean = c(coef_bands$mean),
                     Upper = c(coef_bands$upper))
coef_bands %>%
  ggplot(aes(x = Frequency)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Covariate), alpha = 0.5) +
  geom_line(aes(y = Mean)) +
  geom_hline(yintercept = 0) +
  facet_wrap(Covariate ~., scales = "free") +
  ylab(expression(paste("Normalized power density   (", mu, "V"^2, " / Hz)"))) +
  theme_minimal()

Age_seq <- seq(from = -1.7, to = 2.3, by = .1)
Age_rescaled <- mean(demDat$Age) + sd(demDat$Age) * Age_seq
Ages <- c(50, 70, 90, 110)
Age_scaled <- (Ages - mean(demDat$Age))/sd(demDat$Age)

p <- list()
counter <- 1
for(i in 1:4){
  # Older ASD
  Group <- 1
  Month <- Age_scaled[i]
  # poly_coef <- predict(poly_age, Month)
  # xi <- c(1, Group, poly_coef, Group * poly_coef)
  xi <- c(1, Group, Month, Group * Month)
  coef_bands <- get_posterior_means(mcmc_results, xi, alpha)
  coef_bands <- cbind(t_paf, coef_bands)
  colnames(coef_bands) <- c("Frequency", "Lower", "Mean", "Upper")
  coef_bands <- as_tibble(coef_bands)
  p[[counter]] <- coef_bands %>%
    ggplot(aes(x = Frequency)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, fill = "blue") +
    geom_line(aes(y = Mean)) +
    xlab("Frequency") +
    ylab("Mean") +
    ggtitle(paste(Ages[i], "Months, ", "Group = ASD")) + 
    theme_bw() 
  counter <- counter + 1
}

gridExtra::grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]],
                        p[[5]], p[[6]], p[[7]], p[[8]], nrow = 2)


### Some covariance visualization ###
Group <- 0
Month <- 2.3
# poly_coef <- predict(poly_age, Month)
# zi <- c(1, Group, poly_coef, Group * poly_coef)
zi <- c(1, Group, Month, Group * Month, Month^2, Group * Month^2)
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
par(mfrow = c(1,1))
persp3D(t_paf, t_paf, eigen_bands$surface)


Age_seq <- seq(from = -1.7, to = 2.3, by = .1)
Age_rescaled <- mean(demDat$Age) + sd(demDat$Age) * Age_seq
evals <- 4
Cov_mat <- matrix(0, nrow = length(Age_seq) * 2, ncol = 3)
Group <- 0
counter <- 1
for(Month in Age_seq){
  zi <- c(1, Group, Month, Group * Month)
  #zi <- c(1, Group, Month, Group * Month, Month^2, Group * Month^2)
  #zi <- c(1, Group, Month, Group * Month, Month^2, Group * Month^2, Month^3, Month^3 * Group)
  eigen_bands <- get_posterior_eigen2(mcmc_results, evals, zi, alpha)
  Cov_mat[counter, ] <- eigen_bands$magnitude
  counter <- counter + 1
}

Group <- 1
for(Month in Age_seq){
  zi <- c(1, Group, Month, Group * Month)
  #zi <- c(1, Group, Month, Group * Month, Month^2, Group * Month^2)
  #zi <- c(1, Group, Month, Group * Month, Month^2, Group * Month^2, Month^3, Month^3 * Group)
  eigen_bands <- get_posterior_eigen2(mcmc_results, evals, zi, alpha)
  Cov_mat[counter,] <- eigen_bands$magnitude
  counter <- counter + 1
}

Cov_mat <- tibble(Group = rep(c("TD", "ASD"),each = length(Age_seq)), Age = rep(Age_rescaled, times = 2),
                  Lower = Cov_mat[,1], Mean = Cov_mat[,2],
                  Upper = Cov_mat[,3])

Cov_mat %>%
  ggplot(aes(x = Age)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Group), alpha = 0.5) +
  xlab("Age (months)") +
  ylab("Total variance") +
  theme_bw() +
  ggtitle("Heterogeneity within group by age")



Group <- 0
Month <- 36
Month_t <- (Month - mean(demDat$Age)) / sd(demDat$Age)
xi <- c(1, Group, Month_t, Group * Month_t)
coef_bands <- cbind(get_posterior_means(mcmc_results, xi, alpha), rep(Group, length(t_paf)), rep(Month / 12, length(t_paf)), t_paf)

Group <- 0
Month <- 84
Month_t <- (Month - mean(demDat$Age)) / sd(demDat$Age)
xi <- c(1, Group, Month_t, Group * Month_t)
coef_bands <- rbind(coef_bands, cbind(get_posterior_means(mcmc_results, xi, alpha), rep(Group, length(t_paf)), rep(Month / 12, length(t_paf)), t_paf))

Group <- 0
Month <- 132
Month_t <- (Month - mean(demDat$Age)) / sd(demDat$Age)
xi <- c(1, Group, Month_t, Group * Month_t)
coef_bands <- rbind(coef_bands, cbind(get_posterior_means(mcmc_results, xi, alpha), rep(Group, length(t_paf)), rep(Month / 12, length(t_paf)), t_paf))

Group <- 1
Month <- 36
Month_t <- (Month - mean(demDat$Age)) / sd(demDat$Age)
xi <- c(1, Group, Month_t, Group * Month_t)
coef_bands <- rbind(coef_bands, cbind(get_posterior_means(mcmc_results, xi, alpha), rep(Group, length(t_paf)), rep(Month / 12, length(t_paf)), t_paf))

Group <- 1
Month <- 84
Month_t <- (Month - mean(demDat$Age)) / sd(demDat$Age)
xi <- c(1, Group, Month_t, Group * Month_t)
coef_bands <- rbind(coef_bands, cbind(get_posterior_means(mcmc_results, xi, alpha), rep(Group, length(t_paf)), rep(Month / 12, length(t_paf)), t_paf))

Group <- 1
Month <- 132
Month_t <- (Month - mean(demDat$Age)) / sd(demDat$Age)
xi <- c(1, Group, Month_t, Group * Month_t)
coef_bands <- rbind(coef_bands, cbind(get_posterior_means(mcmc_results, xi, alpha), rep(Group, length(t_paf)), rep(Month / 12, length(t_paf)), t_paf))

colnames(coef_bands) <- c("lower", "mean", "upper", "Group", "Age (years)", "Frequency")
coef_bands <- as_tibble(coef_bands)
coef_bands
age.labs <- c("3 years", "7 years", "11 years")
names(age.labs) <- c("3", "7", "11")

group.labs <- c("TD", "ASD")
names(group.labs) <- c("0", "1")
coef_bands %>%
  ggplot(aes(x = Frequency)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  geom_line(aes(y = mean)) +
  facet_grid(Group ~ `Age (years)`,
             labeller = labeller(Group = group.labs, `Age (years)` = age.labs)) +
  labs(y = expression(paste("Normalized power density   (", mu, "V"^2, " / Hz)"))) +
  theme_bw()

Group <- 0
Month <- 36
evals <- 3
Month_t <- (Month - mean(demDat$Age)) / sd(demDat$Age)
zi <- c(1, Group, Month_t, Group * Month_t)
eigen_bands <- cbind(c(get_posterior_eigen2(mcmc_results, evals, zi, alpha)$mean), rep(1:3, each = length(t_paf)), rep(Group, length(t_paf) * 3), rep(Month / 12, length(t_paf) * 3), rep(t_paf, 3))

Group <- 0
Month <- 84
Month_t <- (Month - mean(demDat$Age)) / sd(demDat$Age)
zi <- c(1, Group, Month_t, Group * Month_t)
eigen_bands <- rbind(eigen_bands, cbind(c(get_posterior_eigen2(mcmc_results, evals, zi, alpha)$mean), rep(1:3, each = length(t_paf)), rep(Group, length(t_paf) * 3), rep(Month / 12, length(t_paf) * 3), rep(t_paf, 3)))

Group <- 0
Month <- 132
Month_t <- (Month - mean(demDat$Age)) / sd(demDat$Age)
zi <- c(1, Group, Month_t, Group * Month_t)
eigen_bands <- rbind(eigen_bands, cbind(c(get_posterior_eigen2(mcmc_results, evals, zi, alpha)$mean), rep(1:3, each = length(t_paf)), rep(Group, length(t_paf) * 3), rep(Month / 12, length(t_paf) * 3), rep(t_paf, 3)))

Group <- 1
Month <- 36
Month_t <- (Month - mean(demDat$Age)) / sd(demDat$Age)
zi <- c(1, Group, Month_t, Group * Month_t)
eigen_bands <- rbind(eigen_bands, cbind(c(get_posterior_eigen2(mcmc_results, evals, zi, alpha)$mean), rep(1:3, each = length(t_paf)), rep(Group, length(t_paf) * 3), rep(Month / 12, length(t_paf) * 3), rep(t_paf, 3)))

Group <- 1
Month <- 84
Month_t <- (Month - mean(demDat$Age)) / sd(demDat$Age)
zi <- c(1, Group, Month_t, Group * Month_t)
eigen_bands <- rbind(eigen_bands, cbind(c(get_posterior_eigen2(mcmc_results, evals, zi, alpha)$mean), rep(1:3, each = length(t_paf)), rep(Group, length(t_paf) * 3), rep(Month / 12, length(t_paf) * 3), rep(t_paf, 3)))

Group <- 1
Month <- 132
Month_t <- (Month - mean(demDat$Age)) / sd(demDat$Age)
zi <- c(1, Group, Month_t, Group * Month_t)
eigen_bands <- rbind(eigen_bands, cbind(c(get_posterior_eigen2(mcmc_results, evals, zi, alpha)$mean), rep(1:3, each = length(t_paf)), rep(Group, length(t_paf) * 3), rep(Month / 12, length(t_paf) * 3), rep(t_paf, 3)))

colnames(eigen_bands) <- c("Power", "Eigenfunction", "Group", "Age (years)", "Frequency")
eigen_bands <- as_tibble(eigen_bands)
eigen_t <- NULL
eigen_last <- numeric(length(t_paf))
for(eigenfn in c(1,2,3)){
  for(grp in c(0, 1)){
    eigen1 <- eigen_bands %>%
      filter(Group == grp, Eigenfunction == eigenfn, `Age (years)` == 3) %>%
      select(Power)
    eigen2 <- eigen_bands %>%
      filter(Group == grp, Eigenfunction == eigenfn, `Age (years)` == 7) %>%
      select(Power)
    eigen3 <- eigen_bands %>%
      filter(Group == grp, Eigenfunction == eigenfn, `Age (years)` == 11) %>%
      select(Power)
    if(sum((eigen1 + eigen_last)^2) < sum((eigen1 - eigen_last)^2)){
      eigen1 <- -eigen1
    }
    
    if(sum((eigen2 + eigen1)^2) < sum((eigen2 - eigen1)^2)){
      eigen2 <- -eigen2
    }
    if(sum((eigen3 + eigen1)^2) < sum((eigen3 - eigen1)^2)){
      eigen3 <- -eigen3
    }
    eigen_last <- eigen1
    eigen_t <- rbind(eigen_t, eigen1, eigen2, eigen3)
  }
}

colnames(eigen_t) <- "eigen_t"
eigen_tibble <- eigen_bands %>%
  arrange(Eigenfunction, Group, `Age (years)`) %>%
  bind_cols(eigen_t)
eigen_tibble$`Age (years)` <- as.character(eigen_tibble$`Age (years)`)
eigen_tibble$Eigenfunction <- as.character(eigen_tibble$Eigenfunction)
eigen.labs <- c("Eigenfunction 1", "Eigenfunction 2", "Eigenfunction 3")
names(eigen.labs) <- c("1", "2", "3")

group.labs <- c("TD", "ASD")
names(group.labs) <- c("0", "1")
eigen_tibble %>%
  ggplot(aes(x = Frequency)) +
  geom_line(aes(y = eigen_t, group = `Age (years)`, color = `Age (years)`)) +
  labs(y = "Power", color = "Age (years)") +
  facet_grid(Group ~ Eigenfunction, labeller = labeller(Group = group.labs, Eigenfunction = eigen.labs)) +
  theme_bw()


# Load data and convert dose to a factor variable
data("ToothGrowth")
ToothGrowth$dose <- as.factor(ToothGrowth$dose)

# Box plot, facet accordding to the variable dose and supp
p <- ggplot(ToothGrowth, aes(x = dose, y = len)) + 
  geom_boxplot(aes(fill = supp), position = position_dodge(0.9)) +
  scale_fill_viridis_d() 
p + facet_grid(dose ~ supp)
p + facet_grid(dose ~ supp, labeller = label_both)

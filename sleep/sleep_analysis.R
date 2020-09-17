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
sleep_tabulated <- read.csv(paste0("/Users/johnshamshoian/Documents/R_projects/",
            "bfcr/sleep/tabulated_data/",
            "shhs1-dataset-0.15.0.csv"), stringsAsFactors = FALSE)

num_epochs <- 120
id_range <- 200001:200600

sleep_tabulated_filtered <- sleep_tabulated %>%
  filter(EEG1qual == 4) %>%
  select(nsrrid, age_s1, bmi_s1) %>%
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
epoch_df <- ceiling(10 / 100 * num_epochs)
age_grid <- sleep_data_filtered %>%
  group_by(nsrrid) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(age_s1)
age_df <- ceiling(10 / 100 * (max(age_grid) - min(age_grid)))

epoch_basis <- ps(epoch_grid, df = epoch_df, intercept = TRUE)
epoch_penalty <- attr(epoch_basis, "S")

# Age adjusted
if (TRUE) {
  age_basis <- ps(age_grid$age_s1, df = age_df, intercept = FALSE)
  age_penalty <- attr(age_basis, "S")
  if (!(attr(age_basis, "intercept"))) {
    model_penalties <- tensor.prod.penalties(list(age_penalty, epoch_penalty))
    
    mean_indices <- c(1, 1, 2)
    var_indices <- c(1, 1, 2)
    var_penalty <- list(model_penalties[[1]], model_penalties[[2]], epoch_penalty)
    # var_penalty <- list(epoch_penalty)
    # var_penalty <- list(epoch_penalty, diag(180), diag(180))
    mean_penalty <- list(model_penalties[[1]], model_penalties[[2]], epoch_penalty)
    # mean_penalty <- list(epoch_penalty, diag(180), diag(180))
    design_mean <- cbind(age_basis, 1)
    design_var <- cbind(age_basis, 1)
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

# Not age adjusted
if (FALSE) {
  design_mean <- cbind(rep(1, num_subjects))
  design_var <- cbind(rep(1, num_subjects))
  mean_indices <- c(1)
  var_indices <- c(1)
  mean_penalty <- list(epoch_penalty)
  var_penalty <- list(epoch_penalty)
}
response <- t(matrix(sleep_data_filtered$psd,
                   nrow = num_epochs,
                   ncol = num_subjects))
k <- 12
iter <- 10000
burnin <- 2500
thin <- 1
mcmc_results <- run_mcmc(response, design_mean,
                  design_var, epoch_basis,
                  epoch_grid,
                  mean_penalty, var_penalty,
                  mean_indices, var_indices,
                  k, iter, burnin, thin = thin,
                  var = "unequal")
# calculate_waic(mcmc_results)
saveRDS(mcmc_results, file = paste0("/Users/johnshamshoian/Documents/R_projects/",
                                 "bfcr/sleep/mcmc_output/mcmc_results",
                                 k,
                                 ".rds"))
#
subject_bands <- get_posterior_subject_bands(mcmc_results)
mean_bands <- get_posterior_means(mcmc_results, mcmc_results$data$design_mean[3,])
evals <- k
eigen_bands <- get_posterior_eigen(mcmc_results, evals, mcmc_results$data$design_var[3,])
subj <- 31:39
subject_bands %>%
  filter(id %in% subj) %>%
  ggplot() +
  geom_point(aes(x = time, y = response), alpha = .5) +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper), alpha = 0.5) +
  facet_wrap(. ~ id) +
  theme_bw()

number.labs <- paste0("Eigenfunction ", 1:evals, ": ", 100 * round(eigen_bands$prop_var_explained[2,],2), "%",
                      " (", 100 * round(eigen_bands$prop_var_explained[1,], 2), "% - ",
                      100 * round(eigen_bands$prop_var_explained[3,], 2), "%)")
names(number.labs) <- c("1":evals)
eigen_bands$eigenfunctions %>%
  ggplot(aes(x=time)) +
  facet_wrap(. ~ number, labeller = labeller(number = number.labs), scales = "free") +
  geom_line(aes(y=mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .3) +
  labs(x = "Epoch", y = "Value") +
  theme_bw()

plot_ly() %>%
  add_surface(z =~ eigen_bands$surface)

L <- list()
L$y <- lapply(1:num_subjects, function(i) response[i,])
L$t <- lapply(1:num_subjects, function(i) 1:num_epochs)
res <- FPCA(L$y, L$t, list(dataType = "Dense", methodMuCovEst = "smooth"))

plot_ly() %>%
  add_surface(z =~ res$fittedCov)
mean_bands %>%
  ggplot(mapping = aes(time)) +
  geom_line(aes(y=mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .3) +
  labs(x = "Epoch", y = "Relative delta power spectral density") +
  theme_bw()
plot(res$mu, type = "l")
plot(res$phi[,1], type = "l")
lines(eigen_bands$eigenfunctions$mean[1:120], col = "blue")
plot(res$phi[,2], type = "l")
lines(eigen_bands$eigenfunctions$mean[121:240], col="blue")
cumsum(eigen(res$fittedCov)$values)[1:20] * 100 / sum(eigen(res$fittedCov)$values)
cumsum(eigen_bands$prop_var_explained[2,]) * 100

(eigen(res$fittedCov)$values[1:20] * 100 / sum(eigen(res$fittedCov)$values))[1:k]
eigen_bands$prop_var_explained * 100

beta_list <- matrix(0, nrow = 7500, ncol = 120)
counter <- 1
for (i in 2501:10000) {
  print(i)
  beta_list[counter,] <- epoch_basis %*% mcmc_results$samples$beta[,,i] %*% design_mean[2,]
  counter <- counter + 1
}
matplot(t(beta_list), type = "l")

eigvec_list <- matrix(0, nrow = 7500, ncol = 120)
mycov <- matrix(0, nrow = 120, ncol = 120)
counter <- 1
for (i in 3000:3000) {
  print(i)
  mycov[] <- 0
  for (this_k in 1:k) {
    mycov <- mycov + tcrossprod(epoch_basis %*% mcmc_results$samples$lambda[[i]][,,this_k] %*% design_var[5,] )
  }
  eigvec_list[counter,] <- eigen(mycov)$vectors[,2]
  counter <- counter + 1
}
eigvec_list_copy <- eigvec_list
matplot(t(eigvec_list_copy), type = "l")
initial <- eigvec_list[1,]
for (i in 2:7500) {
  print(i)
#  if (sum((eigvec_list[i,] + initial)^2) < sum((eigvec_list[i,] - initial)^2)) {
  if (eigvec_list[i, 112] > 0) {
    eigvec_list[i,] <- eigvec_list[i,] * (-1)

  }

  initial <- (i - 1) / i * initial + 1 / i * eigvec_list[i,]
}
matplot(t(eigvec_list), type = "l")
plot(colMeans(eigvec_list), type = "l")
lines(res$phi[,2], col = "blue")
lines(-eigen_bands$eigenfunctions$mean[121:240], col = "red")
# #
# L <- list()
# L$y <- lapply(1:num_subjects, function(i) response[i,])
# L$t <- lapply(1:num_subjects, function(i) 1:num_epochs)
# res <- FPCA(L$y, L$t, list(dataType = "Dense", methodMuCovEst = "smooth"))
# CreateCovPlot(res)
# res <- covPACE(t(response), time = 1:120)
# my_bs <- splines::bs(age_grid$age_s1, df = 8)
# dim(predict(my_bs, newx = seq(from = 40, to = 80, by = 1)))
# predict(age_basis, newx = seq(from = 40, to = 80, by = 1))
# ?splines::predict.bs
#
# mysterybasis <- smoothCon(s(Age, bs = "ps", k = 6), data = data.frame(Age))
# spec_age <- s(Age, bs = "ps", k = 6)
# Predict.matrix(mysterybasis, new_age)
# bsbasis <- splines::bs(Age, df = 8)
# psbasis <- dlnm::ps(Age, df = 8)
# age_basis2 <- smooth.construct(s(Age, k = 6, bs = "ps", m = 1), data=data.frame(Age), knots=NULL)
# age_seq <- seq(from = 40, to = 80, by = 1)
# new_age <- data.frame(Age = age_seq)
# new_age_spline <- Predict.matrix(age_basis2, new_age)
# dim(new_age_spline)
# mean_age <- matrix(0, nrow = length(epoch_grid), ncol = dim(new_age_spline)[1])
# for (i in 1:dim(new_age_spline)[1]) {
#   mean_bands <- get_posterior_means(mcmc_results, new_age_spline[i,])
#   mean_age[,i] <- mean_bands$mean
# }
#
# plot_ly() %>%
#   add_surface(z =~ mean_age)
# magn <- matrix(0, nrow = 3,ncol = dim(new_age_spline)[1])
# pve <- magn
# mean_eigen <- matrix(0, nrow = length(epoch_grid), ncol = dim(new_age_spline)[1])
# for (i in 1:dim(new_age_spline)[1]) {
#   eigen_bands <- get_posterior_eigen(mcmc_results, 1, new_age_spline[i,])
#   if (i > 2) {
#     if (sum((mean_eigen[, i-1] +eigen_bands$eigenfunctions$mean)^2) <
#         sum((mean_eigen[, i-1] - eigen_bands$eigenfunctions$mean)^2)) {
#       eigen_bands$eigenfunctions$mean <- -eigen_bands$eigenfunctions$mean
#     }
#   }
#   magn[,i] <- eigen_bands$magnitude
#   pve[,i] <- eigen_bands$prop_var_explained
#   mean_eigen[,i] <- eigen_bands$eigenfunctions$mean
# }
#
#
# plot_ly() %>%
#   add_surface(z =~ mean_eigen)
#
#
# my_ps <- function (x, df = 10, knots = NULL, degree = 3, intercept = FALSE,
#           fx = FALSE, S = NULL, diff = 2)
# {
#   nx <- names(x)
#   x <- as.vector(x)
#   range <- range(x, na.rm = TRUE)
#   nax <- is.na(x)
#   if (nas <- any(nax))
#     x <- x[!nax]
#   if ((degree <- as.integer(degree)) < 1)
#     stop("'degree' must be integer >= 1")
#   if (is.null(knots) || length(knots) == 2L) {
#     nik <- df - degree + 2 - intercept
#     if (nik <= 1)
#       stop("basis dimension too small for b-spline degree")
#     xl <- (if (length(knots) == 2L)
#       min(knots)
#       else min(x)) - diff(range) * 0.001
#     xu <- (if (length(knots) == 2L)
#       max(knots)
#       else max(x)) + diff(range) * 0.001
#     dx <- (xu - xl)/(nik - 1)
#     cat("xl is", xl, " xr is ", xu)
#     knots <- seq(xl - dx * degree, xu + dx * degree, length = nik +
#                    2 * degree)
#     print(knots)
#   }
#   else {
#     df <- length(knots) - degree - 2 + intercept
#     if (df - degree <= 1)
#       stop("basis dimension too small for b-spline degree")
#   }
#   if (any(x < knots[degree + 1] | knots[length(knots) - degree] <
#           x))
#     warning("all obs expected within inner df-degree+int knots")
#   basis <- splineDesign(knots, x, degree + 1, x * 0, TRUE)
#   if (!intercept)
#     basis <- basis[, -1L, drop = FALSE]
#   if (nas) {
#     nmat <- matrix(NA, length(nax), ncol(basis))
#     nmat[!nax, ] <- basis
#     basis <- nmat
#   }
#   if (diff < 1L)
#     stop("'diff' must be an integer >=1")
#   if (fx) {
#     S <- NULL
#   }
#   else if (is.null(S)) {
#     S <- crossprod(diff(diag(ncol(basis) + !intercept), diff = diff))
#     S <- (S + t(S))/2
#     if (!intercept)
#       S <- S[-1L, -1L, drop = FALSE]
#   }
#   else if (any(dim(S) != ncol(basis)))
#     stop("dimensions of 'S' not compatible")
#   dimnames(basis) <- list(nx, seq(ncol(basis)))
#   attributes(basis) <- c(attributes(basis), list(df = df, knots = knots,
#                                                  degree = degree, intercept = intercept, fx = fx, S = S,
#                                                  diff = diff))
#   class(basis) <- c("ps", "matrix")
#   return(basis)
# }
# ps1 <- my_ps(epoch_grid, df = 24, intercept = FALSE)
# ps
# psbasis$m <- c(2,2)
# predict_pspline(psbasis, new_age)
#
# bspline <- function(x, xl, xr, ndx, bdeg) {
#   dx<-(xr-xl)/ndx
#   knots<-seq(xl-bdeg*dx,xr+bdeg*dx,by=dx)
#   print(knots)
#   B <- spline.des(knots, x, bdeg + 1, 0 * x)$design
#   B
# }
# ps
# epoch_basis_bs <- bspline(epoch_grid, .881, 120.119, 21, 3)
# epoch_basis_bs1 <- bs(epoch_grid, knots = attr(epoch_basis, "knots"), intercept = FALSE)
# ps2 <- spline.des(attr(epoch_basis, "knots"), epoch_grid, 4, 0 * epoch_grid)
# ps3 <- bs(epoch_grid, knots = attr(epoch_basis, "knots"), Boundary.knots = NULL, intercept = TRUE)
# ps4

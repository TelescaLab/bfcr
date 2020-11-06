

metrics_agg <- array(0, c(9, 4, 300, 2))
tags <- c("ca", "base", "nocov", "nocovbase")
sample_size <- c(100, 400)
for (n in 1:2) {
  for (i in 1:300) {
    for (j in 1:4) {
      load(paste0("/Users/johnshamshoian/Documents/R_projects/bfcr/Simulations/Metrics/n", sample_size[n], "_", tags[j],"_seed", i, ".RData"))
      # metrics_agg[, j, i, n] <- unlist(metrics)
      metrics_agg[1, j, i, n] <- metrics$subj_coverage
      metrics_agg[2, j, i, n] <- metrics$subj_error
      metrics_agg[3, j, i, n] <- metrics$subj_width
      metrics_agg[4, j, i, n] <- metrics$mean_coverage
      metrics_agg[5, j, i, n] <- metrics$mean_error
      metrics_agg[6, j, i, n] <- metrics$mean_width
      metrics_agg[7, j, i, n] <- metrics$covariance_coverage
      metrics_agg[8, j, i, n] <- metrics$covariance_error
      metrics_agg[9, j, i, n] <- metrics$covariance_width
    }
  }
}

round(apply(metrics_agg[,,,2], c(1,2), median), 2)

ts_data <- read.csv("/Users/johnshamshoian/Documents/Lucid-Circuit/October/ts_all.csv", header = FALSE)
ts_data <- as.matrix(ts_data)
plot(seq(from = 1, to = 60, length.out = 18000), ts_data[,9], type = "l")

library(tidyverse)
ts_data_tibble <- tibble(node = rep(1:4, each = 18000 * 3),
                  channel = rep(rep(c("Power", "Compute", "Temperature"), each = 18000), 4),
                  Value = c(ts_data),
                  time = rep(seq(from = 0, to = 300, length.out = 18000), 12))
ts_data_tibble %>%
  ggplot() +
  geom_line(aes(time, Value)) +
  facet_wrap(~channel + node, scales = "free") +
  labs(x = "Time (seconds)", title = "Telemetry view from a single node") +
  theme_bw()


relative_psd2 <- tibble(nsrrid = rep(NA, 4180205),
                        epoch = rep(NA, 4180205),
                        psd = rep(NA, 4180205))
start <- 1
for (i in 5792:5804) {
  print(i)
  id <- 200000 + i
  file_name <- paste0("spectrum/shhs1-", id, ".RData")
  tryCatch(load(file_name), error = function(e) NULL)
  num_epochs <- dim(psd_tibble)[1]
  end <- start + num_epochs - 1
  relative_psd2$nsrrid[start:end] <- rep(id, num_epochs)
  relative_psd2$epoch[start:end] <- 1:num_epochs
  relative_psd2$psd[start:end] <- psd_tibble$relative_power
  start <- end + 1
}
file_name <- paste0("/Users/johnshamshoian/Documents/R_projects/",
                    "bfcr/sleep/data/relative_psd2.RData")
save(relative_psd2, file = file_name)

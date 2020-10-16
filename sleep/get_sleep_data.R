library(luna)
library(tidyverse)
library(future.apply)
library(nsrr)
library(filesstrings)

setwd("/Users/johnshamshoian/Documents/R_projects/bfcr/sleep")
rdata <- dir("spectrum")
done <- as.numeric(substr(rdata, 9, 12))
to_do <- setdiff(1:5900, done)

download_shhs_file <- function(x) {
  id <- 200000 + x
  new_dir <- paste0("shhs1-", x)
  dir.create(new_dir)
  edf_path <- paste0("polysomnography/edfs/shhs1/shhs1-", id, ".edf")
  edf_file <- nsrr_download_file("shhs", path = edf_path)
  file.move(edf_file$outfile, new_dir)
  old_name <- str_split(edf_file$outfile, pattern = "/")[[1]][9]
  new_name <- paste0("shhs1-", id, ".edf")
  old_file_name <- paste0("/Users/johnshamshoian/Documents/R_projects/bfcr/sleep/",
                          new_dir, "/", old_name)
  new_file_name <- paste0("/Users/johnshamshoian/Documents/R_projects/bfcr/sleep/",
                          new_dir, "/", new_name)
  file.rename(old_file_name, new_file_name)
  
  xml_path <- paste0("polysomnography/annotations-events-profusion/shhs1/shhs1-",
                     id, "-profusion.xml")
  xml_file <- nsrr_download_file("shhs", path = xml_path)
  file.move(xml_file$outfile, new_dir)
  old_name <- str_split(xml_file$outfile, pattern = "/")[[1]][9]
  new_name <- paste0("shhs1-", id, "-profusion.xml")
  old_file_name <- paste0("/Users/johnshamshoian/Documents/R_projects/bfcr/sleep/",
                          new_dir, "/", old_name)
  new_file_name <- paste0("/Users/johnshamshoian/Documents/R_projects/bfcr/sleep/",
                          new_dir, "/", new_name)
  file.rename(old_file_name, new_file_name)
  
  system(paste0("luna --build shhs1-", x, " -ext=-profusion.xml > lst_files/", x, ".lst"))
  sl <- lsl(paste0("lst_files/", x, ".lst"))
  suppressWarnings(lattach(sl, paste0("shhs1-", id)))
  leval("SIGNALS keep=EEG")
  leval("MASK all")
  leval("MASK unmask-if=${sleep}")
  leval("FILTER bandpass=0.3,35 ripple=0.02 tw=0.5")
  leval("ARTIFACTS mask & SIGSTATS epoch mask threshold=2,2")
  leval("RE")
  k <- leval("PSD epoch-spectrum") 
  k_tibble <- as_tibble(k$PSD$CH_E_F)
  k_tibble <- k_tibble %>%
    mutate(delta_ind = if_else(F >= .5 & F <= 4, 1, 0))
  psd_tibble <- k_tibble %>%
    group_by(E) %>%
    summarise(delta_power = sum(PSD[delta_ind == 1]), total = sum(PSD),
              relative_power = delta_power / total)
  save(psd_tibble, file = 
         paste0("spectrum/shhs1-", id, ".RData"))
  system(paste0("rm -rf ", new_dir))

}

error_checking_download_shhs <- function(i) {
  tryCatch(download_shhs_file(i), error = function(e) NULL)
}
plan(multiprocess, workers = 8)
future_lapply(to_do, error_checking_download_shhs)
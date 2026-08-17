### Setup script for running the code in this repository.
###
### Usage:
### 1. Clone or download the repository.
### 2. Set the working directory to the repository root.
### 3. Run:
###       source("setup.R")
###
### This script copies the synthetic pregnancy dataset to the folders
### expected by the Monolix projects.

if (!exists("REPO")) {
  REPO <- getwd()
}

syn <- file.path(REPO, "data", "pregnancy_synthetic.csv")

if (!file.exists(syn)) {
  stop(
    "Synthetic dataset not found at: ",
    syn,
    "\nCheck that REPO points to the repository root."
  )
}

for (m in c("Model_1", "Model_2", "Model_3")) {
  
  dest_dir <- file.path(REPO, "monolix", "data1", m)
  
  if (!dir.exists(dest_dir)) {
    stop("Monolix project directory not found: ", dest_dir)
  }
  
  dest <- file.path(dest_dir, "jointmodel_data_new.csv")
  
  file.copy(
    from = syn,
    to = dest,
    overwrite = TRUE
  )
  
  cat("Copied to:", dest, "\n")
}

cat("\nSetup completed successfully.\n")
cat("The synthetic pregnancy data are ready for the Monolix data1 projects.\n")
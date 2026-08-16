### Section 3.3 -- Sensitivity of SAEM to the initial values (multistart study).
###
### Replaces the former sensibilidad_valores_iniciales_M{1,2,3}.R. The fitting loop was
### identical across the three; what differed is the model configuration and the scheme used
### to build the starting values, both collected in MODELS below.
###
### Referee point 1 of SMM-24-0694: "How does the selection of initial values impact this
### complex model? ... simulation studies to investigate the impact of initial values on
### convergence?" The project is refitted from n_starts sets of initial values (near, far --
### including sign flips of the classification coefficients -- and the naive values that the
### .mlxtran shipped with), and the final estimates and log-likelihood are compared across runs.
### If all runs converge to the same point, SAEM is robust to the starting point; if they do
### not, there is genuine multimodality / weak identifiability to report.
###
### Usage (REPO = path to your clone, see README):
###   REPO  <- "/path/to/joint-model-nlme-saem"
###   MODEL <- "Model_1"                                  # or "Model_2" / "Model_3"
###   source(file.path(REPO, "R/simulation/02_initial_values.R"))

library(lixoftConnectors)
initializeLixoftConnectors(software = "monolix", force = TRUE)

## --- starting-value schemes (kept verbatim from the original scripts) ---
## Scheme A (Model 1): multiply the fitted values, flipping the sign of alpha_pop/beta_pop
## in the "far" starts. Scheme B (Models 2 and 3): treat the longitudinal and classification
## blocks separately, drawing the classification coefficients freely in the far starts.
starts_scheme_A <- function(cfg, n_starts) {
  fit <- cfg$fit
  starts <- vector("list", n_starts)
  for (k in 1:n_starts) {
    if (k <= 8) {
      starts[[k]] <- fit * runif(length(fit), 0.5, 1.5)          # near: +/-50%
    } else if (k <= 16) {
      f <- runif(length(fit), 0.1, 4)                            # far: wide factor + sign flip
      starts[[k]] <- fit * f
      starts[[k]]["alpha_pop"] <- -starts[[k]]["alpha_pop"]
      starts[[k]]["beta_pop"]  <- -starts[[k]]["beta_pop"]
    } else {
      starts[[k]] <- cfg$naive * runif(length(fit), 0.3, 3)      # naive .mlxtran values
    }
    for (nm in c("omega_a", "omega_b", "omega_c")) starts[[k]][nm] <- abs(starts[[k]][nm])
  }
  starts
}
starts_scheme_B <- function(cfg, n_starts) {
  fit_long <- cfg$fit_long; fit_clas <- cfg$fit_clas
  naive_long <- cfg$naive_long; naive_clas <- cfg$naive_clas
  starts <- vector("list", n_starts)
  for (k in 1:n_starts) {
    if (k <= 8) {
      long_k <- fit_long * runif(length(fit_long), 0.7, 1.3)
      clas_k <- fit_clas + rnorm(length(fit_clas), sd = 0.3 * abs(fit_clas) + 1)
    } else if (k <= 16) {
      long_k <- fit_long * runif(length(fit_long), 0.1, 4)
      clas_k <- runif(length(fit_clas), -150, 150)
      names(clas_k) <- names(fit_clas)
    } else {
      long_k <- naive_long * runif(length(fit_long), 0.3, 3)
      clas_k <- naive_clas + rnorm(length(fit_clas), sd = 10)
    }
    for (nm in c("omega_a", "omega_b")) long_k[nm] <- abs(long_k[nm])
    starts[[k]] <- c(long_k, clas_k)
  }
  starts
}

## --- model configurations ---
MODELS <- list(
  Model_1 = list(
    project = "Model.1.mlxtran", scheme = starts_scheme_A,
    fit = c(a_pop = 4.54, b_pop = 15.6, c_pop = 7.22, alpha_pop = 34.3, beta_pop = -7.94,
            omega_a = 0.502, omega_b = 3.31, omega_c = 1.95),
    naive = c(a_pop = 4, b_pop = 15, c_pop = 7, alpha_pop = -15, beta_pop = 3,
              omega_a = 1, omega_b = 1, omega_c = 1)),
  Model_2 = list(
    project = "Model.2n.mlxtran", scheme = starts_scheme_B,
    fit_long = c(a_pop = 4.54, b_pop = 15.6, c_pop = 7, omega_a = 0.468, omega_b = 4.35),
    fit_clas = c(alpha_pop = 28.3, beta1_pop = -6.57, beta2_pop = 2.63e-07),
    naive_long = c(a_pop = 4, b_pop = 15, c_pop = 7, omega_a = 1, omega_b = 1),
    naive_clas = c(alpha_pop = 5, beta1_pop = -8, beta2_pop = -3)),
  Model_3 = list(
    project = "Model.3ln.mlxtran", scheme = starts_scheme_B,
    fit_long = c(a_pop = 4.55, b_pop = 15, c_pop = 7.18, omega_a = 0.0794, omega_b = 0.283),
    fit_clas = c(alpha_pop = 46.9, beta1_pop = -11.1, beta2_pop = 0.0892),
    naive_long = c(a_pop = 4, b_pop = 15, c_pop = 7, omega_a = 1, omega_b = 1),
    naive_clas = c(alpha_pop = 5, beta1_pop = -8, beta2_pop = -3))
)
if (!exists("MODEL")) MODEL <- "Model_1"
stopifnot(MODEL %in% names(MODELS))
cfg <- MODELS[[MODEL]]
tag <- sub("Model_", "M", MODEL)

setwd(file.path(REPO, "monolix", "data1", MODEL))
project <- cfg$project

set.seed(2026)
n_starts <- 20
starts <- cfg$scheme(cfg, n_starts)

## --- refit from each starting point ---
resultados <- data.frame()
for (k in 1:n_starts) {
  cat("=== ", MODEL, ": start ", k, " of ", n_starts, " ===\n", sep = "")
  res_k <- tryCatch({
    loadProject(project)

    popparams <- getPopulationParameterInformation()
    ini <- starts[[k]]
    for (nm in names(ini)) popparams$initialValue[popparams$name == nm] <- as.numeric(ini[nm])
    setPopulationParameterInformation(popparams)

    scenario <- getScenario()
    scenario$tasks <- c(populationParameterEstimation = TRUE,
                        conditionalModeEstimation = FALSE,
                        conditionalDistributionSampling = FALSE,
                        standardErrorEstimation = TRUE,
                        logLikelihoodEstimation = TRUE)
    setScenario(scenario)

    t0 <- Sys.time()
    runScenario()
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    est <- getEstimatedPopulationParameters()
    ll  <- getEstimatedLogLikelihood()   # NOTE: nested list, flatten with unlist()

    ini_df <- as.data.frame(as.list(setNames(as.numeric(ini), paste0("ini_", names(ini)))))
    est_df <- as.data.frame(as.list(setNames(as.numeric(est[names(ini)]), paste0("est_", names(ini)))))
    ll_df  <- as.data.frame(as.list(unlist(ll)))

    cbind(data.frame(start = k, status = "OK"), ini_df, est_df, ll_df,
          data.frame(elapsed_sec = elapsed))
  }, error = function(e) data.frame(start = k, status = paste("ERROR:", conditionMessage(e))))
  resultados <- dplyr::bind_rows(resultados, res_k)
}

out <- paste0("sensibilidad_valores_iniciales_", tag, "_resultados.csv")
write.csv(resultados, out, row.names = FALSE)
cat("\nDone. Results saved to ", out, "\n", sep = "")

## --- dispersion of the final estimates across successful runs ---
## Small sd/|mean| for every parameter -> SAEM robust to the starting point.
## Large sd for the classification coefficients -> confirms the weak identifiability.
ok <- resultados[resultados$status == "OK", ]
cat("\nSuccessful runs:", nrow(ok), "of", n_starts, "\n")
if (nrow(ok) > 1) {
  est_cols <- grep("^est_", names(ok), value = TRUE)
  resumen <- sapply(ok[est_cols], function(x) c(media = mean(x), sd = sd(x),
                                                min = min(x), max = max(x)))
  cat("\nFinal estimates across the", nrow(ok), "successful runs:\n")
  print(round(t(resumen), 4))
  cat("\nRange of OFV (-2LL) across runs:", paste(range(ok$OFV), collapse = " to "), "\n")
}

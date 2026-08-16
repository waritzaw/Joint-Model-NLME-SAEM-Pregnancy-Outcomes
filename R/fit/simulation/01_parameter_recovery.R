### Section 3.1 -- Parameter recovery (bias, RMSE, empirical 95% CI coverage). Table 1.
###
### Replaces the former recuperacion_parametros_M{1,2,3}.R, which were identical except for
### the model configuration collected in MODELS below.
###
### R datasets are simulated from the fitted ("true") parameters of each model, using the real
### visit design, and refitted by SAEM from reasonable (+/-30% jittered) starting values.
### Coverage is the proportion of replicates whose 95% Wald interval (estimate +/- 1.96 SE,
### with SE from the Fisher information) contains the true value.
###
### Usage (REPO = path to your clone, see README):
###   REPO  <- "/path/to/joint-model-nlme-saem"
###   MODEL <- "Model_1"                                   # or "Model_2" / "Model_3"
###   source(file.path(REPO, "R/simulation/01_parameter_recovery.R"))
###
### Requires Monolix 2024R1 + lixoftConnectors.

library(lixoftConnectors)
initializeLixoftConnectors(software = "monolix", force = TRUE)

## --- model configurations (the only thing that differed between the three scripts) ---
## re_dist  : distribution of the random effects a_i, b_i
## c_random : TRUE if c has a random effect (omega_c), FALSE if c is fixed across subjects
## covars   : random effects entering the GLM, with their coefficient names in `truth`
MODELS <- list(
  Model_1 = list(
    project = "Model.1.mlxtran", re_dist = "normal", c_random = TRUE,
    truth = c(a_pop = 4.54, b_pop = 15.6, c_pop = 7.22, alpha_pop = 34.3, beta_pop = -7.94,
              omega_a = 0.502, omega_b = 3.31, omega_c = 1.95),
    sigma_resid = 0.261, coefs = c(beta_pop = "a")),
  Model_2 = list(
    project = "Model.2n.mlxtran", re_dist = "normal", c_random = FALSE,
    truth = c(a_pop = 4.54, b_pop = 15.6, c_pop = 7.0, alpha_pop = 28.3,
              beta1_pop = -6.57, beta2_pop = 2.63e-07, omega_a = 0.468, omega_b = 4.35),
    sigma_resid = 0.266, coefs = c(beta1_pop = "a", beta2_pop = "b")),
  Model_3 = list(
    project = "Model.3ln.mlxtran", re_dist = "lognormal", c_random = FALSE,
    truth = c(a_pop = 4.55, b_pop = 15.0, c_pop = 7.18, alpha_pop = 46.9,
              beta1_pop = -11.1, beta2_pop = 0.0892, omega_a = 0.0794, omega_b = 0.283),
    sigma_resid = 0.30, coefs = c(beta1_pop = "a", beta2_pop = "b"))
)
if (!exists("MODEL")) MODEL <- "Model_1"
stopifnot(MODEL %in% names(MODELS))
cfg <- MODELS[[MODEL]]
truth <- cfg$truth; sigma_resid <- cfg$sigma_resid
tag <- sub("Model_", "M", MODEL)

setwd(file.path(REPO, "monolix", "data1", MODEL))
project <- cfg$project

## --- 1. Real design: visit times per id, plus the ytype=2 row (time + group) ---
real <- read.csv("jointmodel_data_new.csv", colClasses = c(group = "character"))
ids <- sort(unique(real$id))
design <- lapply(ids, function(i) {
  di <- real[real$id == i, ]
  list(times_long = di$time[di$ytype == 1],
       time_bin   = di$time[di$ytype == 2],
       group      = di$group[1])
})
names(design) <- ids

## --- 2. Simulate one dataset from the model's own "truth" ---
simular_dataset <- function(truth, sigma_resid, design, seed) {
  set.seed(seed)
  n <- length(design)
  if (cfg$re_dist == "lognormal") {
    ## a_i = median * exp(N(0, log-SD))
    a_i <- as.numeric(truth["a_pop"]) * exp(rnorm(n, 0, truth["omega_a"]))
    b_i <- as.numeric(truth["b_pop"]) * exp(rnorm(n, 0, truth["omega_b"]))
  } else {
    a_i <- rnorm(n, truth["a_pop"], truth["omega_a"])
    b_i <- rnorm(n, truth["b_pop"], truth["omega_b"])
  }
  if (cfg$c_random) {
    c_i <- pmax(rnorm(n, truth["c_pop"], truth["omega_c"]), 0.1)  # avoid c<=0
  } else {
    c_i <- rep(as.numeric(truth["c_pop"]), n)                     # c fixed across subjects
  }

  filas <- vector("list", n)
  for (k in seq_len(n)) {
    d <- design[[k]]
    t_long <- d$times_long
    E <- a_i[k] / (1 + exp(-(t_long - b_i[k]) / c_i[k]))
    y_long <- E + rnorm(length(t_long), 0, sigma_resid)

    ## GLM linear predictor: intercept + one term per random effect used as covariate
    re_k <- c(a = a_i[k], b = b_i[k])
    lgp1 <- truth["alpha_pop"] + sum(mapply(function(cf, re) truth[cf] * re_k[[re]],
                                            names(cfg$coefs), cfg$coefs))
    level_sim <- rbinom(1, 1, 1 / (1 + exp(-lgp1)))

    filas[[k]] <- rbind(
      data.frame(id = names(design)[k], time = t_long, y = y_long,
                 group = d$group, ytype = 1),
      data.frame(id = names(design)[k], time = d$time_bin, y = level_sim,
                 group = d$group, ytype = 2)
    )
  }
  do.call(rbind, filas)
}

## --- 3. Simulate + refit loop ---
set.seed(2026)
R <- 50
resultados <- data.frame()

for (r in 1:R) {
  cat("=== ", MODEL, ": replicate ", r, " of ", R, " ===\n", sep = "")
  res_r <- tryCatch({
    sim_data <- simular_dataset(truth, sigma_resid, design, seed = 10000 + r)
    fname <- paste0("sim_data_", r, ".csv")
    write.csv(sim_data, fname, row.names = FALSE, quote = c(4))  # quote "group" only

    loadProject(project)
    BaseData <- getData()
    setData(fname, headerTypes = BaseData$headerTypes, observationTypes = BaseData$observationTypes)

    ## reasonable starting values: +/-30% jitter around the truth (not adversarial)
    popparams <- getPopulationParameterInformation()
    ini <- truth * runif(length(truth), 0.7, 1.3)
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
    se  <- tryCatch(getEstimatedStandardErrors()$stochasticApproximation, error = function(e) NULL)
    ll  <- getEstimatedLogLikelihood()   # NOTE: nested list, flatten with unlist()

    est_df <- as.data.frame(as.list(setNames(as.numeric(est[names(truth)]), paste0("est_", names(truth)))))
    if (!is.null(se) && all(names(truth) %in% se$parameter)) {
      se_vec <- setNames(se$se[match(names(truth), se$parameter)], names(truth))
    } else {
      se_vec <- setNames(rep(NA_real_, length(truth)), names(truth))
    }
    se_df <- as.data.frame(as.list(setNames(as.numeric(se_vec), paste0("se_", names(truth)))))
    ll_df <- as.data.frame(as.list(unlist(ll)))

    cbind(data.frame(rep = r, status = "OK"), est_df, se_df, ll_df,
          data.frame(elapsed_sec = elapsed))
  }, error = function(e) data.frame(rep = r, status = paste("ERROR:", conditionMessage(e))))
  resultados <- dplyr::bind_rows(resultados, res_r)
  unlink(paste0("sim_data_", r, ".csv"))
}

out_res <- paste0("recuperacion_parametros_", tag, "_resultados.csv")
write.csv(resultados, out_res, row.names = FALSE)
cat("\nDone. Results saved to ", out_res, "\n", sep = "")

## --- 4. Bias, RMSE and empirical 95% CI coverage ---
ok <- resultados[resultados$status == "OK", ]
cat("\nSuccessful replicates:", nrow(ok), "of", R, "\n")

resumen <- data.frame()
for (nm in names(truth)) {
  est_col <- ok[[paste0("est_", nm)]]
  se_col  <- ok[[paste0("se_", nm)]]
  sesgo <- mean(est_col, na.rm = TRUE) - truth[nm]
  rmse  <- sqrt(mean((est_col - truth[nm])^2, na.rm = TRUE))
  cobertura <- if (all(!is.na(se_col))) {
    mean(truth[nm] >= (est_col - 1.96 * se_col) & truth[nm] <= (est_col + 1.96 * se_col), na.rm = TRUE)
  } else NA
  resumen <- rbind(resumen, data.frame(parametro = nm, verdad = truth[nm],
                                       media_est = mean(est_col, na.rm = TRUE),
                                       sesgo = sesgo,
                                       sesgo_relativo_pct = 100 * sesgo / truth[nm],
                                       rmse = rmse, cobertura_95 = cobertura))
}
cat("\n--- Parameter recovery (", nrow(ok), " replicates) ---\n", sep = "")
print(resumen, row.names = FALSE)
write.csv(resumen, paste0("recuperacion_parametros_", tag, "_resumen.csv"), row.names = FALSE)

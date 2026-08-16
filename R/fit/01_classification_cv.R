### Leave-one-out cross-validated classification for the pregnancy application (data1).
###
### Replaces the former ClassCV_M1.R / ClassCV_M2.R / ClassCV_M3.R, which were identical
### except for the model configuration collected in MODELS below.
###
### For each woman i: refit the joint model on the other 172 subjects, fix the population
### parameters, sample the conditional distribution of her random effects, and predict her
### probability of an abnormal pregnancy from the GLM component. A confusion matrix and the
### usual accuracy metrics are reported at the end.
###
### Usage (REPO = path to your clone, see README):
###   REPO <- "/path/to/joint-model-nlme-saem"
###   MODEL <- "Model_1"                                  # or "Model_2" / "Model_3"
###   source(file.path(REPO, "R/fit/01_classification_cv.R"))
###
### Requires Monolix 2024R1 + lixoftConnectors (install it once, from your MonolixSuite
### installation; this script does not install anything).

library(dplyr)
library(lixoftConnectors)
initializeLixoftConnectors(software = "monolix", force = TRUE)

## --- model configurations (the only thing that differed between the three scripts) ---
MODELS <- list(
  Model_1 = list(project = "Model.1.mlxtran",   n_slopes = 1),  # logit(P) = eta0 + eta1*a
  Model_2 = list(project = "Model.2n.mlxtran",  n_slopes = 2),  # logit(P) = eta0 + eta1*a + eta2*b, normal REs
  Model_3 = list(project = "Model.3ln.mlxtran", n_slopes = 2)   # idem, lognormal REs
)
if (!exists("MODEL")) MODEL <- "Model_1"
stopifnot(MODEL %in% names(MODELS))
cfg <- MODELS[[MODEL]]
tag <- sub("Model_", "M", MODEL)

## Working directory: by default the model folder inside the clone. Set WORKDIR beforehand to
## run against a different copy of the project (e.g. a private folder holding the real data,
## so that confidential data never has to be placed inside the repository).
if (!exists("WORKDIR")) WORKDIR <- file.path(REPO, "monolix", "data1", MODEL)
setwd(WORKDIR)
project <- cfg$project
data <- read.csv("jointmodel_data_new.csv")
cat(sprintf("Working directory : %s\nData              : %d rows, %d subjects\n",
            WORKDIR, nrow(data), length(unique(data$id))))

loadProject(project)
pred <- NULL
aux.id <- unique(data$id)

for (i in seq_along(aux.id)) {
  cat(sprintf("=== %s: leaving out subject %d/%d ===\n", MODEL, i, length(aux.id)))
  loadProject(project)

  ## --- fit on all subjects except i ---
  data.aux <- data[-which(data$id == i), ]
  write.csv(data.aux, paste0("dataAux", tag, ".csv"), row.names = FALSE)
  BaseData <- getData()
  setData(paste0("dataAux", tag, ".csv"),
          headerTypes = BaseData$headerTypes, observationTypes = BaseData$observationTypes)

  saveProject(paste0("Modelo-lineal-", i, ".mlxtran"))
  scenario <- getScenario()
  scenario$tasks <- c(populationParameterEstimation = TRUE,
                      conditionalModeEstimation = TRUE,
                      conditionalDistributionSampling = TRUE,
                      standardErrorEstimation = FALSE,
                      logLikelihoodEstimation = FALSE)
  setScenario(scenario)
  runScenario()

  ## --- fix the population parameters, then sample subject i's random effects ---
  setInitialEstimatesToLastEstimates()
  popparams <- getPopulationParameterInformation()
  popparams$method <- "FIXED"
  setPopulationParameterInformation(popparams)

  data.new <- data[data$id == i, ]
  write.csv(data.new, paste0("dataNew", tag, ".csv"), row.names = FALSE)
  BaseData <- getData()
  setData(paste0("dataNew", tag, ".csv"),
          headerTypes = BaseData$headerTypes, observationTypes = BaseData$observationTypes)
  setConditionalDistributionSamplingSettings(nbminiterations = 500, ratio = 1,
                                             nbsimulatedparameters = 200)

  runPopulationParameterEstimation()   # mandatory before the other tasks (all parameters fixed)
  runConditionalDistributionSampling() # sample from the posterior of the random effects

  simparams <- getSimulatedIndividualParameters()
  simparams$id  <- row.names(simparams)
  simparams$rep <- NULL
  write.csv(simparams, file = paste0("table_simulated_parameters_", tag, ".csv"), row.names = FALSE)

  ## E[phi_i | Y_ij, theta_hat] with phi_i = (a_i, b_i, c_i)
  re <- colMeans(simparams[, 2:4])
  ## linear predictor: intercept popparams[4,2] + one term per slope
  lp <- popparams[4, 2] + sum(sapply(seq_len(cfg$n_slopes),
                                     function(j) popparams[4 + j, 2] * re[j]))
  pred <- rbind(pred, 1 / (1 + exp(-lp)))
}

write.csv(pred, paste0("pred-", tag, "-CV.csv"), row.names = FALSE)

## --- classification at the 0.5 threshold + accuracy metrics ---
pr    <- read.csv(paste0("pred-", tag, "-CV.csv"))
res.c <- ifelse(pr[[1]] > 0.5, 1, 0)
obs   <- data$group[cumsum(table(data[data$ytype == 1, 1]))]

pred_f <- factor(res.c, levels = c(0, 1), labels = c("Normal", "Abnormal"))
obs_f  <- factor(obs,   levels = c(0, 1), labels = c("Normal", "Abnormal"))

## Raw confusion counts (rows = observed, columns = predicted). Everything below is a
## deterministic function of this table, so it is saved for the record.
tab <- table(Observed = obs_f, Predicted = pred_f)
cat("\n=== Confusion matrix (rows = observed, columns = predicted) ===\n")
print(tab)
write.csv(as.data.frame.matrix(tab), paste0("confusion_", tag, "_CV.csv"))

library(caret)
## caret's convention is confusionMatrix(data = PREDICTED, reference = OBSERVED).
## Passing them the other way round transposes the table, which exchanges sensitivity with
## precision and specificity with the negative predictive value (accuracy, F1 and kappa are
## unaffected). We report both positive-class conventions so the labels are unambiguous.
for (pos in c("Abnormal", "Normal")) {
  CM <- confusionMatrix(data = pred_f, reference = obs_f, positive = pos)
  cat("\n=== Accuracy metrics -- positive class = ", pos, " ===\n", sep = "")
  print(round(c(CM$overall[c("Accuracy", "Kappa")],
                CM$byClass[c("Sensitivity", "Specificity", "Precision", "Recall", "F1")],
                "Error rate" = unname(1 - CM$overall["Accuracy"]),
                "McNemar p" = unname(CM$overall["McnemarPValue"])), 4))
  if (pos == "Abnormal") {
    met <- data.frame(metric = c("Accuracy", "Kappa", "Sensitivity", "Specificity",
                                 "Precision", "F1", "ErrorRate"),
                      value = round(c(CM$overall["Accuracy"], CM$overall["Kappa"],
                                      CM$byClass[c("Sensitivity", "Specificity", "Precision", "F1")],
                                      1 - CM$overall["Accuracy"]), 4))
    write.csv(met, paste0("metrics_", tag, "_CV.csv"), row.names = FALSE)
  }
}
cat("\nSaved: confusion_", tag, "_CV.csv  and  metrics_", tag, "_CV.csv\n", sep = "")

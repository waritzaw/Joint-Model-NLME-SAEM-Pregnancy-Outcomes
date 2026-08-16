### Comparacion de estructuras de error residual REM I-IV -- Modelo I (embarazo)
### Referato SMM-24-0694, punto 7: justificar con AIC/BIC por que se elige REM I.
###
### Ajusta el MISMO modelo (curva logistica de crecimiento, Modelo I con a,b,c
### aleatorios + clasificacion logit(P)=eta0+eta1*a) cambiando SOLO el errorModel
### de la observacion continua y1, a las 4 estructuras del paper:
###   REM I   -> constant       y = mu + a*eps
###   REM II  -> proportional    y = mu + b*mu*eps        (d=1)
###   REM III -> combined1       y = mu + (a + b*mu)*eps
###   REM IV  -> combined2       y = mu + sqrt(a^2+b^2*mu^2)*eps
### y reporta -2LL, AIC y BIC de cada una.
###
### NOTA: no hace saveProject, asi que el proyecto original Model.1.mlxtran NO se
### modifica; solo se sobrescribe la carpeta de resultados 'Model.1' en cada corrida
### (extraemos la verosimilitud inmediatamente despues de cada run).

library(lixoftConnectors)
initializeLixoftConnectors(software = "monolix", force = TRUE)

setwd(file.path(REPO, "monolix", "data1", "Model_1"))   # REPO = path to your clone (define it once; see README)
project <- "Model.1.mlxtran"

## REM del paper -> nombre de errorModel en Monolix
error_models <- c(REM_I = "constant", REM_II = "proportional",
                  REM_III = "combined1", REM_IV = "combined2")

## Valores iniciales buenos para los parametros estructurales + aleatorios
## (del ajuste ya reportado, Model.1/summary.txt) para que las 4 corridas partan
## del mismo lugar razonable y el AIC/BIC sea comparable. Los parametros del error
## se dejan en el valor por defecto que asigna setErrorModel.
ini_struct <- c(a_pop = 4.54, b_pop = 15.6, c_pop = 7.22, alpha_pop = 34.3,
                beta_pop = -7.94, omega_a = 0.502, omega_b = 3.31, omega_c = 1.95)

resultados <- data.frame()

for (k in seq_along(error_models)) {
  nombre <- names(error_models)[k]
  em <- as.character(error_models[k])
  cat("=== ", nombre, " (errorModel = ", em, ") ===\n", sep = "")
  res <- tryCatch({
    loadProject(project)

    ## cambiar SOLO el modelo de error de la observacion continua y1
    setErrorModel(y1 = em)

    if (k == 1) {
      cat("--- diagnostico (solo 1a corrida) ---\n")
      cat("Modelo de observacion tras setErrorModel:\n")
      print(getContinuousObservationModel())
      cat("Parametros poblacionales disponibles:\n")
      print(getPopulationParameterInformation()$name)
    }

    ## fijar valores iniciales de los parametros estructurales/aleatorios presentes
    pp <- getPopulationParameterInformation()
    for (nm in names(ini_struct)) {
      if (nm %in% pp$name) pp$initialValue[pp$name == nm] <- as.numeric(ini_struct[nm])
    }
    setPopulationParameterInformation(pp)

    ## escenario: estimacion poblacional + log-verosimilitud (para AIC/BIC)
    sc <- getScenario()
    sc$tasks <- c(populationParameterEstimation = TRUE,
                  conditionalModeEstimation = FALSE,
                  conditionalDistributionSampling = FALSE,
                  standardErrorEstimation = FALSE,
                  logLikelihoodEstimation = TRUE)
    setScenario(sc)

    t0 <- Sys.time()
    runScenario()
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    ll <- unlist(getEstimatedLogLikelihood())
    getv <- function(suf) {
      v <- ll[grepl(paste0(suf, "$"), names(ll))]
      if (length(v) > 0) as.numeric(v[1]) else NA_real_
    }
    ## nb de parametros estimados (para verificar AIC = -2LL + 2k)
    k_par <- sum(getPopulationParameterInformation()$method == "MLE")

    data.frame(REM = nombre, errorModel = em, status = "OK",
               neg2LL = getv("OFV"), AIC = getv("AIC"), BIC = getv("BIC"),
               n_par = k_par, elapsed_sec = elapsed)
  }, error = function(e) {
    data.frame(REM = nombre, errorModel = em, status = paste("ERROR:", conditionMessage(e)))
  })
  resultados <- dplyr::bind_rows(resultados, res)
}

write.csv(resultados, "REM_comparison_data1_resultados.csv", row.names = FALSE)
cat("\nListo. Guardado en REM_comparison_data1_resultados.csv\n\n")
print(resultados)

## resumen ordenado por AIC (el mejor primero)
ok <- resultados[resultados$status == "OK", ]
if (nrow(ok) > 1) {
  cat("\n--- Ordenado por AIC (menor = mejor) ---\n")
  print(ok[order(ok$AIC), c("REM", "errorModel", "neg2LL", "AIC", "BIC")], row.names = FALSE)
}

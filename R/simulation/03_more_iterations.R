### Prueba: ¿mas iteraciones resuelve la degeneracion, o es un problema de
### identificabilidad? -- Modelo 1
###
### Tomamos 3 semillas ya conocidas del estudio de sensibilidad
### (sensibilidad_valores_iniciales_M1_resultados.csv) y las re-corremos con el
### techo de iteraciones multiplicado por 10x. Si las estimaciones finales no
### cambian (siguen degeneradas/divergentes), confirma que el "autostop" de
### Monolix esta deteniendose por estabilidad, no por falta de presupuesto, y que
### el problema es de identificabilidad (cresta plana), no de convergencia
### insuficiente.
###
### Casos elegidos (valores iniciales EXACTOS de esas corridas, ver CSV):
###  - run 8  ("clasificacion degenerada"): partida original -> alpha=80881.5, beta=-18317.8, OFV=661.40
###  - run 10 ("falla total"): partida original -> OFV = Inf
###  - run 11 ("diverge en todo"): partida original -> OFV=873.98 (vs ~657 del optimo bueno)

library(lixoftConnectors)
initializeLixoftConnectors(software = "monolix", force = TRUE)

setwd(file.path(REPO, "monolix", "data1", "Model_1"))   # REPO = path to your clone (define it once; see README)
project <- "Model.1.mlxtran"

casos <- list(
  caso8_clasificacion_degenerada = c(a_pop = 4.62727826514281, b_pop = 11.4877842290327,
    c_pop = 10.6264668183122, alpha_pop = 38.8451444911771, beta_pop = -11.5786571766809,
    omega_a = 0.695949697401375, omega_b = 3.38673825793201, omega_c = 2.37103548687883),
  caso10_falla_total = c(a_pop = 13.2656067795442, b_pop = 53.7520772833005,
    c_pop = 5.96658245291654, alpha_pop = -88.066027803903, beta_pop = 23.0311746217008,
    omega_a = 1.11094896934, omega_b = 2.97363493212289, omega_c = 7.03941625866224),
  caso11_diverge_todo = c(a_pop = 2.64218413194083, b_pop = 59.4521162133571,
    c_pop = 22.7335345181799, alpha_pop = -132.134928007515, beta_pop = 4.33192519332795,
    omega_a = 1.96405157812629, omega_b = 8.34064079754474, omega_c = 7.24808435458574)
)

## Resultados YA CONOCIDOS con el techo original (1x), para comparar sin re-correrlos:
resultado_1x <- data.frame(
  caso = names(casos),
  OFV_1x = c(661.404567609234, Inf, 873.978037332877),
  alpha_pop_1x = c(80881.5256472564, -88.066027803903, 125382.857704503),
  beta_pop_1x  = c(-18317.7787648133, 23.0311746217008, -15930.0947406202)
)

resultados <- data.frame()

for (nombre in names(casos)) {
  cat("=== Caso:", nombre, "(techo x10) ===\n")
  ini <- casos[[nombre]]
  res_k <- tryCatch({
    loadProject(project)

    popparams <- getPopulationParameterInformation()
    for (nm in names(ini)) {
      popparams$initialValue[popparams$name == nm] <- as.numeric(ini[nm])
    }
    setPopulationParameterInformation(popparams)

    scenario <- getScenario()
    scenario$tasks <- c(populationParameterEstimation = TRUE,
                         conditionalModeEstimation = FALSE,
                         conditionalDistributionSampling = FALSE,
                         standardErrorEstimation = FALSE,
                         logLikelihoodEstimation = TRUE)
    setScenario(scenario)

    ## --- subir el techo de iteraciones x10 (y, si se puede, apagar autostop) ---
    pset <- getPopulationParameterEstimationSettings()
    cat("Campos disponibles en la configuracion de SAEM:\n")
    print(names(pset))
    campos_tocados <- character(0)
    for (campo in c("nbexploratoryiterations", "nbsmoothingiterations")) {
      if (campo %in% names(pset)) {
        pset[[campo]] <- pset[[campo]] * 10
        campos_tocados <- c(campos_tocados, campo)
      }
    }
    for (campo in c("exploratoryautostop", "smoothingautostop")) {
      if (campo %in% names(pset)) pset[[campo]] <- FALSE
    }
    setPopulationParameterEstimationSettings(pset)
    cat("Campos de iteraciones que SI se multiplicaron x10:", paste(campos_tocados, collapse=", "), "\n")
    pset_check <- getPopulationParameterEstimationSettings()
    cat("Valores aplicados luego de setPopulationParameterEstimationSettings:\n")
    print(pset_check[campos_tocados])

    t0 <- Sys.time()
    runScenario()
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    est <- getEstimatedPopulationParameters()
    ll  <- getEstimatedLogLikelihood()
    ll_flat <- unlist(ll)
    ofv_val <- ll_flat[grepl("OFV$", names(ll_flat))]
    ofv_val <- if (length(ofv_val) > 0) as.numeric(ofv_val[1]) else NA_real_

    data.frame(caso = nombre, status = "OK",
               a_pop_10x = as.numeric(est["a_pop"]), b_pop_10x = as.numeric(est["b_pop"]),
               c_pop_10x = as.numeric(est["c_pop"]), alpha_pop_10x = as.numeric(est["alpha_pop"]),
               beta_pop_10x = as.numeric(est["beta_pop"]),
               OFV_10x = ofv_val, elapsed_sec = elapsed,
               campos_iteraciones_multiplicados = paste(campos_tocados, collapse = ";"))
  }, error = function(e) {
    data.frame(caso = nombre, status = paste("ERROR:", conditionMessage(e)))
  })
  resultados <- dplyr::bind_rows(resultados, res_k)
}

comparacion <- merge(resultado_1x, resultados, by = "caso")
write.csv(comparacion, "prueba_mas_iteraciones_M1_comparacion.csv", row.names = FALSE)
cat("\n=== Comparacion techo original (1x) vs techo x10 ===\n")
print(comparacion)

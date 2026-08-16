### Identificabilidad de los coeficientes de clasificacion en funcion del tamano
### muestral N -- Modelo I (embarazo). Complementa la Seccion 3.2/3.3 del paper.
###
### Mensaje que se busca demostrar: la debil identificabilidad de (eta0, eta1) NO es
### un defecto del metodo ni del algoritmo, sino un reflejo de la INFORMACION en los
### datos. Al aumentar N (y por tanto el numero de eventos), la fraccion de ajustes
### que divergen -> 0 y la pendiente eta1 se vuelve bien identificada, mientras que la
### frontera de decision a* = -eta0/eta1 se mantiene estable en TODO el rango.
###
### Diseno: se generan datasets sinteticos desde la "verdad" del Modelo I, remuestreando
### (con reemplazo) los patrones de visitas de los 173 sujetos reales para construir N
### sujetos. Para cada N y cada replica se reajusta con SAEM (partida razonable +/-30%)
### y se guardan eta0, eta1 y a*.

## Set FIT <- FALSE to skip the (long) Monolix sweep and only rebuild the summary and the
## figure from an existing identificabilidad_vs_N_resultados.csv.
if (!exists("FIT")) FIT <- TRUE

library(ggplot2)
if (FIT) {
  library(lixoftConnectors)
  initializeLixoftConnectors(software = "monolix", force = TRUE)
}

setwd(file.path(REPO, "monolix", "data1", "Model_1"))   # REPO = path to your clone (define it once; see README)
project <- "Model.1.mlxtran"

## --- pool de disenos reales (tiempos de visita + fila binaria + group) ---
real <- read.csv("jointmodel_data_new.csv", colClasses = c(group = "character"))
ids <- sort(unique(real$id))
design_pool <- lapply(ids, function(i) {
  di <- real[real$id == i, ]
  list(times_long = di$time[di$ytype == 1],
       time_bin   = di$time[di$ytype == 2],
       group      = di$group[1])
})

## --- verdad (Modelo I, de Model.1/summary.txt) ---
truth <- c(a_pop = 4.54, b_pop = 15.6, c_pop = 7.22, alpha_pop = 34.3, beta_pop = -7.94,
           omega_a = 0.502, omega_b = 3.31, omega_c = 1.95)
sigma_resid <- 0.261

simular_N <- function(truth, sigma_resid, pool, N, seed) {
  set.seed(seed)
  picks <- sample(seq_along(pool), N, replace = TRUE)
  a_i <- rnorm(N, truth["a_pop"], truth["omega_a"])
  b_i <- rnorm(N, truth["b_pop"], truth["omega_b"])
  c_i <- pmax(rnorm(N, truth["c_pop"], truth["omega_c"]), 0.1)
  filas <- vector("list", N)
  for (k in seq_len(N)) {
    d <- pool[[picks[k]]]
    t_long <- d$times_long
    E <- a_i[k] / (1 + exp(-(t_long - b_i[k]) / c_i[k]))
    y_long <- E + rnorm(length(t_long), 0, sigma_resid)
    p1 <- 1 / (1 + exp(-(truth["alpha_pop"] + truth["beta_pop"] * a_i[k])))
    level <- rbinom(1, 1, p1)
    filas[[k]] <- rbind(
      data.frame(id = k, time = t_long,     y = y_long, group = d$group, ytype = 1),
      data.frame(id = k, time = d$time_bin, y = level,  group = d$group, ytype = 2))
  }
  do.call(rbind, filas)
}

## --- barrido sobre N ---
## Corrida completa, 4 valores de N a R=50 (consistente con Secciones 3.1 y 3.5).
Ns <- c(173, 350, 700, 1400)
R  <- 50
if (FIT) {
resultados <- data.frame()

for (N in Ns) {
  for (r in seq_len(R)) {
    cat(sprintf("=== N=%d  rep=%d/%d ===\n", N, r, R))
    res <- tryCatch({
      sim <- simular_N(truth, sigma_resid, design_pool, N, seed = N * 1000 + r)
      write.csv(sim, "simN.csv", row.names = FALSE, quote = c(4))  # quote solo group
      loadProject(project)
      BaseData <- getData()
      setData("simN.csv", headerTypes = BaseData$headerTypes,
              observationTypes = BaseData$observationTypes)
      pp <- getPopulationParameterInformation()
      ini <- truth * runif(length(truth), 0.7, 1.3)
      for (nm in names(ini)) pp$initialValue[pp$name == nm] <- as.numeric(ini[nm])
      setPopulationParameterInformation(pp)
      sc <- getScenario()
      sc$tasks <- c(populationParameterEstimation = TRUE, conditionalModeEstimation = FALSE,
                    conditionalDistributionSampling = FALSE, standardErrorEstimation = FALSE,
                    logLikelihoodEstimation = FALSE)   # solo estimacion puntual -> mas rapido
      setScenario(sc)
      t0 <- Sys.time()
      runScenario()
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      est <- getEstimatedPopulationParameters()
      ea <- as.numeric(est["alpha_pop"]); eb <- as.numeric(est["beta_pop"])
      data.frame(N = N, rep = r, status = "OK",
                 est_alpha = ea, est_beta = eb, a_star = -ea / eb,
                 est_a_pop = as.numeric(est["a_pop"]), elapsed_sec = elapsed)
    }, error = function(e) data.frame(N = N, rep = r, status = paste("ERROR:", conditionMessage(e))))
    resultados <- dplyr::bind_rows(resultados, res)
    unlink("simN.csv")
  }
}

write.csv(resultados, "identificabilidad_vs_N_resultados.csv", row.names = FALSE)
cat("\nSweep done. Saved to identificabilidad_vs_N_resultados.csv\n")
}  # end if (FIT)

## Read back from disk, so the summary/figure below are identical whether or not the sweep
## was just run in this session.
resultados <- read.csv("identificabilidad_vs_N_resultados.csv")

## --- resumen por N ---
ok <- resultados[resultados$status == "OK", ]
cv <- function(x) 100 * sd(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE))
resumen <- do.call(rbind, lapply(sort(unique(ok$N)), function(N) {
  s <- ok[ok$N == N, ]
  div <- abs(s$est_alpha) > 200            # ajuste "divergente" (pendiente disparada)
  nodiv <- s[!div, ]
  data.frame(N = N, n = nrow(s),
             frac_divergente_pct = round(100 * mean(div), 1),
             CV_a_star_pct       = round(cv(s$a_star), 1),        # deberia ser bajo en todo N
             CV_beta_pct_todos   = round(cv(s$est_beta), 1),      # dominado por divergentes
             CV_beta_pct_noDiv   = round(cv(nodiv$est_beta), 1),  # se encoge al crecer N
             media_beta_noDiv    = round(mean(nodiv$est_beta), 2))
}))
cat("\n--- Identificabilidad en funcion de N (verdad: eta1=-7.94, a*=4.32) ---\n")
print(resumen, row.names = FALSE)
write.csv(resumen, "identificabilidad_vs_N_resumen.csv", row.names = FALSE)

## --- Figure (Section 3.4): divergent fits and CV of a*, as a function of N ---
out_fig <- file.path(REPO, "figures")
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)
d <- resultados[resultados$status == "OK", ]
res <- do.call(rbind, lapply(sort(unique(d$N)), function(N) {
  s <- d[d$N == N, ]; div <- abs(s$est_alpha) > 200; nd <- s[!div, ]
  data.frame(N = N, divergent = round(100 * mean(div), 1),
             cv_astar = round(cv(s$a_star), 1),
             cv_beta = round(cv(nd$est_beta), 1),
             mean_astar = round(mean(s$a_star), 2))
}))
write.csv(res, "identificabilidad_vs_N_resumen_4puntos.csv", row.names = FALSE)
print(res, row.names = FALSE)

long <- data.frame(N = rep(res$N, 2),
  serie = factor(rep(c("Divergent fits", "CV of decision boundary a*"), each = nrow(res)),
                 levels = c("Divergent fits", "CV of decision boundary a*")),
  value = c(res$divergent, res$cv_astar))
p <- ggplot(long, aes(N, value, colour = serie, shape = serie)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  scale_x_continuous(trans = "log2", breaks = res$N) +
  scale_colour_manual(values = c("Divergent fits" = "#d73027",
                                 "CV of decision boundary a*" = "#4575b4"), name = NULL) +
  scale_shape_manual(values = c(16, 15), name = NULL) +
  labs(x = "Sample size  N  (log scale)", y = "Percentage (%)") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        plot.background = element_rect(fill = "white", colour = NA))
ggsave(file.path(out_fig, "identif_vs_N.pdf"), p, width = 7, height = 4.4)
cat("Figure written: figures/identif_vs_N.pdf\n")

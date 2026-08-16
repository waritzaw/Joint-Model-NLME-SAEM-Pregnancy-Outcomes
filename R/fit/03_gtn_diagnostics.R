### Diagnosticos OP y VPC para la aplicacion GTN (data2) -- puntos 12 y 13 del referato.
###
### 100% reproducible y SIN conector de Monolix: los OP se leen de predictions.txt
### (ya calculado por Monolix) y los VPC se simulan desde populationParameters.txt.
### Modelo longitudinal GTN: E = a_i + b_i * t, con (a_i,b_i) ~ N normal bivariada,
### error REM I (constant) o REM IV (combined2). Estilo grafico UNIFORME (punto 13).
###
### Incluye ademas la curva ROC / AUC con IC de DeLong (antes roc_GTN.R).
###
### Salidas (a figures/): OPdata2_REMI.pdf, VPCdata2_REMI.pdf,
###                       OPdata2_REMIV.pdf, VPCdata2_REMIV.pdf, ROC-data2.pdf

library(ggplot2)

base_dir <- file.path(REPO, "monolix", "data2")
out_dir  <- file.path(REPO, "figures")

## ---- Tema uniforme para TODAS las figuras (punto 13) ----
theme_paper <- theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = "white", colour = NA),
        axis.title = element_text(size = 15),
        axis.text  = element_text(size = 13),
        legend.position = "top")

models <- list(
  REMI  = list(dir = "Modelo-lineal",                error = "constant"),
  REMIV = list(dir = "Modelo-lineal-errorCombined2", error = "combined2")
)

read_pop <- function(dir) {
  p <- read.csv(file.path(base_dir, dir, "populationParameters.txt"))
  setNames(p$value, p$parameter)
}

## =====================================================================
## 1) OP: observaciones vs predicciones individuales (de predictions.txt)
## =====================================================================
for (nm in names(models)) {
  m <- models[[nm]]
  pred <- read.csv(file.path(base_dir, m$dir, "predictions.txt"))
  wres <- pred$indWRes_mode
  is_out <- abs(wres) > qnorm(0.95)          # fuera del intervalo de prediccion 90%
  outlier_prop <- mean(is_out, na.rm = TRUE)
  df <- data.frame(obs = pred$y1, ip = pred$indivPred_mode, out = is_out)
  lim <- range(c(df$obs, df$ip), na.rm = TRUE)

  p <- ggplot(df, aes(ip, obs)) +
    geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.7) +
    geom_point(aes(colour = out), size = 1.4, alpha = 0.7) +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "orange"),
                        labels = c("within 90% PI", "outlier"), name = NULL) +
    coord_equal(xlim = lim, ylim = lim) +
    labs(x = "Individual predictions", y = "Observations",
         title = sprintf("Observations vs predictions (%s)", sub("REM", "REM ", nm)),
         subtitle = sprintf("Outliers proportion: %.1f%%", 100 * outlier_prop)) +
    theme_paper
  ggsave(file.path(out_dir, sprintf("OPdata2_%s.pdf", nm)), p, width = 6, height = 5.5)
  cat(sprintf("OP %s: outliers = %.1f%%\n", nm, 100 * outlier_prop))
}

## =====================================================================
## 2) VPC por simulacion desde el modelo ajustado
## =====================================================================
dat <- read.csv(file.path(base_dir, "dataHCG.csv"))
dat <- dat[dat$ytype == 1, c("id", "y", "time")]   # solo longitudinal
ids <- unique(dat$id)
row_id <- match(dat$id, ids)                        # fila -> indice de sujeto

simular_y <- function(par, error, Nsim = 500, seed = 2026) {
  set.seed(seed)
  mu <- c(par["a_pop"], par["b_pop"])
  Sig <- matrix(c(par["omega_a"]^2,
                  par["corr_b_a"] * par["omega_a"] * par["omega_b"],
                  par["corr_b_a"] * par["omega_a"] * par["omega_b"],
                  par["omega_b"]^2), 2, 2)
  R <- chol(Sig)                                    # z %*% R ~ N(0, Sig)
  nId <- length(ids)
  simmat <- matrix(NA_real_, nrow(dat), Nsim)
  for (s in seq_len(Nsim)) {
    AB <- matrix(mu, nId, 2, byrow = TRUE) + matrix(rnorm(nId * 2), nId, 2) %*% R
    a_i <- AB[row_id, 1]; b_i <- AB[row_id, 2]
    E <- a_i + b_i * dat$time
    sdev <- if (error == "constant") par["a_"] else sqrt(par["a_"]^2 + par["b_"]^2 * E^2)
    simmat[, s] <- E + rnorm(nrow(dat), 0, sdev)
  }
  simmat
}

vpc_bands <- function(simmat, probs = c(0.1, 0.5, 0.9)) {
  times <- sort(unique(dat$time))
  res <- list()
  for (pr in probs) {
    band <- sapply(times, function(tt) {
      rows <- which(dat$time == tt)
      obs_p <- quantile(dat$y[rows], pr, names = FALSE)
      sim_p <- apply(simmat[rows, , drop = FALSE], 2, quantile, probs = pr, names = FALSE)
      c(obs = obs_p, lo = quantile(sim_p, 0.05, names = FALSE),
        hi = quantile(sim_p, 0.95, names = FALSE))
    })
    res[[as.character(pr)]] <- data.frame(time = times, prob = factor(pr),
                                          obs = band["obs", ],
                                          lo = band["lo", ], hi = band["hi", ])
  }
  do.call(rbind, res)
}

for (nm in names(models)) {
  m <- models[[nm]]
  par <- read_pop(m$dir)
  sm <- simular_y(par, m$error)
  vb <- vpc_bands(sm)

  p <- ggplot(vb, aes(time, group = prob)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = prob), alpha = 0.3) +
    geom_line(aes(y = obs, colour = prob), linewidth = 0.8) +
    geom_point(aes(y = obs, colour = prob), size = 1.6) +
    scale_fill_manual(values = c("0.1" = "#4575b4", "0.5" = "#d73027", "0.9" = "#4575b4"),
                      name = "Percentile") +
    scale_colour_manual(values = c("0.1" = "#4575b4", "0.5" = "#d73027", "0.9" = "#4575b4"),
                        name = "Percentile") +
    labs(x = "Time (weeks since evacuation)",
         y = expression(log[10](beta*"-HCG")),
         title = sprintf("Visual predictive check (%s)", sub("REM", "REM ", nm))) +
    theme_paper
  ggsave(file.path(out_dir, sprintf("VPCdata2_%s.pdf", nm)), p, width = 6.5, height = 5)
  cat(sprintf("VPC %s: exportado\n", nm))
}

cat("\nListo. 4 PDFs escritos en:", out_dir, "\n")


## ================= ROC curve / AUC with DeLong CI (formerly roc_GTN.R) =================
library(pROC)

## edad por id (constante por paciente)
dat <- read.csv(file.path(base_dir, "dataHCG.csv"))
age_by_id <- tapply(dat$age, dat$id, function(x) x[1])

prob_gtn <- function(dir) {
  ip <- read.csv(file.path(base_dir, dir, "IndividualParameters",
                           "estimatedIndividualParameters.txt"))
  ## coeficientes logisticos (constantes entre individuos): del primer registro
  a0 <- ip$alpha_mode[1]; b1 <- ip$beta_mode[1]
  g1 <- ip$gamma_mode[1]; d1 <- ip$delta_mode[1]
  age <- as.numeric(age_by_id[as.character(ip$id)])
  lp <- a0 + b1 * ip$a_mode + g1 * ip$b_mode + d1 * age
  data.frame(id = ip$id, prob = plogis(lp), truth = ip$group)
}

pr1 <- prob_gtn("Modelo-lineal")                 # REM I
pr4 <- prob_gtn("Modelo-lineal-errorCombined2")  # REM IV

roc1 <- roc(pr1$truth, pr1$prob, quiet = TRUE, direction = "<")
roc4 <- roc(pr4$truth, pr4$prob, quiet = TRUE, direction = "<")
ci1 <- as.numeric(ci.auc(roc1, method = "delong"))  # [low, auc, high]
ci4 <- as.numeric(ci.auc(roc4, method = "delong"))

cat(sprintf("REM I : AUC=%.4f  95%% CI [%.4f, %.4f]\n", ci1[2], ci1[1], ci1[3]))
cat(sprintf("REM IV: AUC=%.4f  95%% CI [%.4f, %.4f]\n", ci4[2], ci4[1], ci4[3]))

pdf(file.path(out_dir, "ROC-data2.pdf"), width = 6.5, height = 6.5)
par(mar = c(4.5, 4.5, 1, 1), cex.lab = 1.2, cex.axis = 1.05)
plot(roc4, col = "black", lty = 1, lwd = 2, legacy.axes = TRUE,
     xlab = "False positive rate", ylab = "True positive rate")
plot(roc1, col = "red", lty = 2, lwd = 2, add = TRUE)
legend("bottomright", bty = "n", lwd = 2, col = c("red", "black"), lty = c(2, 1),
       legend = c(sprintf("Joint Model REM I:  AUC = %.3f (95%% CI %.3f-%.3f)", ci1[2], ci1[1], ci1[3]),
                  sprintf("Joint Model REM IV: AUC = %.3f (95%% CI %.3f-%.3f)", ci4[2], ci4[1], ci4[3])))
dev.off()

cat("\nListo. ROC-data2.pdf regenerado en", out_dir, "\n")


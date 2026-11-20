# simular_datos.R ---------------------------------------------------------
# Genera un dataset FICTICIO con la misma estructura que el del estudio.
# Los datos reales no se publican: son de pacientes.
#
# Qué reproduce por construcción:
#   - n = 807, con un 77 % de mujeres y edad mediana en torno a 53 años.
#   - La tabla 2x2 del artículo (27 / 43 / 33 / 704) con los puntos de corte
#     del fabricante (CLIA > 29 UI/mL; RIA > 7 UI/mL).
#   - Los límites de cuantificación: RIA entre 2,5 y 100 UI/mL; CLIA hasta 600.
# Qué NO reproduce exactamente: el resto de cifras (correlación, Bland-Altman,
# AUC, Tobit) saldrán parecidas, pero no idénticas a las publicadas.
#
# Uso: source("R/simular_datos.R") crea data/antidsdna_simulado.csv
# -------------------------------------------------------------------------

set.seed(2019)

n_both      <- 27   # CLIA+ y RIA+
n_clia_only <- 43   # CLIA+ y RIA-
n_ria_only  <- 33   # CLIA- y RIA+
n_neither   <- 704  # CLIA- y RIA-
n <- n_both + n_clia_only + n_ria_only + n_neither

# Límites analíticos y puntos de corte (inserto de cada fabricante)
RIA_LOQ  <- 2.5;  RIA_MAX  <- 100; RIA_CUT  <- 7
CLIA_MAX <- 600;  CLIA_CUT <- 29

# Lognormal truncada entre dos límites (por inversión de la función de distribución)
rlnorm_trunc <- function(n, meanlog, sdlog, lower, upper) {
  p <- runif(n, plnorm(lower, meanlog, sdlog), plnorm(upper, meanlog, sdlog))
  qlnorm(p, meanlog, sdlog)
}

# RIA negativo: la mayoría por debajo del límite de cuantificación (se informa 2,5)
sim_ria_neg <- function(k) {
  below_loq <- runif(k) < 0.67
  ifelse(below_loq, RIA_LOQ, rlnorm_trunc(k, log(3.5), 0.4, RIA_LOQ, RIA_CUT))
}
# RIA positivo: cola larga; lo que pasa de 100 se informa como 100
sim_ria_pos <- function(k) pmin(rlnorm_trunc(k, log(14), 1.3, RIA_CUT, 1e4), RIA_MAX)

# CLIA negativo: muchos ceros y una cola hasta el punto de corte.
# n_grey casos caen en la "zona gris" (19,2-29 UI/mL), justo por debajo del corte.
sim_clia_neg <- function(k, n_grey) {
  zero <- runif(k - n_grey) < 0.33
  low  <- ifelse(zero, 0, rlnorm_trunc(k - n_grey, log(3), 1.1, 0.1, 19.1))
  grey <- runif(n_grey, 19.2, CLIA_CUT)
  sample(round(c(low, grey), 1))
}
# CLIA positivo: cola muy larga; lo que pasa de 600 se informa como 600
sim_clia_pos <- function(k, meanlog) pmin(round(rlnorm_trunc(k, meanlog, 1.1, CLIA_CUT + 0.1, 1e5), 1), CLIA_MAX)

group <- rep(c("both", "clia_only", "ria_only", "neither"),
             times = c(n_both, n_clia_only, n_ria_only, n_neither))

ria  <- numeric(n)
clia <- numeric(n)
ria[group == "both"]       <- sim_ria_pos(n_both)
ria[group == "ria_only"]   <- sim_ria_pos(n_ria_only)
ria[group == "clia_only"]  <- sim_ria_neg(n_clia_only)
ria[group == "neither"]    <- sim_ria_neg(n_neither)
clia[group == "both"]      <- sim_clia_pos(n_both, log(110))
clia[group == "clia_only"] <- sim_clia_pos(n_clia_only, log(50))
clia[group == "ria_only"]  <- sim_clia_neg(n_ria_only, n_grey = 7)
clia[group == "neither"]   <- sim_clia_neg(n_neither, n_grey = 61)

antidsdna <- data.frame(
  id   = sprintf("S%03d", seq_len(n)),
  sex  = sample(c("F", "M"), n, replace = TRUE, prob = c(0.77, 0.23)),
  age  = pmin(pmax(round(rnorm(n, 53, 17)), 15), 95),
  ria  = round(ria, 1),
  clia = clia
)
# Orden aleatorio, para que el grupo no se deduzca de la posición
antidsdna <- antidsdna[sample(n), ]
antidsdna$id <- sprintf("S%03d", seq_len(n))

dir.create("data", showWarnings = FALSE)
write.csv(antidsdna, "data/antidsdna_simulado.csv", row.names = FALSE)

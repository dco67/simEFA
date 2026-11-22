# ===========================================================================
# _common.R — réglages et objets partagés par les chapitres du manuel.
# Chaque chapitre le charge au début : source("_common.R").
# Les données sont simulées EXACTEMENT comme dans simEFA (mêmes fonctions,
# copiées dans _moteur.R) : les nombres cités sont ceux de l'écran.
# ===========================================================================
suppressPackageStartupMessages({ library(ggplot2); library(lavaan); library(psych) })
source("_moteur.R")

knitr::opts_chunk$set(
  echo = FALSE, message = FALSE, warning = FALSE,
  fig.align = "center", out.width = "88%", fig.width = 6.4, fig.height = 3.6,
  dev = "ragg_png", dpi = if (knitr::is_latex_output()) 300 else 150
)

# --- Nombres à la française -------------------------------------------------
fr <- function(x, d = 2) {
  x <- as.numeric(x); x[!is.na(x) & abs(x) < 0.5 * 10^(-d)] <- 0
  sub("^-", "−", formatC(x, format = "f", digits = d, decimal.mark = ",", big.mark = " "))
}
pc <- function(x, d = 1) paste0(fr(100 * x, d), " %")
frm <- function(x, d = 2) gsub(",", "{,}", sub("−", "-", fr(x, d), fixed = TRUE), fixed = TRUE)
pv <- function(p) if (p < 0.001) "< 0,001" else paste0("= ", fr(p, 3))
axe_fr <- function(x) sub("^-", "−", format(x, decimal.mark = ",", big.mark = " ",
                                                   trim = TRUE, scientific = FALSE, drop0trailing = TRUE))

# --- Couleurs et thème ---------------------------------------------------------
pal <- list(bleu = "#2c7fb8", orange = "#d95f0e", vert = "#1b9e77", rouge = "#c0392b",
            gris = "grey45", violet = "#7570b3", or = "#d9a53b")
theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(), legend.position = "bottom",
                  strip.text = element_text(face = "bold"),
                  plot.title = element_text(face = "bold", size = 12)))

# --- Les quatre situations du laboratoire (préréglages de simEFA) ------------
SCENARIOS <- list(
  clear      = list(n = 300, nf = 2, nv = 4, main = .70, cross = .05, prop = .25, noise = 1, corr = .30),
  correlated = list(n = 300, nf = 2, nv = 4, main = .65, cross = .05, prop = .25, noise = 1, corr = .55),
  cross      = list(n = 300, nf = 2, nv = 4, main = .60, cross = .30, prop = .50, noise = 1, corr = .30),
  fragile    = list(n = 100, nf = 3, nv = 3, main = .50, cross = .15, prop = .35, noise = 1, corr = .30)
)
simuler <- function(scenario = "clear", graine = 1234, ...) {
  p <- modifyList(SCENARIOS[[scenario]], list(...))
  generate_efa_simulation(p$n, p$nf, p$nv, p$main, p$cross, p$prop, p$noise, TRUE, p$corr, graine)
}
pop_de <- function(d, scenario = "clear", ...) {
  p <- modifyList(SCENARIOS[[scenario]], list(...))
  phi <- matrix(p$corr, p$nf, p$nf); diag(phi) <- 1
  population_simulee(attr(d, "loading_pattern"), phi, p$noise)
}
parallele <- function(d, fm = "ml", graine = 2026) {
  set.seed(graine)
  utils::capture.output(pa <- suppressWarnings(suppressMessages(psych::fa.parallel(d, fa = "both", fm = fm, n.iter = 100, plot = FALSE))))
  list(fa_obs = pa$fa.values, fa_sim = pa$fa.sim, pc_obs = pa$pc.values, pc_sim = pa$pc.sim,
       sugg_fa = suggestion_parallele(pa$fa.values, pa$fa.sim),
       sugg_pc = suggestion_parallele(pa$pc.values, pa$pc.sim), kaiser = sum(pa$pc.values > 1))
}

# --- Le cas fil rouge : « Structure claire », graine 1234 -----------------
fil <- simuler("clear")
fil_pop <- pop_de(fil, "clear")
fil_afe <- suppressWarnings(ajuster_afe(fil, 2, "oblimin", "ML"))

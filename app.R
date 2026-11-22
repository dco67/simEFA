library(shiny)
library(bslib)
library(DT)
library(plotly)
library(shinyjs)
library(readxl)
library(writexl)
library(haven)
library(stats)
library(psych) 
library(corrplot)
library(lavaan)
library(ggplot2)

# Les tableaux DT prennent leur hauteur naturelle (pas de hauteur fixe de 400 px
# qui tronquerait les tableaux et leur légende à l'intérieur des cartes).
DTOutput <- function(outputId, width = "100%", height = "auto") {
  DT::DTOutput(outputId, width = width, height = height, fill = FALSE)
}
# Les jeux classiques d'EFA.dimensions et d'EFAtools sont lus à la demande
# avec data(); ces paquets n'ont donc pas besoin d'être chargés au démarrage.

# Fonction pour convertir un modèle EFA en syntaxe CFA
convert_efa_to_cfa <- function(efa_result, loading_threshold = 0.3,
                               crossloading_rule = c("primary", "all")) {
  crossloading_rule <- match.arg(crossloading_rule)
  # Extraire les charges factorielles 
  loadings_matrix <- efa_result$loadings
  
  # Identifier les facteurs
  factors <- colnames(loadings_matrix)
  
  # Initialiser la syntaxe CFA
  cfa_syntax <- c()
  
  primary_factor <- apply(abs(loadings_matrix), 1, which.max)
  primary_value <- apply(abs(loadings_matrix), 1, max)

  # Par défaut, la CFA candidate retient seulement l'appartenance principale.
  # Une charge croisée peut être conservée explicitement, mais elle ne doit pas
  # s'ajouter silencieusement au modèle sous prétexte qu'elle dépasse un seuil.
  for (factor_idx in 1:length(factors)) {
    factor <- factors[factor_idx]
    factor_items <- if (crossloading_rule == "all") {
      rownames(loadings_matrix)[abs(loadings_matrix[, factor_idx]) >= loading_threshold]
    } else {
      rownames(loadings_matrix)[primary_factor == factor_idx & primary_value >= loading_threshold]
    }
    
    if (length(factor_items) > 0) {
      # Construire la syntaxe du modèle CFA pour ce facteur
      cfa_syntax <- c(cfa_syntax, paste(factor, "=~", paste(factor_items, collapse=" + ")))
    }
  }
  
  # Combiner en une seule chaîne
  cfa_model <- paste(cfa_syntax, collapse = "\n")
  
  return(cfa_model)
}

# Extrait les charges entièrement standardisées (Std.all) du modèle lavaan.
# Les coefficients bruts dépendent de l'échelle des variables et ne doivent
# pas être comparés directement à un seuil de lecture comme .30 ou .40.
extract_standardized_efa_loadings <- function(efa_fit, item_names = NULL) {
  fit <- if (inherits(efa_fit, "efaList")) efa_fit[[1]] else efa_fit
  std <- tryCatch(lavaan::standardizedSolution(fit), error = function(e) NULL)
  if (is.null(std) || !"est.std" %in% names(std)) return(NULL)

  std <- std[std$op == "=~" & !is.na(std$est.std), , drop = FALSE]
  if (nrow(std) == 0) return(NULL)

  if (is.null(item_names)) item_names <- unique(std$rhs)
  std <- std[std$rhs %in% item_names, , drop = FALSE]
  factor_names <- unique(std$lhs)
  if (length(factor_names) == 0) return(NULL)

  loadings <- matrix(0, nrow = length(item_names), ncol = length(factor_names),
                     dimnames = list(item_names, paste0("F", seq_along(factor_names))))
  for (i in seq_len(nrow(std))) {
    row <- match(std$rhs[i], item_names)
    col <- match(std$lhs[i], factor_names)
    if (!is.na(row) && !is.na(col)) loadings[row, col] <- std$est.std[i]
  }
  residuals <- tryCatch(lavaan::standardizedSolution(fit), error = function(e) NULL)
  if (!is.null(residuals)) {
    residuals <- residuals[residuals$op == "~~" & residuals$lhs == residuals$rhs &
                             residuals$lhs %in% item_names, , drop = FALSE]
    communalities <- rep(NA_real_, length(item_names))
    names(communalities) <- item_names
    communalities[match(residuals$lhs, item_names)] <- 1 - residuals$est.std
    attr(loadings, "communalities") <- communalities
  }
  loadings
}

# Génère un nouvel échantillon issu de la même population simulée. La graine
# change, mais tous les paramètres générateurs restent identiques.
generate_efa_simulation <- function(n, num_factors, vars_per_factor,
                                    loading_strength, cross_loading,
                                    cross_proportion, error_multiplier,
                                    correlated, factor_cor, seed,
                                    loading_pattern = NULL) {
  factor_cov <- if (correlated) {
    tmp <- matrix(factor_cor, nrow = num_factors, ncol = num_factors)
    diag(tmp) <- 1
    tmp
  } else {
    diag(num_factors)
  }
  set.seed(seed)
  total_vars <- num_factors * vars_per_factor
  if (is.null(loading_pattern)) {
    loading_pattern <- matrix(0, nrow = total_vars, ncol = num_factors)
    for (f in seq_len(num_factors)) {
      variable_indices <- ((f - 1) * vars_per_factor + 1):(f * vars_per_factor)
      loading_pattern[variable_indices, f] <- loading_strength
      if (num_factors > 1 && cross_loading > 0 && cross_proportion > 0) {
        ambiguous <- variable_indices[runif(length(variable_indices)) < cross_proportion]
        for (var_idx in ambiguous) {
          # Attention : sample(x, 1) avec un x de longueur 1 tire dans 1:x.
          # On choisit donc explicitement parmi les AUTRES facteurs, sans
          # jamais écraser la charge principale de l'item.
          autres <- setdiff(seq_len(num_factors), f)
          other <- if (length(autres) == 1) autres else sample(autres, size = 1)
          loading_pattern[var_idx, other] <- cross_loading
        }
      }
    }
  }
  if (!is.matrix(loading_pattern) || nrow(loading_pattern) != total_vars ||
      ncol(loading_pattern) != num_factors) {
    stop("Le patron de charges simulées est incompatible avec les réglages demandés.")
  }
  common_variance <- diag(loading_pattern %*% factor_cov %*% t(loading_pattern))
  if (any(!is.finite(common_variance)) || any(common_variance >= 0.98)) {
    stop("Cette combinaison de charges et de corrélation est trop forte.")
  }

  if (correlated) {
    factors <- matrix(rnorm(n * num_factors), nrow = n, ncol = num_factors) %*%
      chol(factor_cov)
  } else {
    factors <- matrix(rnorm(n * num_factors), nrow = n, ncol = num_factors)
  }
  factors <- scale(factors)

  data <- matrix(0, nrow = n, ncol = total_vars)
  for (var_idx in seq_len(total_vars)) {
    signal <- factors %*% loading_pattern[var_idx, ]
    noise_sd <- sqrt(1 - common_variance[var_idx]) * error_multiplier
    data[, var_idx] <- scale(signal + rnorm(n, 0, noise_sd))
  }
  data <- as.data.frame(data)
  names(data) <- paste0("X", seq_len(total_vars))
  attr(data, "common_variance") <- common_variance
  attr(data, "loading_pattern") <- loading_pattern
  data
}

# Fonction pour générer des interprétations automatiques des facteurs
generate_factor_interpretation <- function(loadings, threshold = 0.3) {
  # Obtenir le nombre de facteurs
  n_factors <- ncol(loadings)
  
  # Structure pour stocker les interprétations
  interpretations <- list()
  
  # Pour chaque facteur
  for (f in 1:n_factors) {
    # Trouver les variables qui chargent significativement sur ce facteur
    significant_loadings <- loadings[, f]
    
    # Trier les charges par ordre décroissant de valeur absolue
    sorted_idx <- order(abs(significant_loadings), decreasing = TRUE)
    sorted_loadings <- significant_loadings[sorted_idx]
    var_names <- rownames(loadings)[sorted_idx]
    
    # Filtrer les charges au-dessus du seuil
    above_threshold <- abs(sorted_loadings) >= threshold
    strong_vars <- var_names[above_threshold]
    strong_loadings <- sorted_loadings[above_threshold]
    
    # Séparer les charges positives et négatives
    positive_vars <- strong_vars[strong_loadings > 0]
    negative_vars <- strong_vars[strong_loadings < 0]
    
    # Générer l'interprétation pour ce facteur
    factor_name <- paste0("Facteur ", f)
    
    # Analyser la force et le nombre de charges
    n_significant <- sum(above_threshold)
    avg_loading <- mean(abs(strong_loadings))
    
    # Évaluer la clarté du facteur
    if (n_significant == 0) {
      clarity <- "Ce facteur n'a pas de charges significatives au seuil défini."
      factor_quality <- "insuffisant"
    } else if (n_significant < 3) {
      clarity <- paste0("Ce facteur a seulement ", n_significant, " variable(s) avec des charges significatives, ce qui est généralement insuffisant pour une interprétation robuste.")
      factor_quality <- "faible"
    } else if (avg_loading < 0.4) {
      clarity <- "Ce facteur a des charges modérées, rendant son interprétation moins claire."
      factor_quality <- "modéré"
    } else if (avg_loading >= 0.6) {
      clarity <- "Ce facteur a des charges fortes et clairement définies."
      factor_quality <- "excellent"
    } else {
      clarity <- "Ce facteur a des charges de force moyenne à bonne."
      factor_quality <- "bon"
    }
    
    # Résumé des variables principales
    if (length(positive_vars) > 0) {
      pos_summary <- paste0("Variables avec charges positives : ", paste(positive_vars, collapse=", "))
    } else {
      pos_summary <- "Aucune variable avec charge positive significative"
    }
    
    if (length(negative_vars) > 0) {
      neg_summary <- paste0("Variables avec charges négatives : ", paste(negative_vars, collapse=", "))
    } else {
      neg_summary <- "Aucune variable avec charge négative significative"
    }
    
    # Suggestion d'interprétation
    suggestion <- ""
    if (n_significant >= 3) {
      suggestion <- paste0("Suggestion : Examinez les ", min(5, n_significant), " variables ayant les charges les plus élevées pour identifier le concept latent commun. ",
                           "Envisagez de nommer ce facteur en fonction du contenu théorique ou sémantique partagé par ces variables.")
    }
    
    # Regrouper toutes les informations
    interpretations[[f]] <- list(
      name = factor_name,
      quality = factor_quality,
      clarity = clarity,
      positive_vars = positive_vars,
      negative_vars = negative_vars,
      pos_summary = pos_summary,
      neg_summary = neg_summary,
      suggestion = suggestion,
      n_vars = n_significant,
      mean_loading = avg_loading
    )
  }
  
  return(interpretations)
}

# Diagnostics simples de structure. Ils ne remplacent jamais le jugement
# théorique, mais rendent visibles les charges croisées et les items faibles.
summarise_simple_structure <- function(loadings, communalities, threshold = 0.30,
                                       cross_gap = 0.15) {
  q <- ncol(loadings)
  ordre <- t(apply(abs(loadings), 1, order, decreasing = TRUE))
  if (q == 1) ordre <- matrix(1, nrow = nrow(loadings))
  idx1 <- ordre[, 1]
  prim <- loadings[cbind(seq_len(nrow(loadings)), idx1)]
  sec <- if (q > 1) loadings[cbind(seq_len(nrow(loadings)), ordre[, 2])] else rep(0, nrow(loadings))
  gap <- abs(prim) - abs(sec)
  statut <- ifelse(abs(prim) < threshold, "charge faible",
            ifelse(abs(sec) >= threshold & gap < cross_gap, "charge croisée",
            ifelse(abs(sec) >= threshold, "seconde charge notable", "structure simple")))
  data.frame(
    Item = rownames(loadings),
    Facteur_principal = paste0("F", idx1),
    Charge_principale = round(prim, 3),
    Seconde_charge = round(sec, 3),
    Écart = round(gap, 3),
    Communalité = round(communalities, 3),
    Statut = statut,
    stringsAsFactors = FALSE
  )
}

# Les solutions factorielles sont indéterminées quant au signe et à l'ordre
# des facteurs. Avant de résumer un bootstrap, on les aligne donc sur la
# solution originale; sans cela, des intervalles de charges peuvent mélanger
# artificiellement des facteurs permutés ou inversés.
all_factor_permutations <- function(x) {
  if (length(x) <= 1) return(matrix(x, nrow = 1))
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[i], all_factor_permutations(x[-i]))
  }))
}

align_bootstrap_loadings <- function(reference, candidate) {
  if (!identical(dim(reference), dim(candidate))) return(candidate)
  q <- ncol(reference)
  permutations <- all_factor_permutations(seq_len(q))
  scores <- apply(permutations, 1, function(p) {
    sum(abs(colSums(reference * candidate[, p, drop = FALSE], na.rm = TRUE)))
  })
  aligned <- candidate[, permutations[which.max(scores), ], drop = FALSE]
  for (j in seq_len(q)) {
    if (sum(reference[, j] * aligned[, j], na.rm = TRUE) < 0) {
      aligned[, j] <- -aligned[, j]
    }
  }
  colnames(aligned) <- colnames(reference)
  aligned
}

extract_factor_correlation <- function(efa_fit, nfactors) {
  fit <- if (inherits(efa_fit, "efaList")) efa_fit[[1]] else efa_fit
  psi <- tryCatch(lavaan::inspect(fit, "est")$psi, error = function(e) NULL)
  if (is.null(psi) || !is.matrix(psi) || nrow(psi) != nfactors) {
    return(diag(nfactors))
  }
  out <- tryCatch(stats::cov2cor(psi), error = function(e) diag(nfactors))
  diag(out) <- 1
  colnames(out) <- rownames(out) <- paste0("F", seq_len(nfactors))
  out
}

# ---------------------------------------------------------------------------
# Une seule fonction ajuste l'AFE et rassemble tout ce que l'application
# affiche. Elle sert à l'analyse initiale, au changement de rotation, au
# recalcul sans certains items et au bootstrap : les quatre parcours restent
# ainsi strictement cohérents (mêmes données, mêmes charges standardisées).
# ---------------------------------------------------------------------------
ajuster_afe <- function(data_complete, nfactors, rotation, method,
                        n_total = nrow(data_complete), maxit = 1000) {
  if (ncol(data_complete) <= nfactors) {
    stop("Il faut plus d'items que de facteurs demandés.")
  }
  modele <- lavaan::efa(
    data = data_complete, nfactors = nfactors, rotation = rotation,
    estimator = method,
    # Optimiseur par défaut de lavaan (nlminb) : nettement plus rapide et plus
    # stable que BFGS lorsque l'on demande un facteur de trop.
    control = list(iter.max = maxit)
  )
  fit <- if (inherits(modele, "efaList")) modele[[1]] else modele
  L <- extract_standardized_efa_loadings(modele, names(data_complete))
  if (is.null(L)) stop("Charges standardisées indisponibles : la solution n'a pas convergé.")
  h2 <- attr(L, "communalities"); attr(L, "communalities") <- NULL
  if (is.null(h2) || all(is.na(h2))) h2 <- rowSums(L^2)
  names(h2) <- rownames(L)
  phi <- extract_factor_correlation(modele, nfactors)
  ss <- colSums(L^2)
  list(
    efa = modele, loadings = L, communality = h2, uniquenesses = 1 - h2,
    method = method, rotation = rotation, nfactors = nfactors,
    # Somme des carrés des charges (patrons) par facteur. Avec une rotation
    # oblique, ces quantités se chevauchent : leur total n'est pas la variance
    # commune. La variance commune réelle est la moyenne des communalités.
    Vaccounted = data.frame(SS_loadings = ss, Proportion = ss / nrow(L),
                            Cumulative = cumsum(ss / nrow(L))),
    variance_commune = mean(h2),
    fit_indices = tryCatch(lavaan::fitMeasures(fit), error = function(e) NULL),
    converged = isTRUE(tryCatch(lavaan::lavInspect(fit, "converged"), error = function(e) FALSE)),
    factor_correlation = phi, structure_loadings = L %*% phi,
    data_complete = data_complete, variables = names(data_complete),
    n_total = n_total, n_used = nrow(data_complete)
  )
}

# Charges standardisées de POPULATION d'une simulation : la charge réglée
# est divisée par l'écart-type de l'item avant standardisation. Elle n'égale
# la valeur du curseur que si le multiplicateur de bruit vaut 1.
population_simulee <- function(loading_pattern, factor_cor, error_multiplier) {
  q <- ncol(loading_pattern)
  phi <- if (q == 1) matrix(1) else factor_cor
  commune <- diag(loading_pattern %*% phi %*% t(loading_pattern))
  var_item <- commune + (1 - commune) * error_multiplier^2
  L <- loading_pattern / sqrt(var_item)
  dimnames(L) <- list(paste0("X", seq_len(nrow(L))), paste0("F", seq_len(q)))
  list(loadings = L, communality = setNames(commune / var_item, rownames(L)), phi = phi)
}

# Oméga total d'une sous-échelle : modèle à un facteur (congénérique) ajusté
# sur les seuls items de la sous-échelle, ω = (Σλ)² / [(Σλ)² + Σψ].
omega_sous_echelle <- function(d) {
  m <- tryCatch(suppressWarnings(psych::fa(d, nfactors = 1, fm = "ml")), error = function(e) NULL)
  if (is.null(m)) return(list(omega = NA_real_, lambda = rep(NA_real_, ncol(d))))
  l <- as.numeric(m$loadings[, 1]); psi <- 1 - l^2
  list(omega = sum(l)^2 / (sum(l)^2 + sum(psi)), lambda = setNames(l, colnames(d)))
}

# Coefficient de congruence de Tucker entre deux colonnes de charges.
congruence_tucker <- function(a, b) sum(a * b) / sqrt(sum(a^2) * sum(b^2))

# Nombre de facteurs suggéré par l'analyse parallèle : on compte les valeurs
# propres observées qui dépassent le hasard AVANT le premier croisement.
suggestion_parallele <- function(observees, simulees) {
  sous <- which(observees <= simulees)
  if (length(sous) == 0) length(observees) else sous[1] - 1
}

# Nombres à la française pour les textes et les tableaux.
fr <- function(x, d = 2) {
  x <- as.numeric(x)
  out <- formatC(x, format = "f", digits = d, decimal.mark = ",", big.mark = " ")
  out[is.na(x)] <- "—"
  sub("^-", "−", out)
}
fr_sci <- function(x) {
  if (!is.finite(x) || x == 0) return(fr(x, 2))
  e <- floor(log10(abs(x)))
  if (e >= -2) return(fr(x, 4))
  paste0(fr(x / 10^e, 2), " × 10<sup>", sub("^-", "−", e), "</sup>")
}
fr_p <- function(p) ifelse(is.na(p), "—", ifelse(p < 0.001, "< 0,001", fr(p, 3)))

# Tableaux DT en français, sans dépendre d'un fichier de traduction en ligne.
options(DT.options = list(language = list(
  search = "Rechercher :", lengthMenu = "Afficher _MENU_ lignes",
  info = "Lignes _START_ à _END_ sur _TOTAL_", infoEmpty = "Aucune ligne",
  infoFiltered = "(filtrées sur _MAX_)", zeroRecords = "Aucune ligne correspondante",
  emptyTable = "Aucune donnée", paginate = list(previous = "Précédent", `next` = "Suivant"))))

# Interface utilisateur
ui <- page_sidebar(
  fillable = FALSE,
  title = "Analyse Factorielle Exploratoire",
  theme = bs_theme(bootswatch = "flatly"),
  useShinyjs(),
  tags$style(HTML("
    #loadings_table table.dataTable {
      width: auto !important;
      margin-left: 0 !important;
    }
    #loadings_table th, #loadings_table td {
      width: auto !important;
      white-space: nowrap;
      padding: 0.35rem 0.55rem;
    }
  ")),
  
  # Barre latérale avec options d'importation et paramètres d'analyse
  sidebar = sidebar(
    id = "main_sidebar",
    width = "350px",
    
    # Section d'importation des données
    accordion(
      id = "sidebar_sections",
      accordion_panel(
        "Importation des données",
        value = "import_data",
        radioButtons("import_method", "Méthode d'importation :",
                     choices = c(
                       "Charger un fichier" = "file",
                       "Données simulées" = "example",
                       "Données classiques" = "classic" 
                     ),
                     selected = "example",
                     inline = TRUE),
        
        conditionalPanel(
          condition = "input.import_method == 'file'",
          div(
            fileInput("data_file", NULL,
                      accept = c(".csv", ".sav", ".dta", ".xlsx", ".xls"),
                      buttonLabel = "Parcourir...",
                      placeholder = "Aucun fichier sélectionné"),
            
            conditionalPanel(
              condition = "input.data_file != null && input.data_file.name.toLowerCase().endsWith('.csv')",
              layout_columns(
                col_widths = c(6, 6),
                checkboxInput("header", "Noms des variables en 1ère ligne", TRUE),
                radioButtons("sep", "Séparateur (CSV)",
                             choices = c(Virgule = ",",
                                         "Point-virgule" = ";",
                                         Tabulation = "\\t"),
                             selected = ",",
                             inline = TRUE)
              ),
              layout_columns(
                col_widths = c(12),
                radioButtons("dec", "Point décimal (CSV)",
                             choices = c(Point = ".", Virgule = ","),
                             selected = ".",
                             inline = TRUE)
              )
            ),
            
            div(
              class = "border rounded p-3 mt-3",
              style = "background-color: #f8f9fa; width: 100%; min-height: 90px;",
              uiOutput("file_summary_compact")
            )
          )
        ),
        conditionalPanel(
          condition = "input.import_method == 'example'",
          h4("Laboratoire de simulation"),
          selectInput("sim_scenario", "Situation à explorer",
                      choices = c(
                        "Structure claire" = "clear",
                        "Facteurs corrélés" = "correlated",
                        "Items ambigus (charges croisées)" = "cross",
                        "Petit échantillon / structure fragile" = "fragile"
                      ), selected = "clear"),
          tags$p(class = "text-muted small",
                 "Choisissez d'abord une situation, puis modifiez les réglages si vous souhaitez expérimenter."),
          numericInput("sim_nobs", "Nombre d'observations", value = 300, min = 50, max = 1000),
          numericInput("sim_nfactors", "Nombre de facteurs latents", value = 2, min = 1, max = 5),
          numericInput("sim_nvars", "Nombre de variables par facteur", value = 4, min = 2, max = 10),
          sliderInput("sim_loading", "Force des charges principales", min = 0.5, max = 0.9, value = 0.70, step = 0.05),
          sliderInput("sim_cross_loading", "Force des charges secondaires", min = 0, max = 0.4, value = 0.05, step = 0.05),
          sliderInput("sim_cross_proportion", "Proportion d’items ambigus", min = 0, max = 1, value = 0.25, step = 0.05),
          sliderInput("sim_error", "Multiplicateur du bruit résiduel (1 = charges telles que réglées)", min = 0.1, max = 1.0, value = 1.00, step = 0.05),
          checkboxInput("sim_correlated", "Facteurs corrélés", value = TRUE),
          conditionalPanel(
            condition = "input.sim_correlated == true",
            sliderInput("sim_factor_cor", "Corrélation entre facteurs", min = 0.1, max = 0.9, value = 0.3, step = 0.05)
          ),
          numericInput("sim_seed", "Graine aléatoire (pour reproduire l'exemple)",
                       value = 1234, min = 1, max = 999999),
          layout_columns(
            col_widths = c(12),
            actionButton("generate_example", "Générer", 
                         icon = icon("random"),
                         class = "btn-success",
                         width = "100%")
          )
        ),
        conditionalPanel(
          condition = "input.import_method == 'classic'",
          div(
            h4("Sélection de jeux de données classiques"),
            selectInput("classic_dataset", "Choisir un jeu de données :",
                        choices = list(
                          "-- Sélectionner --" = "",
                          "Structures nettes" = c(
                            "Holzinger et Swineford (1939) : 9 tests cognitifs" = "hs1939",
                            "Besoins psychologiques (SDT) : 12 items" = "sdt",
                            "Propension au risque (GRiPS) : 8 items, un facteur" = "grips"
                          ),
                          "Préparer les items" = c(
                            "Estime de soi de Rosenberg : items négatifs NON recodés" = "rse",
                            "BFI : 25 items de personnalité, valeurs manquantes" = "bfi"
                          ),
                          "Structures plus exigeantes" = c(
                            "Anxiété face aux statistiques (SAQ, Field) : 23 items" = "field",
                            "Prise de risque par domaine (DOSPERT) : 30 items" = "dospert",
                            "NEO-PI-R : 30 facettes des cinq grands traits" = "neo"
                          )
                        )),
            tags$div(
              id = "dataset_description",
              style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px; margin-top: 10px;",
              uiOutput("dataset_info")
            ),
            layout_columns(
              col_widths = c(12),
              actionButton("load_classic", "Charger", 
                           icon = icon("database"),
                           class = "btn-primary",
                           width = "100%")
            )
          )
        )
      ),
      
      # Section de filtrage des données
      accordion_panel(
        "Filtrage des données",
        value = "filtering",
        conditionalPanel(
          condition = "output.data_loaded",
          selectInput("filter_var", "Variable de filtre",
                      choices = NULL),
          uiOutput("filter_controls"),
          actionButton("apply_filter", "Appliquer le filtre", 
                       icon = icon("filter"),
                       class = "btn-info",
                       width = "100%"),
          actionButton("reset_filter", "Réinitialiser", 
                       icon = icon("undo"),
                       width = "100%")
        )
      ),

      accordion_panel(
        "Préparation des items",
        value = "item_preparation",
        conditionalPanel(
          condition = "output.data_loaded",
          tags$p(class = "text-muted small",
                 "Inverser un item recode x en min + max − x. Utilisez cette étape seulement lorsque le libellé de l'item le justifie, jamais pour « améliorer » un résultat."),
          selectizeInput("reverse_items", "Items à inverser", choices = NULL,
                         multiple = TRUE, options = list(plugins = list("remove_button"))),
          layout_columns(
            col_widths = c(6, 6),
            numericInput("reverse_min", "Minimum de l'échelle", value = NA),
            numericInput("reverse_max", "Maximum de l'échelle", value = NA)
          ),
          tags$p(class = "text-muted small",
                 "Laissez vide pour utiliser le minimum et le maximum observés. Indiquez les bornes théoriques (p. ex. 1 et 5) si une catégorie extrême n'a jamais été choisie."),
          actionButton("reverse_selected_items", "Inverser les items sélectionnés",
                       icon = icon("repeat"), class = "btn-outline-primary", width = "100%"),
          actionButton("reset_preparation", "Restaurer les données importées",
                       icon = icon("undo"), class = "btn-outline-secondary mt-2", width = "100%"),
          uiOutput("preparation_summary")
        )
      ),
      
      # Nouvelle section pour les paramètres d'analyse
      accordion_panel(
        "Paramètres de l'analyse",
        value = "analysis_parameters",
        conditionalPanel(
          condition = "output.data_loaded",
          selectInput("efa_variables", "Variables à inclure", 
                      choices = NULL, multiple = TRUE),
          selectInput("efa_method", "Méthode d'extraction",
                      choices = c("Maximum de vraisemblance (défaut)" = "ML",
                                  "Maximum de vraisemblance robuste" = "MLR",
                                  "Moindres carrés non pondérés" = "ULS"),
                      selected = "ML"),

          radioButtons("measurement_level", "Nature des items",
                       choices = c("Continus ou Likert approximativement continus" = "continuous",
                                   "Likert nettement ordinal / très asymétrique" = "ordinal"),
                       selected = "continuous"),

          tags$div(class = "alert alert-info small",
                   tags$b("Repère. "),
                   "Cette version réalise l’AFE sur les corrélations de Pearson. Pour des items Likert très asymétriques ou peu catégorisés, une AFE sur corrélations polychorïques constitue une analyse distincte; ne la remplacez pas par un simple test de normalité."),
          conditionalPanel(
            condition = "input.measurement_level == 'ordinal'",
            tags$div(class = "alert alert-warning small",
                     tags$b("Limite à documenter. "),
                     "L’AFE affichée demeure fondée sur Pearson. Interprétez-la comme une approximation ou réalisez, dans un outil approprié, une analyse polychorïque/WLSMV complémentaire.")
          ),

          layout_column_wrap(
            width = 1,
            card(
              card_header("Rotation dynamique"),
              card_body(
                selectInput("dynamic_rotation", "Changer la méthode de rotation",
                            choices = c("Oblimin — facteurs possiblement corrélés" = "oblimin",
                                        "Varimax — facteurs supposés indépendants" = "varimax",
                                        "Aucune rotation (à des fins de comparaison)" = "none"),
                            selected = "oblimin"),
                tags$p(class = "text-muted small",
                       "La rotation facilite l'interprétation; elle ne décide pas à elle seule du nombre de facteurs. Après une AFE, tout changement est appliqué automatiquement.")
              )
            )
          ),
          numericInput("efa_nfactors", "Nombre de facteurs", 
                       value = 2, min = 1, max = 10),
          actionButton("run_efa", "Exécuter l'analyse", 
                       icon = icon("play"),
                       class = "btn-success",
                       width = "100%"),

        )
      )
    ),
    conditionalPanel(
      condition = "output.data_loaded",
      div(
        class = "mt-3",
        actionButton(
          "open_analysis_settings",
          "Étape suivante : définir les paramètres",
          icon = icon("sliders"),
          class = "btn-outline-secondary",
          width = "100%"
        ),
        tags$p(
          class = "text-muted small mt-2 mb-0",
          "Choisissez d'abord les variables et le nombre de facteurs dans « Paramètres de l'analyse ». Le bouton vert d'exécution s'y trouve ensuite."
        )
      )
    ),
    hr(),
    actionButton("help_button", "Aide", icon = icon("question-circle"), width = "100%")
  ),
  
  # Contenu principal avec système d'onglets
  navset_card_tab(
    title = "Analyse des données",
    id = "main_tabs",
    # Onglet d'introduction (ajouter au début des onglets principaux)
    nav_panel(
      title = "Introduction",
      icon = icon("info-circle"),
      
      layout_column_wrap(
        width = 1,
        card(
          card_header("Analyse Factorielle Exploratoire"),
          card_body(
            tags$h4("Bienvenue dans cette application d'analyse factorielle exploratoire"),
            tags$p("Cette application vous permet de réaliser des analyses factorielles exploratoires (AFE) 
               sur vos données de manière simple et interactive, tout en offrant des fonctionnalités avancées."),
            
            tags$h5("Fonctionnalités principales :"),
            tags$ul(
              tags$li(tags$b("Importation de données"), " - Chargez vos fichiers CSV, Excel, SPSS ou Stata, 
                  ou générez des données simulées pour tester l'application"),
              tags$li(tags$b("Tests préliminaires"), " - Vérifiez si vos données sont adaptées à l'AFE 
                  (KMO, Bartlett, déterminant)"),
              tags$li(tags$b("Analyse factorielle guidée"), " - Déterminez le nombre de facteurs, puis comparez une rotation Varimax et une rotation Oblimin"),
              tags$li(tags$b("Visualisation"), " - Visualisez les résultats avec des graphiques interactifs"),
              tags$li(tags$b("Interprétation raisonnée"), " - Distinguez charges principales, charges croisées et communalités")
            ),
            
            tags$h5("Comment utiliser cette application :"),
            tags$ol(
              tags$li("Commencez par importer vos données via la barre latérale à gauche"),
              tags$li("Explorez vos données dans l'onglet 'Aperçu des données'"),
              tags$li("Exécutez les tests préliminaires pour vérifier l'adéquation de vos données"),
              tags$li("Utilisez l'analyse parallèle et le graphique d'éboulis pour justifier un nombre de facteurs"),
              tags$li("Lancez l'analyse, puis interprétez les charges, les charges croisées et les communalités"),
              tags$li("Rédigez une décision : ce que les données appuient, ce qui demeure ambigu et ce qu'il faudrait vérifier ailleurs")
            )
          )
        ),
        
        card(
          card_header("Ressources et références"),
          card_body(
            tags$h5("Principes de l'analyse factorielle exploratoire"),
            tags$p("L'AFE est une technique statistique permettant d'identifier des facteurs latents 
               (non directement observables) qui expliquent les corrélations entre variables observées."),
            
            tags$h5("Étapes clés de l'AFE :"),
            tags$ol(
              tags$li("Vérification de l'adéquation des données"),
              tags$li("Détermination du nombre de facteurs à extraire"),
              tags$li("Choix d'une méthode d'extraction des facteurs"),
              tags$li("Sélection d'une méthode de rotation"),
              tags$li("Interprétation des facteurs")
            ),
            
            tags$h5("Références utiles :"),
            tags$ul(
              tags$li(tags$a(href = "https://doi.org/10.1037/1082-989X.4.3.272", target = "_blank",
                             "Fabrigar, L. R., Wegener, D. T., MacCallum, R. C. et Strahan, E. J. (1999). Evaluating the use of exploratory factor analysis in psychological research. Psychological Methods, 4(3), 272-299.")),
              tags$li(tags$a(href = "https://doi.org/10.1037/0033-2909.99.3.432", target = "_blank",
                             "Zwick, W. R. et Velicer, W. F. (1986). Comparison of five rules for determining the number of components to retain. Psychological Bulletin, 99(3), 432-442.")),
              tags$li(tags$a(href = "https://doi.org/10.1007/s12144-019-00300-2", target = "_blank",
                             "Goretzko, D., Pham, T. T. H. et Bühner, M. (2021). Exploratory factor analysis: Current use, methodological developments and recommendations for good practice. Current Psychology, 40, 3510-3521."))
            ),
            
            tags$div(
              style = "margin-top: 20px; padding: 10px; background-color: #f8f9fa; border-radius: 5px;",
              tags$p(style = "font-style: italic;", 
                     "Si vous rencontrez des difficultés dans l'utilisation de cette application, 
                 cliquez sur le bouton 'Aide' dans la barre latérale pour obtenir des informations supplémentaires.")
            )
          )
        )
      )
    ),
    # Onglet de visualisation des données
    nav_panel(
      title = "Aperçu des données",
      icon = icon("table"),
      
      # Utiliser des sous-onglets pour organiser l'aperçu des données
      navset_card_tab(
        id = "data_overview_tabs",
        
        # Sous-onglet pour l'affichage des données
        nav_panel(
          title = "Données",
          icon = icon("table"),
          
          layout_column_wrap(
            width = 1,
            card(
              card_header("Aperçu des données"),
              card_body(
                conditionalPanel(
                  condition = "output.data_loaded",
                  layout_columns(
                    col_widths = c(12),
                    conditionalPanel(
                      condition = "input.import_method == 'example'",
                      div(
                        style = "display: flex; align-items: center; gap: 15px;",
                        downloadButton("download_simulated_data", "Télécharger", 
                                       icon = icon("download"),
                                       class = "btn-info"),
                        div(
                          style = "margin-top: 5px;",
                          radioButtons("download_format", "Format :",
                                       choices = c("CSV" = "csv",
                                                   "SPSS" = "sav",
                                                   "Stata" = "dta",
                                                   "Excel" = "xlsx"),
                                       selected = "csv",
                                       inline = TRUE)
                        )
                      )
                    )
                  ),
                  DTOutput("data_preview")
                )
              )
            )
          )
        ),
        
        # Sous-onglet : la « vérité » d'une simulation, pour juger ce que l'AFE retrouve
        nav_panel(
          title = "Population simulée",
          icon = icon("bullseye"),
          card(
            card_header("Ce que la simulation a réellement généré"),
            card_body(
              uiOutput("population_info"),
              DTOutput("population_table"),
              tags$hr(),
              uiOutput("population_phi")
            )
          )
        ),

        # Sous-onglet pour les statistiques descriptives
        nav_panel(
          title = "Statistiques descriptives",
          icon = icon("chart-bar"),
          
          # Utiliser des sous-onglets pour organiser les statistiques descriptives
          navset_card_tab(
            id = "stats_tabs",
            
            # Sous-onglet pour la matrice de corrélation (optimisé pour utiliser tout l'espace)
            nav_panel(
              title = "Matrice de corrélation",
              icon = icon("table"),
              
              # Utiliser tout l'espace disponible pour la matrice
              card(
                card_header("Matrice de corrélation"),
                card_body(
                  radioButtons("corr_display_type", "Affichage :",
                               choices = c("Graphique" = "plot", "Texte" = "text"),
                               selected = "plot",
                               inline = TRUE),
                  uiOutput("correlation_output")
                )
              )
            ),
            
            # Sous-onglet pour le résumé statistique et tests de normalité (automatiques)
            nav_panel(
              title = "Résumé et distributions",
              icon = icon("info-circle"),
              
              layout_column_wrap(
                width = 1,
                card(
                  card_header("Résumé statistique"),
                  card_body(DTOutput("data_summary"))
                )
              ),

              layout_column_wrap(
                width = 1,
                card(
                  card_header("Valeurs manquantes"),
                  card_body(
                    tags$p(class = "text-muted small",
                           "L’AFE utilise actuellement les cas complets pour les variables sélectionnées. Consultez cette information avant d’interpréter l’échantillon retenu."),
                    DTOutput("missingness_table")
                  )
                )
              ),
              
              layout_column_wrap(
                width = 1,
                card(
                  card_header("Indicateurs de distribution (informatifs)"),
                  card_body(
                    tags$div(class = "alert alert-light small",
                             "Un test de Shapiro significatif ne décide pas à lui seul si une AFE est possible. Examinez aussi les graphiques et la nature des items."),
                    DTOutput("normality_results_table")
                  )
                )
              )
            ),
            
            # Sous-onglet pour la matrice de visualisations
            nav_panel(
              title = "Visualisation",
              icon = icon("chart-column"),
              
              layout_column_wrap(
                width = 1,
                card(
                  card_header("Visualisation des variables"),
                  card_body(
                    fluidRow(
                      column(6,
                             selectInput("viz_variables", "Sélectionner des variables (max 8) :", 
                                         choices = NULL, multiple = TRUE, selected = NULL)
                      ),
                      column(3,
                             radioButtons("matrix_plot_type", "Type de graphique :",
                                          choices = c("Histogrammes" = "hist", 
                                                      "Boxplots" = "box",
                                                      "QQ-Plots" = "qq"),
                                          selected = "hist",
                                          inline = TRUE)
                      ),
                      column(3,
                             checkboxInput("matrix_density", "Ajouter courbes de densité", value = TRUE)
                      )
                    ),
                    plotOutput("matrix_plot", height = "600px")
                  )
                )
              )
            )
          )
        )
      )
    ),
    
    # Onglet pour l'analyse factorielle (caché par défaut)
    nav_panel(
      title = "Analyse factorielle",
      icon = icon("chart-line"),
      
      # Utiliser des sous-onglets pour organiser le contenu
      navset_card_tab(
        id = "efa_tabs",
        
        # Sous-onglet des tests préliminaires
        nav_panel(
          title = "Tests préliminaires",
          icon = icon("check-circle"),
          tags$p(class = "text-muted small mt-2",
                 "Ces diagnostics décrivent la matrice de corrélation des items retenus (cas complets). Ils informent la décision; ils ne la prennent pas à votre place."),
          layout_columns(
            col_widths = c(4, 4, 4),
            card(card_header("Test de sphéricité de Bartlett"), card_body(uiOutput("bartlett_test"))),
            card(card_header("Déterminant de la matrice"), card_body(uiOutput("determinant_test"))),
            card(card_header("KMO global"), card_body(uiOutput("kmo_global")))
          ),
          layout_columns(
            col_widths = c(4, 8),
            card(
              card_header("KMO par item (MSA)"),
              card_body(DTOutput("kmo_items"))
            ),
            card(
              card_header("Nombre de facteurs : éboulis et analyse parallèle"),
              card_body(
                radioButtons("scree_type", NULL,
                             choices = c("Valeurs propres de l'AFE (matrice réduite)" = "fa",
                                         "Valeurs propres de l'ACP (matrice complète)" = "pc"),
                             selected = "fa", inline = TRUE),
                plotOutput("scree_plot", height = "380px"),
                uiOutput("parallel_summary")
              )
            )
          )
        ),

        # Sous-onglet des résultats principaux
        nav_panel(
          title = "Résultats principaux",
          icon = icon("list"),
          
          layout_column_wrap(
            width = 1,
            card(
              card_header("Lecture guidée de la solution"),
              card_body(
                uiOutput("efa_summary_guided")
              )
            ),
            card(
              card_header("Détails techniques"),
              card_body(
                verbatimTextOutput("efa_results")
              )
            )
          )
        ),
        
        # Sous-onglet de visualisation
        nav_panel(
          title = "Visualisation",
          icon = icon("chart-bar"),
          
          layout_columns(
            col_widths = c(3, 9),
            
            # Colonne de gauche avec tous les inputs regroupés dans une card
            card(
              card_header("Contrôles de visualisation"),
              card_body(
                selectInput("plot_factor_x", "Facteur en X:", 
                            choices = NULL,
                            selected = 1,
                            width = "150px"),
                selectInput("plot_factor_y", "Facteur en Y:", 
                            choices = NULL,
                            selected = 2,
                            width = "150px"),
                checkboxInput("auto_recalculate", "Recalculer à chaque clic", value = FALSE),
                actionButton("manual_recalculate", "Recalculer sans les items désactivés",
                             icon = icon("refresh"),
                             class = "btn-info",
                             width = "100%"),
                selectizeInput("disabled_variables", "Variables désactivées :", 
                               choices = NULL, multiple = TRUE, 
                               options = list(plugins = list("remove_button"))),
                # Ajout d'une boîte pour les informations après recalcul
                tags$hr(),
                tags$h5("Informations sur l'analyse recalculée"),
                htmlOutput("recalculated_info")
              )
            ),
            
            # Colonne de droite avec le graphique
            card(
              card_header("Visualisation des facteurs"),
              card_body(
                plotlyOutput("factor_plot", height = "500px")
              )
            )
          )
        ),
        
        # Sous-onglet de comparaison des modèles
        nav_panel(
          title = "Comparaison des modèles",
          icon = icon("balance-scale"),
          
          layout_column_wrap(
            width = 1,
            card(
              card_header("Comparer des solutions voisines"),
              card_body(
                tags$div(class = "alert alert-light small",
                         "Utilisez ce tableau pour comparer les solutions suggérées par l’analyse parallèle. Aucun indice ne décide seul : examinez aussi les items, les charges croisées et la cohérence théorique."),
                DTOutput("model_comparison_output"),
                tags$div(class = "small text-muted mt-2",
                         tags$b("Lecture. "), "Le χ² teste l'ajustement exact : avec un grand n, il rejette presque toujours. ",
                         "CFI et TLI (près de 1) comparent au modèle d'indépendance; RMSEA et SRMR (près de 0) mesurent l'écart résiduel. ",
                         "Le BIC le plus petit favorise la solution la plus parcimonieuse compatible avec les données. ",
                         "Ajouter un facteur améliore toujours un peu l'ajustement : la question est de savoir si le nouveau facteur est défini par plusieurs items nets et interprétables.")
              )
            )
          )
        ),
        # Sous-onglet pour l'analyse de fiabilité
        nav_panel(
          title = "Fiabilité des facteurs",
          icon = icon("chart-line"),
          
          layout_columns(
            col_widths = c(6, 6),  # Répartition 50-50 sur une ligne
            
            # Card 1 - Tableau des indices de fiabilité
            card(
              card_header("Indices de fiabilité"),
              card_body(
                selectInput("reliability_method", "Méthodes de fiabilité à afficher :",
                            choices = c("Alpha de Cronbach" = "alpha",
                                        "Oméga de McDonald" = "omega",
                                        "Les deux" = "both"),
                            selected = "both"),
                tags$p(class = "text-muted small",
                       "Ces indices sont exploratoires : ils dépendent des items retenus, de l’échantillon et du modèle. Une valeur élevée n’est pas à elle seule une preuve de validité."),
                DTOutput("reliability_table")
              )
            ),
            
            # Card 2 - Analyse détaillée de fiabilité
            card(
              card_header("Analyse détaillée de fiabilité"),
              card_body(
                selectInput("factor_for_analysis", "Sélectionner un facteur :", 
                            choices = NULL),
                uiOutput("detailed_reliability")
              )
            )
          )
        )
      )
    ),

    
    # Onglet pour l'interprétation des résultats (caché par défaut)
    nav_panel(
      title = "Interprétation",
      icon = icon("search"),

      # Chaque étape de l'interprétation dispose de sa page. Cela évite les
      # cartes à hauteur fixe et les défilements internes qui cachent le
      # raisonnement de l'étudiant.
      navset_card_tab(
        id = "interpretation_tabs",
        nav_panel(
          title = "1. Charges",
          icon = icon("table"),
          card(
            card_header("Charges factorielles"),
            card_body(
              layout_columns(
                col_widths = c(4, 8),
                sliderInput("loadings_threshold", "Seuil de lecture", 
                            min = 0, max = 1, value = 0.3, step = 0.05),
                tags$div(
                  class = "alert alert-light mb-0",
                  "Le signe d'une charge est arbitraire : interprétez son ampleur et le contenu de l'item. Une seconde charge élevée signale une ambiguïté possible."
                )
              ),
              uiOutput("loadings_warning"),
              DTOutput("loadings_table")
            )
          )
        ),
        nav_panel(
          title = "2. Structure des items",
          icon = icon("list-check"),
          card(
            card_header("À vérifier avant de nommer les facteurs"),
            card_body(
              tags$p("Commencez ici : un item à charge faible ou croisée ne distingue pas nettement les facteurs. C'est une invitation à relire l'item et la théorie, non une règle de suppression automatique."),
              DTOutput("simple_structure_table")
            )
          )
        ),
        nav_panel(
          title = "3. Facteurs corrélés",
          icon = icon("link"),
          card(
            card_header("Patrons, structure et corrélations"),
            card_body(
              tags$p("Avec Oblimin, les patrons décrivent les contributions propres; la structure décrit les corrélations item–facteur. Les deux tableaux coïncident lorsque les facteurs sont indépendants."),
              DTOutput("factor_correlations_table"),
              tags$hr(),
              DTOutput("structure_loadings_table")
            )
          )
        ),
        nav_panel(
          title = "4. Qualité de la solution",
          icon = icon("chart-column"),
          layout_columns(
            col_widths = c(6, 6),
            card(
              card_header("Variance expliquée"),
              card_body(plotlyOutput("variance_plot"), uiOutput("variance_note"))
            ),
            card(
              card_header("Communalités"),
              card_body(DTOutput("communalities_table"))
            )
          )
        ),
        nav_panel(
          title = "5. Nommer les facteurs",
          icon = icon("lightbulb"),
          card(
            card_header("Pistes d'interprétation — à confronter à la théorie"),
            card_body(
              p("Ces pistes sont générées à partir des charges. Elles ne peuvent pas nommer un construit à votre place : utilisez le libellé et le contenu des items."),
              uiOutput("auto_interpretation")
            )
          )
        )
      )
    ),
    
# Analyse bootstrap
nav_panel(
  title = "Bootstrap",
  icon = icon("random"),
  
  layout_columns(
    col_widths = c(12, 12),
    fill = FALSE,
    
    # Tableau des intervalles de confiance pour les charges
    card(
      fill = FALSE,
      card_header("Intervalles de confiance des charges factorielles"),
      card_body(
        layout_columns(
          col_widths = c(6, 6),
          numericInput("bootstrap_samples", "Nombre d'échantillons", 
                       value = 100, min = 50, max = 1000, step = 50),
          numericInput("bootstrap_ci", "Niveau de confiance (%)", 
                       value = 95, min = 80, max = 99, step = 1)
        ),
        actionButton("run_bootstrap", "Exécuter l'analyse bootstrap", 
                     icon = icon("play"),
                     class = "btn-success",
                     width = "100%"),
        tags$div(
          class = "mt-3",
          style = "font-size: 0.9em; color: #6c757d;",
          "Le bootstrap estime la stabilité interne des charges par rééchantillonnage. Les solutions sont alignées sur la solution originale avant de résumer les charges; il ne s’agit pas d’une validation externe."
        ),
        tags$hr(),
        selectInput("bootstrap_display_factor", "Afficher le facteur :", choices = NULL),
        DTOutput("bootstrap_loadings_table")
      )
    ),
    
    # Visualisation de la stabilité
    card(
      fill = FALSE,
      card_header("Stabilité des charges"),
      card_body(
        selectInput("bootstrap_stability_var", "Sélectionner une variable :", choices = NULL),
        plotlyOutput("bootstrap_stability_plot")
      )
    )
  )
),

    # Validation CFA : disponible seulement après une AFE sur des données
    # simulées. La CFA porte obligatoirement sur un nouvel échantillon.
    nav_panel(
      title = "Validation CFA",
      icon = icon("check-circle"),
      layout_sidebar(
        fillable = FALSE,
        fill = FALSE,
        sidebar = sidebar(
          title = "Contrôles CFA",
          width = "300px",
          tags$h6("1. Nouvel échantillon"),
          numericInput("cfa_validation_seed", "Graine du nouvel échantillon",
                       value = 1235, min = 1, max = 999999),
          actionButton("generate_validation_sample", "Générer le nouvel échantillon",
                       icon = icon("shuffle"), class = "btn-primary", width = "100%"),
          tags$hr(),
          tags$h6("2. Modèle candidat"),
          sliderInput("cfa_threshold", "Seuil de charge retenu", 
                      min = 0, max = 1, value = 0.3, step = 0.05),
          radioButtons("cfa_crossloading_rule", "Items à charge croisée",
                       choices = c("Conserver seulement la charge principale" = "primary",
                                   "Conserver toutes les charges au-dessus du seuil" = "all"),
                       selected = "primary"),
          actionButton("generate_cfa", "Générer la syntaxe", 
                       icon = icon("cogs"), class = "btn-outline-primary", width = "100%"),
          selectInput("cfa_method", "Méthode d'estimation",
                      choices = c("Maximum de vraisemblance" = "ML", 
                                  "Maximum de vraisemblance robuste" = "MLR",
                                  "Moindres carrés non pondérés" = "ULS"),
                      selected = "MLR"),
          conditionalPanel(
            condition = "output.validation_sample_available",
            actionButton("run_cfa", "Évaluer le modèle", 
                         icon = icon("play"), class = "btn-success", width = "100%")
          ),
          conditionalPanel(
            condition = "!output.validation_sample_available",
            tags$p(class = "text-muted small",
                   "Générez d'abord le nouvel échantillon.")
          )
        ),
        navset_card_tab(
          id = "cfa_main_tabs",
          nav_panel(
            title = "Consignes",
            icon = icon("map-signs"),
            card(
              card_header("De l'AFE à la CFA"),
              card_body(
                tags$div(
                  class = "alert alert-info",
                  tags$b("La logique. "),
                  "L'AFE a servi à proposer une structure. La CFA évalue cette structure sur de nouvelles données. L'échantillon de validation est tiré de la même population simulée : les paramètres générateurs sont identiques, mais les personnes et l'erreur aléatoire sont nouvelles."
                ),
                tags$ol(
                  tags$li("Générez le nouvel échantillon dans le panneau de gauche. Une graine différente est obligatoire."),
                  tags$li("Générez la syntaxe proposée par l'AFE, puis relisez-la à la lumière de la théorie."),
                  tags$li("Évaluez le modèle et interprétez l'ajustement avant de consulter les indices de modification.")
                ),
                uiOutput("validation_sample_info")
              )
            )
          ),
          nav_panel(
            title = "Modèle lavaan",
            icon = icon("code"),
            card(
              card_header("Syntaxe à confirmer"),
              card_body(
                tags$p("Cette syntaxe est une proposition issue de l'AFE. Une modification est permise seulement si elle peut être justifiée par le contenu des items ou par la théorie, et non simplement parce qu'elle améliore un indice d'ajustement."),
                uiOutput("cfa_crossloading_warning"),
                textAreaInput("cfa_syntax", "Syntaxe lavaan", rows = 7, width = "100%",
                              placeholder = "Générez la syntaxe à partir de l'AFE après avoir créé l'échantillon de validation.")
              )
            )
          ),
          nav_panel(
            title = "Résultats et interprétation",
            icon = icon("chart-line"),
            card(
              card_header("Lecture guidée de l'ajustement"),
              card_body(uiOutput("cfa_interpretation"))
            ),
            navset_card_tab(
              nav_panel(title = "Charges standardisées", DTOutput("cfa_loadings_table"), uiOutput("cfa_phi")),
              nav_panel(title = "Sortie lavaan", verbatimTextOutput("cfa_summary")),
              nav_panel(title = "Ajustement détaillé", verbatimTextOutput("cfa_fit_indices")),
              nav_panel(title = "Diagramme", plotOutput("cfa_diagram", height = "520px")),
              nav_panel(title = "Indices de modification", DTOutput("modification_indices"))
            )
          )
        )
      )
    )
  )
)

# Serveur
server <- function(input, output, session) {
  
  # Définir l'onglet Introduction comme onglet actif au démarrage
  updateTabsetPanel(session, "main_tabs", selected = "Introduction")
  
  # Réactives pour stocker les données
  original_data <- reactiveVal(NULL)
  raw_data <- reactiveVal(NULL)
  filtered_data <- reactiveVal(NULL)
  preparation_log <- reactiveVal(data.frame(
    Item = character(), Ancien_min = numeric(), Ancien_max = numeric(),
    Transformation = character(), stringsAsFactors = FALSE
  ))
  efa_result <- reactiveVal(NULL)
  model_comparison <- reactiveVal(NULL)
  simulation_spec <- reactiveVal(NULL)
  validation_data <- reactiveVal(NULL)
  cfa_result <- reactiveVal(NULL)
  
  # Fonction pour nettoyer les objets réactifs
  cleanup_reactive_objects <- function() {
    # Réinitialiser les objets réactifs principaux
    efa_result(NULL)
    model_comparison(NULL)
    bootstrap_results(NULL)
    cfa_result(NULL)
    validation_data(NULL)
    normality_results_data(NULL)
    hideTab(inputId = "main_tabs", target = "Interprétation")
    hideTab(inputId = "main_tabs", target = "Bootstrap")
    hideTab(inputId = "main_tabs", target = "Validation CFA")
    
    # Libérer la mémoire explicitement
    gc()
    
    # Message de débogage en console (pas visible par l'utilisateur)
    message("Nettoyage de la mémoire effectué")
  }
  
  # Toute modification des données (inversion, filtre) rend la solution
  # affichée obsolète : on la retire plutôt que de laisser des résultats
  # qui ne correspondent plus aux données.
  invalider_resultats <- function() {
    efa_result(NULL); model_comparison(NULL); bootstrap_results(NULL); cfa_result(NULL)
    hideTab(inputId = "main_tabs", target = "Interprétation")
    hideTab(inputId = "main_tabs", target = "Bootstrap")
    hideTab(inputId = "main_tabs", target = "Validation CFA")
  }

  # Indicateur de chargement des données
  output$data_loaded <- reactive({
    return(!is.null(raw_data()))
  })
  outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)

  output$validation_sample_available <- reactive({
    !is.null(validation_data())
  })
  outputOptions(output, "validation_sample_available", suspendWhenHidden = FALSE)
  
  # Cacher les onglets d'analyse au démarrage
  observe({
    hideTab(inputId = "main_tabs", target = "Analyse factorielle")
    hideTab(inputId = "main_tabs", target = "Interprétation")
    # Ces deux modules demandent des précautions qui dépassent le parcours
    # d'initiation; ils demeurent dans le code, mais ne sont pas exposés ici.
    hideTab(inputId = "main_tabs", target = "Bootstrap")
    hideTab(inputId = "main_tabs", target = "Validation CFA")
  })

  # Compatibilité avec les versions récentes et antérieures de bslib.
  set_main_sidebar <- function(open) {
    if ("toggle_sidebar" %in% getNamespaceExports("bslib")) {
      bslib::toggle_sidebar(id = "main_sidebar", open = open, session = session)
    } else {
      bslib::sidebar_toggle(id = "main_sidebar", open = open, session = session)
    }
  }

  # La validation CFA possède son propre panneau de commandes. Lorsque cet
  # onglet est actif, le panneau général de l'AFE se replie pour libérer toute
  # la largeur utile au modèle et à son interprétation.
  observeEvent(input$main_tabs, {
    set_main_sidebar(!identical(input$main_tabs, "Validation CFA"))
  }, ignoreInit = FALSE)

  # Raccourci pédagogique placé sous l'accordéon : le bouton ne lance rien,
  # il rend simplement l'étape de préparation impossible à manquer.
  observeEvent(input$open_analysis_settings, {
    accordion_panel_open("sidebar_sections", "analysis_parameters", session = session)
    showNotification(
      "Paramètres ouverts : sélectionnez au moins deux variables, puis cliquez sur « Exécuter l'analyse ».",
      type = "message", duration = 6
    )
  })
  
  # Chargement automatique lorsqu'un fichier est sélectionné
  observeEvent(input$data_file, {
    req(input$data_file)
    
    # Nettoyer les objets existants avant de charger de nouvelles données
    cleanup_reactive_objects()
    
    # Afficher une notification pendant le chargement
    id <- showNotification("Chargement des données en cours...", type = "message", duration = NULL)
    
    # Obtenir l'extension et le chemin du fichier
    file_ext <- tools::file_ext(input$data_file$name)
    file_path <- input$data_file$datapath
    
    # Importer les données selon le format
    tryCatch({
      data <- NULL
      
      if (file_ext == "csv") {
        # Pour CSV, assurer que les options sont bien initialisées
        Sys.sleep(0.3) # court délai pour assurer que les widgets sont initialisés
        
        # S'assurer que les contrôles sont disponibles
        req(input$header, input$sep, input$dec)
        
        # Afficher un message de débogage
        message("Chargement CSV avec: sep=", input$sep, ", dec=", input$dec, ", header=", input$header)
        
        data <- read.csv(file_path, 
                         header = input$header, 
                         sep = input$sep, 
                         dec = input$dec,
                         stringsAsFactors = FALSE)
      } else if (file_ext %in% c("xls", "xlsx")) {
        data <- as.data.frame(readxl::read_excel(file_path))
      } else if (file_ext == "sav") {
        data <- as.data.frame(haven::read_spss(file_path))
      } else if (file_ext == "dta") {
        data <- as.data.frame(haven::read_stata(file_path))
      } else {
        stop("Format de fichier non pris en charge")
      }
      
      # Correction des noms de colonnes vides
      if (any(names(data) == "")) {
        empty_names <- which(names(data) == "")
        new_names <- names(data)
        new_names[empty_names] <- paste0("Colonne_", empty_names)
        names(data) <- new_names
      }
      
      # Conversion des variables labellées
      if (file_ext %in% c("sav", "dta")) {
        for (col in names(data)) {
          if (inherits(data[[col]], "haven_labelled")) {
            if (is.numeric(haven::zap_labels(data[[col]]))) {
              data[[col]] <- as.numeric(haven::zap_labels(data[[col]]))
            } else {
              data[[col]] <- as.factor(haven::as_factor(data[[col]]))
            }
          }
        }
      }
      
      # Stocker les données
      simulation_spec(NULL)
      validation_data(NULL)
      cfa_result(NULL)
      hideTab(inputId = "main_tabs", target = "Validation CFA")
      original_data(data)
      raw_data(data)
      filtered_data(data)
      preparation_log(preparation_log()[0, ])
      
      # Mettre à jour les sélecteurs
      updateSelectInput(session, "filter_var", 
                        choices = c("Aucun" = "", names(data)))
      
      numeric_cols <- names(data)[sapply(data, is.numeric)]
      updateSelectInput(session, "efa_variables", 
                        choices = numeric_cols)
      updateSelectizeInput(session, "reverse_items", choices = numeric_cols,
                           selected = character(0))
      
      # Afficher l'onglet d'analyse
      showTab(inputId = "main_tabs", target = "Analyse factorielle")
      
      removeNotification(id)
      showNotification("Données importées avec succès", type = "message")
    }, error = function(e) {
      removeNotification(id)
      showNotification(paste("Lecture impossible :", conditionMessage(e),
                             "Vérifiez le séparateur et le point décimal (CSV)."), type = "error", duration = 12)
    })
  })
  
  # Créer une réactive pour stocker les résultats des tests de normalité
  normality_results_data <- reactiveVal(NULL)
  
  # Changer la rotation réajuste la MÊME solution (mêmes items, mêmes cas,
  # même nombre de facteurs) dans une nouvelle orientation.
  observeEvent(input$dynamic_rotation, ignoreInit = TRUE, {
    req(efa_result())
    ancien <- isolate(efa_result())
    bootstrap_results(NULL)
    cfa_result(NULL)
    withProgress(message = "Application de la nouvelle rotation...", {
      tryCatch({
        efa_result(ajuster_afe(ancien$data_complete, ancien$nfactors,
                               input$dynamic_rotation, ancien$method,
                               n_total = ancien$n_total))
        showNotification(paste("Rotation", input$dynamic_rotation, "appliquée."), type = "message")
      }, error = function(e) {
        showNotification(paste("Erreur lors de l'application de la rotation :", e$message),
                         type = "error", duration = 10)
      })
    })
  })

  # Visualisation des facteurs
  output$factor_plot <- renderPlotly({
    req(efa_result())
    efa <- efa_result()
    
    # Définir les indices des facteurs à afficher
    factor_idx_x <- 1
    factor_idx_y <- 2
    
    # Si plus de 2 facteurs et des sélecteurs sont utilisés
    if (efa$nfactors > 2 && !is.null(input$plot_factor_x) && !is.null(input$plot_factor_y)) {
      factor_idx_x <- as.numeric(input$plot_factor_x)
      factor_idx_y <- as.numeric(input$plot_factor_y)
    }
    
    # Extraire les charges factorielles pour les facteurs sélectionnés
    loadings <- efa$loadings
    
    if (ncol(loadings) < 2) {
      # Si un seul facteur est extrait
      var_names <- rownames(loadings)
      factor1 <- loadings[, 1]
      factor2 <- rep(0, nrow(loadings))
    } else {
      var_names <- rownames(loadings)
      factor1 <- loadings[, factor_idx_x]
      factor2 <- loadings[, factor_idx_y]
    }
    
    # Vérifier s'il y a des variables désactivées
    disabled_vars <- input$disabled_variables
    if (is.null(disabled_vars)) {
      disabled_vars <- c()
    }
    
    # Créer un data frame pour le graphique
    plot_data <- data.frame(
      Variable = var_names,
      Factor1 = factor1,
      Factor2 = factor2,
      Status = ifelse(var_names %in% disabled_vars, "Désactivée", "Active")
    )
    
    # Obtenir les pourcentages de variance pour l'étiquetage
    variance_percent <- round(efa$Vaccounted$Proportion * 100, 1)
    
    # Créer le graphique
    p <- plot_ly(plot_data, x = ~Factor1, y = ~Factor2, source = "biplot",
                 text = ~Variable,
                 type = "scatter", 
                 mode = "markers+text", 
                 textposition = "top center",
                 color = ~Status,
                 colors = c("Active" = "blue", "Désactivée" = "gray"),
                 marker = list(size = 10),
                 customdata = ~Variable) %>%
      layout(
        title = "Charges standardisées (patrons) des items<br><span style='font-size:12px'>Un item près d'un axe est spécifique à ce facteur; un item entre deux axes est ambigu. Cliquez sur un point pour le désactiver.</span>",
        separators = ", ",
        xaxis = list(title = paste0("Charge sur F", factor_idx_x), range = c(-1.05, 1.05), zeroline = FALSE),
        yaxis = list(title = paste0("Charge sur F", factor_idx_y), range = c(-1.05, 1.05), zeroline = FALSE,
                     scaleanchor = "x"),
        shapes = list(
          list(
            type = "line", line = list(color = "grey", dash = "dash"), 
            x0 = 0, x1 = 0, y0 = -1, y1 = 1
          ),
          list(
            type = "line", line = list(color = "grey", dash = "dash"), 
            y0 = 0, y1 = 0, x0 = -1, x1 = 1
          )
        ),
        legend = list(title = list(text = "Statut de la variable"))
      ) %>%
      event_register("plotly_click")
    
    return(p)
  })
      
      # Tests préliminaires : un seul calcul, réutilisé par toutes les cartes.
      prelim <- reactive({
        req(filtered_data(), input$efa_variables, length(input$efa_variables) >= 2)
        d <- na.omit(filtered_data()[, input$efa_variables, drop = FALSE])
        validate(need(nrow(d) > ncol(d), "Trop peu de cas complets pour ces items."))
        R <- cor(d)
        list(d = d, R = R, n = nrow(d),
             bartlett = psych::cortest.bartlett(R, n = nrow(d)),
             det = det(R),
             kmo = tryCatch(psych::KMO(R), error = function(e) NULL))
      })

      encadre <- function(classe, ...) tags$div(class = paste("alert", classe, "small mb-0 mt-2"), ...)

      output$bartlett_test <- renderUI({
        b <- prelim()$bartlett
        tagList(
          tags$p(HTML(paste0("χ²(", b$df, ") = ", fr(b$chisq, 1), "<br>p ", 
                             ifelse(b$p.value < .001, "< 0,001", paste0("= ", fr(b$p.value, 3)))))),
          if (b$p.value < .05)
            encadre("alert-light", "Les items ne sont pas indépendants : une AFE est envisageable. Avec un grand n, ce résultat est presque attendu; il ne dit pas si la structure est utile.")
          else
            encadre("alert-warning", "On ne peut pas rejeter une matrice identité : les corrélations sont trop faibles pour qu'une structure commune soit plausible.")
        )
      })

      output$determinant_test <- renderUI({
        dt <- prelim()$det
        tagList(
          tags$p(HTML(paste0("|R| = ", fr_sci(dt)))),
          if (dt > 1e-5)
            encadre("alert-light", "Aucune redondance extrême : aucun item n'est une combinaison quasi parfaite des autres.")
          else
            encadre("alert-warning", "Déterminant très proche de 0 : certains items sont presque redondants. L'estimation peut devenir instable.")
        )
      })

      output$kmo_global <- renderUI({
        k <- prelim()$kmo
        validate(need(!is.null(k), "KMO non calculable (matrice singulière)."))
        m <- k$MSA
        lib <- cut(m, c(-Inf, .5, .6, .7, .8, .9, Inf),
                   c("inacceptable", "médiocre", "moyen", "bon", "très bon", "excellent"))
        tagList(
          tags$p(HTML(paste0("KMO = ", fr(m, 3), " (", lib, ")"))),
          encadre(if (m < .6) "alert-warning" else "alert-light",
                  "Le KMO compare les corrélations simples aux corrélations partielles. Près de 1 : les corrélations s'expliquent bien par des facteurs communs. Les libellés de Kaiser (1974) sont des repères, pas des seuils.")
        )
      })

      output$kmo_items <- renderDT({
        k <- prelim()$kmo
        req(k)
        tab <- data.frame(Item = names(k$MSAi), MSA = round(k$MSAi, 3))
        tab <- tab[order(tab$MSA), ]
        tab$Repère <- as.character(cut(tab$MSA, c(-Inf, .5, .6, .7, .8, Inf),
                                       c("à examiner", "médiocre", "moyen", "bon", "excellent")))
        datatable(tab, rownames = FALSE, options = list(dom = if (nrow(tab) > 10) "tp" else "t", pageLength = 10)) %>%
          formatRound("MSA", 3, dec.mark = ",") %>%
          formatStyle("MSA", color = styleInterval(.5, c("#c0392b", "black")),
                      fontWeight = styleInterval(.5, c("bold", "normal")))
      })

      # Analyse parallèle (Horn, 1965). La graine fixe rend la courbe simulée
      # identique d'un affichage à l'autre.
      parallele <- reactive({
        p <- prelim()
        fm <- switch(input$efa_method, ML = "ml", ULS = "uls", "minres")
        set.seed(2026)
        pa <- suppressWarnings(suppressMessages(
          psych::fa.parallel(p$d, fa = "both", fm = fm, n.iter = 100, plot = FALSE)))
        list(fa_obs = pa$fa.values, fa_sim = pa$fa.sim, pc_obs = pa$pc.values, pc_sim = pa$pc.sim,
             sugg_fa = suggestion_parallele(pa$fa.values, pa$fa.sim),
             sugg_pc = suggestion_parallele(pa$pc.values, pa$pc.sim),
             kaiser = sum(pa$pc.values > 1))
      })

      output$scree_plot <- renderPlot({
        pa <- parallele()
        type <- if (is.null(input$scree_type)) "fa" else input$scree_type
        obs <- if (type == "fa") pa$fa_obs else pa$pc_obs
        sim <- if (type == "fa") pa$fa_sim else pa$pc_sim
        sugg <- if (type == "fa") pa$sugg_fa else pa$sugg_pc
        k <- seq_along(obs)
        par(mar = c(4.2, 4.5, 2.5, 1), las = 1)
        yl <- range(c(obs, sim, if (type == "pc") 1, 0)); yl[2] <- yl[2] * 1.25
        plot(k, obs, type = "b", pch = 19, col = "#2c7fb8", lwd = 2, ylim = yl, xaxt = "n", yaxt = "n",
             xlab = "Rang de la valeur propre", ylab = "Valeur propre",
             main = if (type == "fa") "Matrice réduite (communalités sur la diagonale)"
                    else "Matrice complète (1 sur la diagonale)", font.main = 1, cex.main = 1)
        axis(1, at = k)
        at <- pretty(yl); axis(2, at = at, labels = sub("\\.", ",", format(at)))
        abline(h = 0, col = "grey80")
        lines(k, sim, type = "b", pch = 17, col = "#d95f0e", lty = 2, lwd = 2)
        if (type == "pc") abline(h = 1, col = "grey40", lty = 3)
        if (sugg > 0) {
          rect(0.5, par("usr")[3], sugg + 0.5, par("usr")[4], col = adjustcolor("#2c7fb8", 0.08), border = NA)
        }
        legend("topright", bty = "n", cex = 0.85, inset = c(0.01, 0),
               legend = c("Données observées", "Données aléatoires (moyenne de 100 simulations)",
                          if (type == "pc") "Règle de Kaiser (valeur propre = 1)"),
               col = c("#2c7fb8", "#d95f0e", if (type == "pc") "grey40"),
               lty = c(1, 2, if (type == "pc") 3), pch = c(19, 17, NA), lwd = 2)
      })

      output$parallel_summary <- renderUI({
        pa <- parallele()
        tags$div(class = "alert alert-info small mb-0",
          tags$b("Repères. "),
          "Analyse parallèle (AFE) : ", tags$b(pa$sugg_fa), " facteur(s); ",
          "analyse parallèle (ACP) : ", pa$sugg_pc, " composante(s); ",
          paste0("règle de Kaiser (valeurs propres de l'ACP > 1) : ", pa$kaiser, ". "),
          "La zone ombrée couvre les valeurs propres qui dépassent le hasard. Ces repères s'examinent avec le contenu des items et les solutions voisines (onglet « Comparaison des modèles »).")
      })

      # EFA Results
      output$efa_results <- renderPrint({
        req(efa_result())
        efa <- efa_result()
        cat("Extraction :", efa$method, " | Rotation :", efa$rotation,
            " | Facteurs :", efa$nfactors, " | Cas analysés :", efa$n_used, "sur", efa$n_total, "\n")
        if (!efa$converged) cat("ATTENTION : l'estimation n'a pas convergé.\n")
        cat("\nSortie lavaan (charges non masquées sous 0,20 remplacées par des points) :\n\n")
        print(summary(efa$efa, nd = 3, cutoff = 0.2, dot.cutoff = 0.05))
      })

      # Loadings Table
      output$loadings_warning <- renderUI({
        req(efa_result())
        max_loading <- max(abs(efa_result()$loadings), na.rm = TRUE)
        if (is.finite(max_loading) && max_loading > 1.001) {
          tags$div(
            class = "alert alert-warning py-2",
            tags$b("Solution à examiner. "),
            "Au moins une charge standardisée dépasse 1. Cela indique généralement une solution impropre (par exemple une unicité négative), et non une très forte relation à célébrer. Revoyez le nombre de facteurs, les items redondants et les diagnostics."
          )
        }
      })

      output$efa_summary_guided <- renderUI({
        req(efa_result())
        efa <- efa_result()
        st <- summarise_simple_structure(efa$loadings, efa$communality,
                                         threshold = input$loadings_threshold, cross_gap = .15)
        corr <- efa$factor_correlation
        max_r <- if (ncol(corr) > 1) max(abs(corr[upper.tri(corr)])) else 0
        f <- efa$fit_indices
        g <- function(n) if (!is.null(f) && n %in% names(f)) unname(f[n]) else NA
        tags$div(
          class = "alert alert-info mb-0",
          tags$p(tags$b("Solution estimée. "), efa$nfactors, " facteur(s), extraction ",
                 efa$method, ", rotation ", efa$rotation, "."),
          tags$ul(
            tags$li(paste0("Échantillon : ", efa$n_used, " cas complets sur ", efa$n_total,
                           " observations; ", length(efa$variables), " items.")),
            tags$li(paste0("Variance commune reproduite (moyenne des communalités) : ",
                           fr(100 * efa$variance_commune, 1), " %.")),
            tags$li(paste0(sum(st$Statut == "structure simple"), " item(s) à structure simple; ",
                           sum(st$Statut == "seconde charge notable"), " avec une seconde charge notable; ",
                           sum(st$Statut == "charge croisée"), " charge(s) croisée(s); ",
                           sum(st$Statut == "charge faible"), " charge(s) faible(s) (seuil de lecture ",
                           fr(input$loadings_threshold, 2), ").")),
            if (efa$nfactors > 1) tags$li(paste0("Plus forte corrélation absolue entre facteurs : ",
                                                 fr(max_r, 2),
                                                 if (efa$rotation == "varimax") " (nulle par construction avec Varimax)." else ".")),
            if (!is.na(g("rmsea"))) tags$li(paste0("Ajustement : RMSEA = ", fr(g("rmsea"), 3),
                                                    ", SRMR = ", fr(g("srmr"), 3),
                                                    ", CFI = ", fr(g("cfi"), 3), ".")),
            if (!efa$converged) tags$li(tags$b("L'estimation n'a pas convergé : n'interprétez pas cette solution."))
          ),
          tags$p(class = "mb-0", "La décision demeure à confronter à l'analyse parallèle, au contenu des items et aux solutions voisines.")
        )
      })

      output$loadings_table <- renderDT({
        req(efa_result())
        loadings <- efa_result()$loadings
        threshold <- input$loadings_threshold
        
        # Appliquer le seuil
        loadings_thresholded <- loadings
        loadings_thresholded[abs(loadings_thresholded) < threshold] <- NA
        
        # Formater pour l'affichage
        loadings_formatted <- as.data.frame(round(loadings_thresholded, 3))
        needs_horizontal_scroll <- ncol(loadings_formatted) > 5
        
        # Créer un tableau avec coloration selon la valeur
        datatable(
          loadings_formatted,
          options = list(dom = "t", pageLength = 200, ordering = FALSE, scrollX = needs_horizontal_scroll, autoWidth = TRUE),
          rownames = TRUE,
          class = "compact nowrap",
          width = "auto"
        ) %>%
          formatRound(columns = colnames(loadings_formatted), digits = 3, dec.mark = ",") %>%
          formatStyle(
            columns = colnames(loadings_formatted),
            backgroundColor = styleInterval(c(-threshold, threshold), c("#fde2e1", "white", "#e3f1da")),
            color = styleInterval(c(-threshold, threshold), c("#b03a2e", "black", "#1e6b2f")),
            fontWeight = styleInterval(c(-threshold, threshold), c("bold", "normal", "bold"))
          )
      })

      # Tableau explicitant les ambiguïtés de la solution. Le seuil de seconde
      # charge et l'écart minimal sont volontairement visibles dans le texte,
      # afin d'éviter qu'un code couleur soit interprété comme un verdict.
      output$simple_structure_table <- renderDT({
        req(efa_result())
        efa <- efa_result()
        structure <- summarise_simple_structure(
          efa$loadings, efa$communality,
          threshold = input$loadings_threshold,
          cross_gap = .15
        )

        names(structure) <- c("Item", "Facteur principal", "Charge principale", "Seconde charge",
                              "Écart", "Communalité", "Statut")
        datatable(
          structure,
          options = list(pageLength = 15, scrollX = TRUE, dom = if (nrow(structure) > 15) "tip" else "t"),
          rownames = FALSE,
          caption = paste0("Statuts (seuil ", fr(input$loadings_threshold, 2), ", écart minimal 0,15) : ",
                           "charge faible = charge principale sous le seuil; charge croisée = seconde charge au-dessus du seuil ",
                           "et écart inférieur à 0,15; seconde charge notable = au-dessus du seuil, mais nettement plus faible que la principale. ",
                           "Ce sont des invitations à relire l'item, pas des verdicts.")
        ) %>%
          formatRound(c("Charge principale", "Seconde charge", "Écart", "Communalité"), 3, dec.mark = ",") %>%
          formatStyle(
            "Statut",
            target = "row",
            backgroundColor = styleEqual(
              c("structure simple", "seconde charge notable", "charge croisée", "charge faible"),
              c("#e8f5e9", "#fdf6e3", "#fff3cd", "#f8d7da")
            )
          )
      })

      output$factor_correlations_table <- renderDT({
        req(efa_result())
        corr <- efa_result()$factor_correlation
        if (is.null(corr) || ncol(corr) < 2) {
          return(datatable(data.frame(Message = "Une corrélation entre facteurs exige au moins deux facteurs."),
                           options = list(dom = "t"), rownames = FALSE))
        }
        corr[upper.tri(corr)] <- NA
        datatable(as.data.frame(round(corr, 3)), options = list(dom = "t", ordering = FALSE), rownames = TRUE,
                  caption = if (efa_result()$rotation == "varimax")
                    "Varimax impose des facteurs non corrélés : ces corrélations sont nulles par construction, pas parce que les données le montrent."
                  else "Corrélations estimées entre les facteurs (Φ).") %>%
          formatRound(colnames(corr), 3, dec.mark = ",")
      })

      output$structure_loadings_table <- renderDT({
        req(efa_result())
        efa <- efa_result()
        pattern <- efa$loadings
        structure <- efa$structure_loadings
        if (is.null(structure)) structure <- pattern
        out <- cbind(
          setNames(as.data.frame(round(pattern, 3)), paste("Patron", colnames(pattern))),
          setNames(as.data.frame(round(structure, 3)), paste("Structure", colnames(structure)))
        )
        datatable(out, options = list(pageLength = 15, scrollX = TRUE, dom = if (nrow(out) > 15) "tip" else "t"), rownames = TRUE,
                  caption = "Patrons (contributions propres, P) à gauche; structure (corrélations item–facteur, S = PΦ) à droite.") %>%
          formatRound(names(out), 3, dec.mark = ",")
      })
      
      # Variance Plot
      output$variance_plot <- renderPlotly({
        req(efa_result())
        efa <- efa_result()
        part <- 100 * efa$Vaccounted$Proportion
        totale <- 100 * efa$variance_commune
        oblique <- efa$rotation == "oblimin" && efa$nfactors > 1
        p <- plot_ly(x = paste0("F", seq_len(efa$nfactors)), y = part, type = "bar",
                     marker = list(color = "steelblue"), name = "Part propre au facteur",
                     text = paste0(sub("\\.", ",", sprintf("%.1f", part)), " %"),
                     textposition = "outside", hoverinfo = "text") %>%
          layout(
            title = list(text = "Variance des items associée à chaque facteur", font = list(size = 14)),
            xaxis = list(title = ""),
            yaxis = list(title = "% de la variance totale des items", ticksuffix = " %",
                         range = c(0, max(c(part, totale)) * 1.25)),
            separators = ", ",
            shapes = list(list(type = "line", xref = "paper", x0 = 0, x1 = 1, y0 = totale, y1 = totale,
                               line = list(dash = "dash", color = "#c0392b"))),
            annotations = list(list(xref = "paper", x = 1, y = totale, xanchor = "right", yanchor = "bottom",
                                    showarrow = FALSE, font = list(color = "#c0392b"),
                                    text = paste0("Variance commune totale : ",
                                                  sub("\\.", ",", sprintf("%.1f", totale)), " %"))),
            showlegend = FALSE
          )
        p
      })

      output$variance_note <- renderUI({
        req(efa_result())
        efa <- efa_result()
        if (efa$rotation == "oblimin" && efa$nfactors > 1) {
          tags$div(class = "alert alert-light small mb-0",
                   "Avec une rotation oblique, les facteurs corrélés se partagent une partie de la variance : les barres (sommes des carrés des patrons) ne s'additionnent pas pour donner la variance commune totale. Seule la ligne rouge, moyenne des communalités, mesure la variance reproduite par l'ensemble des facteurs.")
        } else {
          tags$div(class = "alert alert-light small mb-0",
                   "Avec des facteurs non corrélés, les barres s'additionnent exactement pour donner la variance commune totale (ligne rouge), égale à la moyenne des communalités.")
        }
      })

      # Table des communalités
      output$communalities_table <- renderDT({
        req(efa_result())
        comm <- efa_result()$communality
        uniq <- efa_result()$uniquenesses
        
        # Préparer les données pour le tableau
        comm_data <- data.frame(
          Variable = names(comm),
          Communalité = round(comm, 3),
          Unicité = round(uniq, 3)
        )
        
        # Créer le tableau
        datatable(
          comm_data,
          options = list(pageLength = 15, scrollX = TRUE, dom = if (nrow(comm_data) > 15) "tip" else "t"),
          rownames = FALSE,
          caption = "Communalité h² = part de la variance de l'item reproduite par les facteurs retenus; unicité = 1 − h²."
        ) %>%
          formatRound(c("Communalité", "Unicité"), 3, dec.mark = ",") %>%
          formatStyle(
            "Communalité",
            background = styleColorBar(c(0, 1), "steelblue"),
            backgroundSize = "100% 90%",
            backgroundRepeat = "no-repeat",
            backgroundPosition = "center"
          )
      })
      
      # Résumé statistique
      output$data_summary <- renderDT({
        req(filtered_data())
        d <- filtered_data()
        num <- names(d)[vapply(d, is.numeric, logical(1))]
        validate(need(length(num) > 0, "Aucune variable numérique."))
        tab <- do.call(rbind, lapply(num, function(v) {
          x <- d[[v]]; x <- x[!is.na(x)]
          data.frame(Variable = v, n = length(x), Moyenne = mean(x), `Écart-type` = sd(x),
                     Min = min(x), Max = max(x), `Valeurs distinctes` = length(unique(x)),
                     Asymétrie = psych::skew(x), Aplatissement = psych::kurtosi(x), check.names = FALSE)
        }))
        datatable(tab, rownames = FALSE, options = list(dom = "tp", pageLength = 15, scrollX = TRUE),
                  caption = "Aplatissement : excès par rapport à la loi normale (0 = normal). Peu de valeurs distinctes signale un item ordinal (Likert).") %>%
          formatRound(c("Moyenne", "Écart-type", "Min", "Max", "Asymétrie", "Aplatissement"), 2, dec.mark = ",")
      })

      output$missingness_table <- renderDT({
        req(filtered_data())
        data <- filtered_data()
        selected_input <- input$efa_variables
        if (is.null(selected_input)) selected_input <- character(0)
        selected <- intersect(selected_input, names(data))
        if (length(selected) == 0) selected <- names(data)[vapply(data, is.numeric, logical(1))]
        if (length(selected) == 0) {
          return(datatable(data.frame(Message = "Aucune variable numérique sélectionnée."),
                           options = list(dom = "t"), rownames = FALSE))
        }
        out <- data.frame(
          Variable = selected,
          Manquantes = vapply(data[selected], function(x) sum(is.na(x)), numeric(1)),
          Pourcentage = vapply(data[selected], function(x) mean(is.na(x)) * 100, numeric(1)),
          stringsAsFactors = FALSE
        )
        datatable(out, options = list(pageLength = 15, dom = if (nrow(out) > 15) "tip" else "t"), rownames = FALSE) %>%
          formatRound("Pourcentage", digits = 1)
      })
      
      # Application du filtre
      observeEvent(input$apply_filter, {
        req(input$filter_var, raw_data())
        
        if (input$filter_var == "") {
          filtered_data(raw_data())
          return()
        }
        
        data <- raw_data()
        var <- data[[input$filter_var]]
        
        if (is.numeric(var) || is.integer(var)) {
          req(input$filter_range)
          filtered <- data[var >= input$filter_range[1] & var <= input$filter_range[2], ]
        } else if (is.factor(var) || is.character(var)) {
          req(input$filter_categories)
          filtered <- data[var %in% input$filter_categories, ]
        } else {
          filtered <- data
        }
        
        filtered_data(filtered)
        showNotification(paste0(nrow(filtered), " observations après filtrage"), type = "message")
        invalider_resultats()
      })
      
      # Réinitialisation du filtre
      observeEvent(input$reset_filter, {
        filtered_data(raw_data())
        invalider_resultats()
        updateSelectInput(session, "filter_var", selected = "")
        showNotification("Filtres réinitialisés", type = "message")
      })
      
      # Aide
      observeEvent(input$help_button, {
        showModal(modalDialog(
          title = "Aide : le parcours d'une AFE dans l'application",
          tags$div(
            tags$ol(
              tags$li(tags$b("Importer ou simuler "), "les données (barre latérale). Pour apprendre, commencez par « Données simulées » : l'onglet ",
                      tags$em("Aperçu des données ▸ Population simulée"), " montre alors la structure vraie."),
              tags$li(tags$b("Préparer "), "les items : inversez ceux qui sont formulés à l'envers, notez les valeurs manquantes."),
              tags$li(tags$b("Examiner "), "la matrice de corrélation et les distributions (", tags$em("Aperçu des données"), ")."),
              tags$li(tags$b("Diagnostiquer "), "la matrice et choisir un nombre de facteurs (", tags$em("Analyse factorielle ▸ Tests préliminaires"),
                      ") : Bartlett, déterminant, KMO, éboulis et analyse parallèle."),
              tags$li(tags$b("Estimer "), "la solution (", tags$em("Paramètres de l'analyse ▸ Exécuter l'analyse"), "), puis comparer les solutions voisines."),
              tags$li(tags$b("Interpréter "), "dans l'ordre des sous-onglets de l'onglet ", tags$em("Interprétation"), "."),
              tags$li(tags$b("Éprouver "), "la solution : fiabilité, bootstrap et, pour des données simulées, validation CFA sur un nouvel échantillon.")
            ),
            tags$p(class = "text-muted", "Les seuils usuels (0,30 ou 0,40 pour une charge, 0,60 pour le KMO) sont des repères de lecture, pas des règles de décision.")
          ),
          easyClose = TRUE,
          footer = modalButton("Fermer")
        ))
      })
      
      # Exécution de l'analyse factorielle avec lavaan
      observeEvent(input$run_efa, {
        req(filtered_data(), input$efa_variables, length(input$efa_variables) >= 2)
        bootstrap_results(NULL)
        cfa_result(NULL)

        withProgress(message = "Calcul de l'analyse factorielle...", {
          tryCatch({
            data_subset <- filtered_data()[, input$efa_variables, drop = FALSE]
            data_complete <- na.omit(data_subset)
            if (nrow(data_complete) < nrow(data_subset)) {
              showNotification(
                paste0(nrow(data_subset) - nrow(data_complete),
                       " observations avec valeurs manquantes ont été écartées (cas complets)."),
                type = "warning")
            }
            if (input$efa_nfactors >= ncol(data_complete)) {
              stop("Le nombre de facteurs doit être inférieur au nombre d'items.")
            }
            vp <- eigen(cor(data_complete), only.values = TRUE)$values
            if (max(vp) / min(vp) > 1000) {
              showNotification("Attention : forte multicolinéarité; l'estimation peut être instable.",
                               type = "warning", duration = 10)
            }

            efa_result(ajuster_afe(data_complete, input$efa_nfactors,
                                   input$dynamic_rotation, input$efa_method,
                                   n_total = nrow(data_subset)))

            # Solutions voisines : de 1 facteur jusqu'à 5 (ou p − 1).
            incProgress(0.6, message = "Comparaison des solutions voisines...")
            # Au plus 5 facteurs, et seulement des modèles identifiés (ddl > 0) :
            # ddl = [(p − q)² − (p + q)] / 2.
            p_items <- ncol(data_complete)
            q_ok <- which(sapply(1:5, function(q) ((p_items - q)^2 - (p_items + q)) / 2 > 0))
            max_factors <- max(1, max(c(1, q_ok)))
            voisines <- lapply(seq_len(max_factors), function(q) {
              tryCatch(ajuster_afe(data_complete, q, input$dynamic_rotation, input$efa_method),
                       error = function(e) NULL)
            })
            model_comparison(list(solutions = voisines, nfactors = seq_len(max_factors)))

            showTab(inputId = "main_tabs", target = "Interprétation")
            showTab(inputId = "main_tabs", target = "Bootstrap")
            if (!is.null(simulation_spec())) {
              showTab(inputId = "main_tabs", target = "Validation CFA")
            }
            showNotification("Analyse factorielle terminée", type = "message")
          }, error = function(e) {
            showNotification(paste("Erreur lors de l'analyse factorielle :", e$message),
                             type = "error", duration = 10)
          })
        })
      })

      observeEvent(input$import_method, {
        if(input$import_method == "example") {
          # Le scénario affiché au départ est réellement le scénario clair.
          updateNumericInput(session, "sim_nobs", value = 300)
          updateNumericInput(session, "sim_nfactors", value = 2)
          updateNumericInput(session, "sim_nvars", value = 4)
          updateSliderInput(session, "sim_loading", value = 0.70)
          updateSliderInput(session, "sim_cross_loading", value = 0.05)
          updateSliderInput(session, "sim_cross_proportion", value = 0.25)
          updateSliderInput(session, "sim_error", value = 1.00)
          updateCheckboxInput(session, "sim_correlated", value = TRUE)
        }
      })

      # Quatre situations délibérément contrastées. Les réglages restent
      # modifiables après le choix du scénario.
      observeEvent(input$sim_scenario, ignoreInit = TRUE, {
        presets <- list(
          clear = list(n = 300, nf = 2, nv = 4, main = .70, cross = .05, cross_prop = .25, noise = 1, corr = .30),
          correlated = list(n = 300, nf = 2, nv = 4, main = .65, cross = .05, cross_prop = .25, noise = 1, corr = .55),
          cross = list(n = 300, nf = 2, nv = 4, main = .60, cross = .30, cross_prop = .50, noise = 1, corr = .30),
          fragile = list(n = 100, nf = 3, nv = 3, main = .50, cross = .15, cross_prop = .35, noise = 1.00, corr = .30)
        )
        p <- presets[[input$sim_scenario]]
        updateNumericInput(session, "sim_nobs", value = p$n)
        updateNumericInput(session, "sim_nfactors", value = p$nf)
        updateNumericInput(session, "sim_nvars", value = p$nv)
        updateSliderInput(session, "sim_loading", value = p$main)
        updateSliderInput(session, "sim_cross_loading", value = p$cross)
        updateSliderInput(session, "sim_cross_proportion", value = p$cross_prop)
        updateSliderInput(session, "sim_error", value = p$noise)
        updateCheckboxInput(session, "sim_correlated", value = TRUE)
        updateSliderInput(session, "sim_factor_cor", value = p$corr)
      })
      
      # Gestion des données d'exemple
      observeEvent(input$generate_example, {
        # Nettoyer les objets existants avant de générer de nouvelles données
        cleanup_reactive_objects()
        withProgress(message = "Génération des données simulées...", {
          # Récupérer les paramètres 
          n <- input$sim_nobs
          num_factors <- input$sim_nfactors
          vars_per_factor <- input$sim_nvars
          loading_strength <- input$sim_loading
          cross_loading <- input$sim_cross_loading
          cross_proportion <- input$sim_cross_proportion
          error_var <- input$sim_error
          correlated <- input$sim_correlated
          factor_cor <- input$sim_factor_cor

          # La même fonction sert ensuite à créer l'échantillon CFA : la
          # population simulée reste donc strictement identique.
          example_data <- tryCatch(
            generate_efa_simulation(
              n, num_factors, vars_per_factor, loading_strength,
              cross_loading, cross_proportion, error_var, correlated, factor_cor, input$sim_seed
            ),
            error = function(e) e
          )
          if (inherits(example_data, "error")) {
            showNotification(
              paste0(example_data$message, " Réduisez les charges secondaires, la corrélation ou le nombre de facteurs."),
              type = "error", duration = 12
            )
            return()
          }
          
          # Conserver les paramètres de la population. Ils seront réutilisés
          # tels quels pour le nouvel échantillon de validation CFA.
          simulation_spec(list(
            n = n, num_factors = num_factors, vars_per_factor = vars_per_factor,
            loading_strength = loading_strength, cross_loading = cross_loading,
            cross_proportion = cross_proportion,
            error_multiplier = error_var, correlated = correlated,
            factor_cor = factor_cor, discovery_seed = input$sim_seed,
            loading_pattern = attr(example_data, "loading_pattern")
          ))
          validation_data(NULL)
          cfa_result(NULL)
          updateNumericInput(session, "cfa_validation_seed",
                             value = if (input$sim_seed >= 999999) 1 else input$sim_seed + 1)
          hideTab(inputId = "main_tabs", target = "Validation CFA")

          # Stocker les données
          original_data(example_data)
          raw_data(example_data)
          filtered_data(example_data)
          preparation_log(preparation_log()[0, ])
          
          # Mettre à jour les sélecteurs
          updateSelectInput(session, "filter_var",
                            choices = c("Aucun" = "", names(example_data)))
          
          updateSelectInput(session, "efa_variables",
                            choices = names(example_data),
                            selected = names(example_data))
          updateSelectizeInput(session, "reverse_items", choices = names(example_data),
                               selected = character(0))
          
          # Mettre à jour le nombre de facteurs recommandé
          updateNumericInput(session, "efa_nfactors", 
                             value = num_factors)
          
          # Afficher l'onglet d'analyse
          showTab(inputId = "main_tabs", target = "Analyse factorielle")
          
          # Ajouter une notification avec un résumé des données simulées
          showNotification(
            HTML(paste0(
              "<b>Données simulées générées</b><br/>",
              "Observations: ", n, "<br/>",
              "Structure: ", num_factors, " facteurs avec ", vars_per_factor, " variables chacun<br/>",
              "Charges principales: ", loading_strength, "<br/>",
              "Facteurs ", ifelse(correlated, "corrélés", "non corrélés"), "<br/>",
              "Graine reproductible : ", input$sim_seed
            )),
            type = "message",
            duration = 5
          )
        })
      })
      
      
      # Résumé du fichier - Chaque élément sur une ligne séparée
      output$file_summary_compact <- renderUI({
        req(input$data_file)
        
        # Déterminer le type de fichier
        file_ext <- tools::file_ext(input$data_file$name)
        file_size <- round(input$data_file$size / 1024, 1)
        
        # Style CSS pour les labels et les valeurs
        label_style <- "font-weight: bold; margin-right: 5px;"
        
        # Création d'un résumé avec chaque élément sur une ligne séparée
        div(
          div(
            span(style = label_style, "Fichier:"), 
            span(input$data_file$name)
          ),
          div(
            span(style = label_style, "Type:"), 
            span(toupper(file_ext))
          ),
          div(
            span(style = label_style, "Taille:"), 
            span(paste0(file_size, " ko"))
          )
        )
      })
      
      # Affichage de l'aperçu des données - limité à 4 décimales
      output$data_preview <- renderDT({
        if(is.null(filtered_data())) {
          return(datatable(data.frame(Message = "Aucune donnée chargée. Veuillez importer un fichier."),
                           options = list(dom = 't'),
                           rownames = FALSE))
        }
        
        # Créer une copie des données pour formater les nombres
        formatted_data <- as.data.frame(filtered_data())
        
        # Formater les colonnes numériques pour afficher au maximum 4 décimales
        for(col in names(formatted_data)) {
          if(is.numeric(formatted_data[[col]])) {
            formatted_data[[col]] <- round(formatted_data[[col]], 4)
          }
        }
        
        datatable(
          formatted_data,
          options = list(
            pageLength = 10,
            scrollX = TRUE,
            dom = 'ftip'
          ),
          rownames = FALSE
        )
      })
      
      # Contrôles de filtre dynamiques
      output$filter_controls <- renderUI({
        req(input$filter_var, raw_data())
        
        if (input$filter_var == "") {
          return(NULL)
        }
        
        var <- raw_data()[[input$filter_var]]
        
        if (is.numeric(var) || is.integer(var)) {
          min_val <- min(var, na.rm = TRUE)
          max_val <- max(var, na.rm = TRUE)
          
          sliderInput("filter_range", 
                      label = paste("Plage de valeurs pour", input$filter_var),
                      min = min_val, 
                      max = max_val,
                      value = c(min_val, max_val),
                      step = (max_val - min_val) / 100)
        } else if (is.factor(var) || is.character(var)) {
          unique_vals <- unique(var)
          checkboxGroupInput("filter_categories",
                             label = paste("Catégories pour", input$filter_var),
                             choices = unique_vals,
                             selected = unique_vals)
        } else {
          return(tags$div("Type de variable non supporté pour le filtrage"))
        }
      })

      output$preparation_summary <- renderUI({
        log <- preparation_log()
        if (is.null(log) || nrow(log) == 0) {
          return(tags$p(class = "text-muted small mt-2 mb-0",
                        "Aucune transformation appliquée."))
        }
        tags$div(
          class = "alert alert-light small mt-3 mb-0",
          tags$b("Transformations consignées : "),
          paste(paste0(log$Item, " (", log$Transformation, ")"), collapse = "; ")
        )
      })

      observeEvent(input$reverse_selected_items, {
        req(raw_data(), input$reverse_items, length(input$reverse_items) > 0)
        data <- raw_data()
        items <- intersect(input$reverse_items, names(data))
        items <- items[vapply(data[items], is.numeric, logical(1))]
        if (length(items) == 0) {
          showNotification("Sélectionnez au moins un item numérique à inverser.",
                           type = "warning")
          return()
        }
        log_rows <- vector("list", length(items))
        for (i in seq_along(items)) {
          item <- items[i]
          values <- data[[item]]
          lo <- if (isTRUE(is.finite(input$reverse_min))) input$reverse_min else min(values, na.rm = TRUE)
          hi <- if (isTRUE(is.finite(input$reverse_max))) input$reverse_max else max(values, na.rm = TRUE)
          if (!is.finite(lo) || !is.finite(hi) || lo == hi) {
            showNotification(paste0("L'item ", item, " ne présente pas une échelle numérique inversable."),
                             type = "error")
            return()
          }
          data[[item]] <- lo + hi - values
          log_rows[[i]] <- data.frame(
            Item = item, Ancien_min = lo, Ancien_max = hi,
            Transformation = paste0(lo, " + ", hi, " − x"),
            stringsAsFactors = FALSE
          )
        }
        raw_data(data)
        filtered_data(data)
        preparation_log(rbind(preparation_log(), do.call(rbind, log_rows)))
        invalider_resultats()
        updateSelectizeInput(session, "reverse_items", selected = character(0))
        showNotification(
          "Items inversés. Le filtre a été réinitialisé; réappliquez-le au besoin.",
          type = "message", duration = 7
        )
      })

      observeEvent(input$reset_preparation, {
        req(original_data())
        raw_data(original_data())
        invalider_resultats()
        filtered_data(original_data())
        preparation_log(preparation_log()[0, ])
        updateSelectizeInput(session, "reverse_items", selected = character(0))
        updateSelectInput(session, "filter_var", selected = "")
        showNotification("Données originales restaurées; les filtres et transformations ont été retirés.",
                         type = "message", duration = 7)
      })
      
      # Output pour afficher la comparaison des modèles
      # Comparaison des solutions voisines : ajustement global ET lisibilité.
      # Aucun indice ne décide seul; le tableau sert à voir ce que l'on gagne
      # (ou perd) en ajoutant un facteur.
      output$model_comparison_output <- renderDT({
        req(model_comparison())
        mc <- model_comparison()
        retenu <- efa_result()$nfactors
        lignes <- lapply(seq_along(mc$solutions), function(k) {
          sol <- mc$solutions[[k]]
          if (is.null(sol)) {
            return(data.frame(Facteurs = mc$nfactors[k], `χ²` = "—", ddl = "—", p = "—",
                              CFI = "—", TLI = "—", RMSEA = "—", SRMR = "—", BIC = "—",
                              `Items nets` = "—", `Facteurs bien définis` = "—", `Croisés / faibles` = "non convergé",
                              check.names = FALSE))
          }
          f <- sol$fit_indices
          g <- function(n) if (!is.null(f) && n %in% names(f)) unname(f[n]) else NA
          st <- summarise_simple_structure(sol$loadings, sol$communality,
                                           threshold = input$loadings_threshold)
          data.frame(
            Facteurs = mc$nfactors[k],
            `χ²` = fr(g("chisq"), 1), ddl = as.integer(g("df")), p = fr_p(g("pvalue")),
            CFI = fr(g("cfi"), 3), TLI = fr(g("tli"), 3), RMSEA = fr(g("rmsea"), 3),
            SRMR = fr(g("srmr"), 3), BIC = fr(g("bic"), 1),
            `Items nets` = sum(st$Statut == "structure simple"),
            `Facteurs bien définis` = paste0(sum(table(factor(st$Facteur_principal[st$Statut == "structure simple"],
                                                              levels = paste0("F", seq_len(mc$nfactors[k])))) >= 3),
                                             " sur ", mc$nfactors[k]),
            `Croisés / faibles` = paste0(sum(st$Statut == "charge croisée"), " / ",
                                         sum(st$Statut == "charge faible")),
            check.names = FALSE)
        })
        tab <- do.call(rbind, lignes)
        datatable(tab, rownames = FALSE,
                  options = list(dom = "t", ordering = FALSE,
                                 columnDefs = list(list(className = "dt-right", targets = "_all"))),
                  caption = paste0("Ligne surlignée : la solution actuellement retenue (", retenu,
                                   " facteur(s)). « Items nets » : charge principale au-dessus du seuil de lecture de l'onglet Interprétation, sans seconde charge notable. « Facteurs bien définis » : facteurs ayant au moins trois items nets.")) %>%
          formatStyle("Facteurs", target = "row",
                      backgroundColor = styleEqual(retenu, "#e8f4fd"))
      })

      output$download_simulated_data <- downloadHandler(
        filename = function() {
          ext <- input$download_format
          paste("donnees_simulees_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext, sep = "")
        },
        content = function(file) {
          req(raw_data())
          data <- raw_data()
          
          # Exporter selon le format sélectionné
          if (input$download_format == "csv") {
            write.csv(data, file, row.names = FALSE, quote = TRUE)
          } else if (input$download_format == "xlsx") {
            writexl::write_xlsx(data, path = file)
          } else if (input$download_format == "sav") {
            haven::write_sav(data, file)
          } else if (input$download_format == "dta") {
            haven::write_dta(data, file)
          }
        }
      )
      
      # Matrice de corrélation : items retenus pour l'AFE (ou toutes les
      # variables numériques si aucune sélection), cas complets par paire.
      cor_items <- reactive({
        req(filtered_data())
        d <- filtered_data()
        v <- intersect(input$efa_variables, names(d))
        if (length(v) < 2) v <- names(d)[vapply(d, is.numeric, logical(1))]
        validate(need(length(v) >= 2, "Il faut au moins deux variables numériques."))
        cor(d[, v, drop = FALSE], use = "pairwise.complete.obs")
      })

      output$correlation_output <- renderUI({
        req(filtered_data())
        if (identical(input$corr_display_type, "plot")) {
          n <- ncol(cor_items())
          plotOutput("correlation_plot", height = paste0(max(380, min(900, 34 * n + 120)), "px"))
        } else DTOutput("correlation_text")
      })

      output$correlation_plot <- renderPlot({
        R <- cor_items(); v <- colnames(R); n <- length(v)
        d <- expand.grid(i = seq_len(n), j = seq_len(n))
        d$r <- R[cbind(d$i, d$j)]
        d <- d[d$i > d$j, ]
        d$x <- factor(v[d$j], levels = v); d$y <- factor(v[d$i], levels = rev(v))
        ggplot(d, aes(x, y, fill = r)) +
          geom_tile(colour = "white") +
          { if (n <= 25) geom_text(aes(label = sub("^-", "−", sub("\\.", ",", sprintf("%.2f", r)))),
                                   size = if (n <= 12) 3.8 else 2.8) } +
          scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", limits = c(-1, 1),
                               labels = function(x) sub("\\.", ",", format(x)), name = "r") +
          labs(x = NULL, y = NULL, title = paste0("Corrélations entre les ", n, " items (triangle inférieur)")) +
          theme_minimal(base_size = 13) +
          theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
      })

      output$correlation_text <- renderDT({
        R <- round(cor_items(), 3)
        R[upper.tri(R)] <- NA
        datatable(as.data.frame(R), rownames = TRUE,
                  options = list(dom = "t", pageLength = 500, ordering = FALSE, scrollX = TRUE)) %>%
          formatRound(colnames(R), 3, dec.mark = ",")
      })

      # Mettre à jour automatiquement les tests de normalité quand les données changent
      observe({
        req(filtered_data())
        
        # Appliquer les tests de normalité automatiquement à toutes les variables numériques
        data <- filtered_data()
        numeric_vars <- names(data)[sapply(data, is.numeric)]
        
        # Limiter à 20 variables pour des raisons de performance
        if (length(numeric_vars) > 20) {
          numeric_vars <- numeric_vars[1:20]
        }
        
        # Créer un tableau pour stocker les résultats
        results <- data.frame(
          Variable = character(),
          Shapiro_Stat = numeric(),
          Shapiro_p = numeric(),
          Normal = character(),
          stringsAsFactors = FALSE
        )
        
        # Exécuter les tests pour chaque variable
        for (var_name in numeric_vars) {
          var_data <- data[[var_name]]
          
          # Vérifier que la variable a suffisamment de données
          if (sum(!is.na(var_data)) >= 3) {
            # Limiter à 5000 observations pour le test de Shapiro (limitation de la méthode)
            if (length(var_data) > 5000) {
              var_data <- sample(var_data[!is.na(var_data)], 5000)
            } else {
              var_data <- var_data[!is.na(var_data)]
            }
            
            # Test de Shapiro-Wilk
            tryCatch({
              shapiro_test <- shapiro.test(var_data)
              
              # Ajouter les résultats
              results <- rbind(results, data.frame(
                Variable = var_name,
                Shapiro_Stat = round(shapiro_test$statistic, 4),
                Shapiro_p = round(shapiro_test$p.value, 4),
                Normal = ifelse(shapiro_test$p.value > 0.05, "Oui", "Non"),
                stringsAsFactors = FALSE
              ))
            }, error = function(e) {
              # En cas d'erreur, ajouter une ligne avec message d'erreur
              results <- rbind(results, data.frame(
                Variable = var_name,
                Shapiro_Stat = NA,
                Shapiro_p = NA,
                Normal = "Erreur de test",
                stringsAsFactors = FALSE
              ))
            })
          }
        }
        
        # Stocker les résultats dans la réactive
        normality_results_data(results)
      })
      
      # Afficher les résultats des tests de normalité dans un tableau interactif
      output$normality_results_table <- renderDT({
        req(normality_results_data())
        results <- normality_results_data()
        
        if (nrow(results) == 0) {
          return(datatable(data.frame(Message = "Aucune variable numérique trouvée pour les tests de normalité"),
                           options = list(dom = 't'),
                           rownames = FALSE))
        }
        
        # Préparer un tableau formaté pour l'affichage
        results_table <- data.frame(
          Variable = results$Variable,
          "Stat. Shapiro" = results$Shapiro_Stat,
          "p-value" = results$Shapiro_p,
          "Distribution normale" = ifelse(results$Normal == "Oui", 
                                          "Oui (p > 0,05)", 
                                          ifelse(results$Normal == "Non", 
                                                 "Non (p < 0,05)", 
                                                 results$Normal)),
          check.names = FALSE
        )
        
        # Créer le tableau interactif avec coloration des p-values
        datatable(
          results_table,
          options = list(
            pageLength = 15,
            dom = 'ftip'
          ),
          rownames = FALSE
        ) %>%
          formatRound(c("Stat. Shapiro", "p-value"), 3, dec.mark = ",") %>%
          formatStyle(
            "p-value",
            color = styleInterval(0.05, c("#b03a2e", "#1e6b2f")),
            fontWeight = styleInterval(0.05, c("bold", "normal"))
          )
      })
      
      # Mettre à jour le selectInput pour la matrice de visualisations
      observe({
        req(filtered_data())
        data <- filtered_data()
        
        # Mettre à jour les sélecteurs avec les variables numériques
        numeric_vars <- names(data)[sapply(data, is.numeric)]
        updateSelectInput(session, "viz_variables", choices = numeric_vars)
      })
      
      # Créer la matrice de visualisations
      output$matrix_plot <- renderPlot({
        req(filtered_data(), input$viz_variables)
        
        # Vérifier qu'il y a des variables sélectionnées
        if (length(input$viz_variables) == 0) {
          plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
          text(0, 0, "Sélectionnez au moins une variable à visualiser", cex = 1.5)
          return()
        }
        
        # Limiter à 8 variables maximum
        vars_to_plot <- input$viz_variables
        if (length(vars_to_plot) > 8) {
          vars_to_plot <- vars_to_plot[1:8]
        }
        
        # Récupérer les données pour ces variables
        data <- filtered_data()
        plot_data <- data[, vars_to_plot, drop = FALSE]
        
        # Déterminer la disposition de la matrice (nombre de lignes et colonnes)
        n_vars <- length(vars_to_plot)
        n_cols <- min(4, n_vars)  # Maximum 4 colonnes
        n_rows <- ceiling(n_vars / n_cols)
        
        # Configurer une mise en page en matrice
        par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
        
        # Créer un graphique pour chaque variable
        for (var in vars_to_plot) {
          var_data <- data[[var]]
          var_data <- var_data[!is.na(var_data)]
          
          if (length(var_data) == 0) {
            plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
            text(0, 0, "Pas de données", cex = 1.2)
            next
          }
          
          if (input$matrix_plot_type == "hist") {
            # Histogramme
            h <- hist(var_data, plot = FALSE, breaks = min(30, max(10, round(length(var_data)/30))))
            
            hist(var_data, 
                 breaks = h$breaks,
                 main = var,
                 xlab = "",
                 ylab = "Fréquence",
                 col = "skyblue",
                 border = "white",
                 prob = input$matrix_density)
            
            if (input$matrix_density) {
              d <- density(var_data, na.rm = TRUE)
              lines(d, col = "red", lwd = 2)
              x <- seq(min(var_data), max(var_data), length.out = 100)
              lines(x, dnorm(x, mean = mean(var_data), sd = sd(var_data)), 
                    col = "blue", lty = 2, lwd = 2)
            }
            
          } else if (input$matrix_plot_type == "box") {
            # Boxplot
            boxplot(var_data, 
                    main = var,
                    xlab = "",
                    ylab = "",
                    col = "lightgreen",
                    horizontal = TRUE)
            
          } else if (input$matrix_plot_type == "qq") {
            # Q-Q plot
            qqnorm(var_data, 
                   main = var,
                   col = "darkblue",
                   pch = 19)
            qqline(var_data, col = "red", lwd = 2)
          }
        }
        
        # Ajouter une légende commune si c'est un histogramme avec densité
        if (input$matrix_plot_type == "hist" && input$matrix_density) {
          # Ajouter une légende dans la marge externe supérieure
          par(fig = c(0, 1, 0, 1), oma = c(0, 0, 3, 0), mar = c(0, 0, 0, 0), new = TRUE)
          plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
          legend("top", 
                 legend = c("Densité observée", "Distribution normale"), 
                 col = c("red", "blue"), 
                 lty = c(1, 2), 
                 lwd = 2,
                 horiz = TRUE,
                 bty = "n")
        }
      })
      
      # Génération indépendante de l'échantillon CFA à partir de la même
      # population simulée. Une graine identique reproduirait l'échantillon
      # d'exploration; elle est donc refusée.
      observeEvent(input$generate_validation_sample, {
        spec <- simulation_spec()
        req(spec)
        if (input$cfa_validation_seed == spec$discovery_seed) {
          showNotification(
            "Choisissez une graine différente : sinon l'échantillon de validation serait identique à celui de l'AFE.",
            type = "error", duration = 8
          )
          return()
        }

        new_sample <- tryCatch(
          generate_efa_simulation(
            n = spec$n,
            num_factors = spec$num_factors,
            vars_per_factor = spec$vars_per_factor,
            loading_strength = spec$loading_strength,
            cross_loading = spec$cross_loading,
            cross_proportion = spec$cross_proportion,
            error_multiplier = spec$error_multiplier,
            correlated = spec$correlated,
            factor_cor = spec$factor_cor,
            seed = input$cfa_validation_seed,
            loading_pattern = spec$loading_pattern
          ),
          error = function(e) e
        )
        if (inherits(new_sample, "error")) {
          showNotification(new_sample$message, type = "error", duration = 10)
          return()
        }
        validation_data(new_sample)
        cfa_result(NULL)
        showNotification(
          "Nouvel échantillon généré. Vous pouvez maintenant évaluer le modèle CFA.",
          type = "message", duration = 6
        )
      })

      output$validation_sample_info <- renderUI({
        spec <- simulation_spec()
        if (is.null(validation_data()) || is.null(spec)) {
          return(tags$p(class = "text-muted mt-3 mb-0",
                        "Aucun échantillon de validation n'a encore été généré."))
        }
        tags$div(
          class = "alert alert-success mt-3 mb-0",
          tags$b("Échantillon de validation prêt. "),
          nrow(validation_data()), " observations, ", spec$num_factors,
          " facteurs et ", spec$vars_per_factor,
          " items par facteur; graine ", input$cfa_validation_seed, "."
        )
      })

      output$cfa_crossloading_warning <- renderUI({
        req(efa_result())
        loadings <- efa_result()$loadings
        threshold <- input$cfa_threshold
        n_hits <- rowSums(abs(loadings) >= threshold)
        ambiguous <- rownames(loadings)[n_hits > 1]
        if (length(ambiguous) == 0) return(NULL)
        if (identical(input$cfa_crossloading_rule, "all")) {
          tags$div(class = "alert alert-warning small",
                   tags$b("Charges croisées incluses. "),
                   "Les items suivants seront liés à plus d’un facteur : ",
                   paste(ambiguous, collapse = ", "),
                   ". Cette décision rend le modèle candidat partiellement exploratoire et doit être justifiée.")
        } else {
          tags$div(class = "alert alert-info small",
                   tags$b("Items ambigus repérés. "),
                   paste(ambiguous, collapse = ", "),
                   ". La syntaxe proposée les affecte seulement à leur facteur principal; relisez leur contenu avant de confirmer ce choix.")
        }
      })

      # ---------------------------------------------------------------
      # Population simulée : charges standardisées vraies, puis, après une
      # AFE comparable, les charges estimées alignées (ordre et signe).
      # ---------------------------------------------------------------
      population <- reactive({
        spec <- simulation_spec()
        req(spec)
        phi <- if (spec$correlated && spec$num_factors > 1) {
          m <- matrix(spec$factor_cor, spec$num_factors, spec$num_factors); diag(m) <- 1; m
        } else diag(spec$num_factors)
        population_simulee(spec$loading_pattern, phi, spec$error_multiplier)
      })

      estimation_alignee <- reactive({
        pop <- population(); efa <- efa_result()
        req(efa)
        if (efa$nfactors != ncol(pop$loadings) || !all(efa$variables %in% rownames(pop$loadings))) return(NULL)
        est <- align_bootstrap_loadings(pop$loadings[efa$variables, , drop = FALSE], efa$loadings)
        list(L = est, h2 = efa$communality, phi = efa$factor_correlation)
      })

      output$population_info <- renderUI({
        if (is.null(simulation_spec())) {
          return(tags$div(class = "alert alert-secondary",
                          "Cet onglet n'a de contenu que pour des données simulées : avec des données réelles, la structure vraie est inconnue."))
        }
        spec <- simulation_spec()
        tags$div(class = "alert alert-info small",
          tags$p(tags$b("Charges standardisées de population. "),
                 "Ce sont les valeurs qu'une AFE retrouverait avec un échantillon infini. Elles diffèrent du curseur « Force des charges principales » dès que le multiplicateur de bruit n'est pas 1 (ici ",
                 fr(spec$error_multiplier, 2), ") : un bruit réduit rend les items plus fiables, donc leurs charges standardisées plus fortes."),
          tags$p(class = "mb-0",
                 if (is.null(efa_result())) "Exécutez une AFE avec le même nombre de facteurs pour voir, à côté, les charges estimées."
                 else "Les colonnes « estimé » proviennent de l'AFE affichée, après alignement de l'ordre et du signe des facteurs.")
        )
      })

      output$population_table <- renderDT({
        pop <- population()
        q <- ncol(pop$loadings)
        tab <- data.frame(Item = rownames(pop$loadings), check.names = FALSE)
        for (k in seq_len(q)) tab[[paste0("F", k, " population")]] <- pop$loadings[, k]
        tab[["h² population"]] <- pop$communality
        est <- tryCatch(estimation_alignee(), error = function(e) NULL)
        cap <- NULL
        if (!is.null(est)) {
          for (k in seq_len(q)) tab[[paste0("F", k, " estimé")]] <- est$L[tab$Item, k]
          tab[["h² estimé"]] <- est$h2[tab$Item]
          cong <- sapply(seq_len(q), function(k) congruence_tucker(pop$loadings[rownames(est$L), k], est$L[, k]))
          cap <- paste0("Congruence de Tucker entre facteurs vrais et estimés : ",
                        paste0("F", seq_len(q), " = ", fr(cong, 3), collapse = "; "),
                        ". Repères usuels : ≥ 0,95 quasi identiques; 0,85 à 0,94 semblables.")
        }
        num <- names(tab)[-1]
        datatable(tab, rownames = FALSE, caption = cap,
                  options = list(dom = "t", pageLength = 50, ordering = FALSE, scrollX = TRUE)) %>%
          formatRound(num, 3, dec.mark = ",")
      })

      output$population_phi <- renderUI({
        pop <- population()
        q <- ncol(pop$phi)
        if (q < 2) return(NULL)
        est <- tryCatch(estimation_alignee(), error = function(e) NULL)
        r_pop <- pop$phi[upper.tri(pop$phi)]
        txt <- paste0("Corrélation entre facteurs dans la population : ", fr(unique(r_pop), 2), ".")
        if (!is.null(est)) {
          r_est <- est$phi[upper.tri(est$phi)]
          txt <- paste0(txt, " Estimation (", efa_result()$rotation, ") : ",
                        paste(fr(r_est, 2), collapse = "; "), ".",
                        if (efa_result()$rotation == "varimax") " Varimax force cette corrélation à 0, quelle que soit la population." else "")
        }
        tags$p(txt)
      })

      # Fonction pour générer le modèle CFA à partir de l'EFA
      observeEvent(input$generate_cfa, {
        req(efa_result())
        if (is.null(validation_data())) {
          showNotification(
            "Générez d'abord le nouvel échantillon dans l'étape 1.",
            type = "warning", duration = 6
          )
          return()
        }
        
        # Générer la syntaxe CFA
        cfa_syntax <- convert_efa_to_cfa(
          efa_result(), input$cfa_threshold, input$cfa_crossloading_rule
        )
        
        # Mettre à jour le champ de texte avec la syntaxe
        updateTextAreaInput(session, "cfa_syntax", value = cfa_syntax)
      })
      
      # Fonction pour exécuter l'analyse confirmatoire
      observeEvent(input$run_cfa, {
        req(input$cfa_syntax, validation_data())
        
        withProgress(message = "Exécution de l'analyse confirmatoire...", {
          tryCatch({
            # Obtenir la syntaxe et les données
            cfa_model_syntax <- input$cfa_syntax
            data <- validation_data()
            
            # Sélectionner seulement les variables utilisées dans le modèle
            used_vars <- lavaan::lavParseModelString(cfa_model_syntax)
            used_vars <- unique(c(used_vars$lhs, used_vars$rhs))
            used_vars <- intersect(used_vars, names(data))
            
            # Vérifier la factoriabilité avant l'analyse
            if (length(used_vars) > 1) {
              data_subset <- data[, used_vars, drop = FALSE]
              data_subset <- na.omit(data_subset)
              
              # Calculer le déterminant de la matrice de corrélation
              cor_matrix <- cor(data_subset)
              det_cor <- det(cor_matrix)
              
              # Avertir si le déterminant est trop proche de zéro
              if (det_cor < 0.00001) {
                showNotification(
                  HTML(paste0(
                    "<b>Attention: Forte dépendance linéaire détectée</b><br/>",
                    "Le déterminant de la matrice de corrélation est très proche de zéro (", round(det_cor, 8), ").<br/>",
                    "Cela peut indiquer que certaines variables sont presque parfaitement corrélées, ",
                    "ce qui peut causer des problèmes d'estimation."
                  )),
                  type = "warning",
                  duration = 10
                )
              }
            }
            
            # Exécuter le modèle CFA avec la méthode sélectionnée
            cfa_model <- lavaan::cfa(
              model = cfa_model_syntax,
              data = data,
 #             missing = "fiml",
              std.lv = TRUE,
              estimator = input$cfa_method,  # Utiliser la méthode sélectionnée
              se = ifelse(input$cfa_method %in% c("MLR", "MLM", "MLMV"), "robust", "standard"),
              control = list(
                iter.max = 10000,
                optimizer = "NLMINB",
                ridge = 0.1
              )
            )
            
            # Vérifier si des problèmes existent après ajustement
            fit_measures <- fitMeasures(cfa_model)
            
            # Vérifier les variances négatives
            param_est <- parameterEstimates(cfa_model)
            var_params <- param_est[param_est$op == "~~" & param_est$lhs == param_est$rhs, ]
            
            # S'il y a des variances négatives
            if (any(var_params$est < 0)) {
              negative_vars <- var_params$lhs[var_params$est < 0]
              
              showNotification(
                HTML(paste0(
                  "<b>Variances négatives détectées</b><br/>",
                  "Variables: ", paste(negative_vars, collapse=", "), "<br/>",
                  "Cela peut indiquer un sur-ajustement ou un problème de spécification.<br/>",
                  "Suggestions:<br/>",
                  "- Simplifier votre modèle<br/>",
                  "- Essayer une méthode d'estimation différente (ex: ULS)<br/>",
                  "- Ajouter des contraintes de non-négativité"
                )),
                type = "warning",
                duration = 10
              )
            }
            
            # Stocker les résultats
            cfa_result(cfa_model)
            
            # Notification de succès
            showNotification(paste0("Analyse confirmatoire terminée (Méthode: ", input$cfa_method, ")"), type = "message")
            
          }, error = function(e) {
            # En cas d'erreur, essayer avec ULS
            tryCatch({
              showNotification(
                paste0("Échec avec la méthode ", input$cfa_method, ". Tentative avec ULS..."),
                type = "warning",
                duration = 5
              )
              
              # Obtenir la syntaxe et les données
              cfa_model_syntax <- input$cfa_syntax
              data <- validation_data()
              
              # Essayer avec ULS qui est moins sensible aux problèmes numériques
              cfa_model <- lavaan::cfa(
                model = cfa_model_syntax,
                data = data,
 #               missing = "fiml",
                std.lv = TRUE,
                estimator = "ULS",
                control = list(
                  iter.max = 15000,
                  optimizer = "NLMINB",
                  ridge = 0.5
                )
              )
              
              # Stocker les résultats
              cfa_result(cfa_model)
              
              showNotification(
                "Analyse confirmatoire terminée avec méthode alternative (ULS)",
                type = "warning",
                duration = 7
              )
              
            }, error = function(e2) {
              # Si même l'ULS échoue
              showNotification(
                HTML(paste0(
                  "<b>Erreur lors de l'analyse confirmatoire</b><br/>",
                  e$message, "<br/><br/>",
                  "Votre modèle pourrait avoir des problèmes d'identification ou de spécification."
                )),
                type = "error",
                duration = 15
              )
            })
          })
        })
      })
      
      output$cfa_loadings_table <- renderDT({
        req(cfa_result())
        st <- lavaan::standardizedSolution(cfa_result())
        ch <- st[st$op == "=~", c("lhs", "rhs", "est.std", "se", "ci.lower", "ci.upper")]
        h2 <- ch$est.std^2
        tab <- data.frame(Facteur = ch$lhs, Item = ch$rhs, `Charge standardisée` = ch$est.std,
                          `Erreur type` = ch$se, `IC inf.` = ch$ci.lower, `IC sup.` = ch$ci.upper,
                          `R² de l'item` = h2, check.names = FALSE)
        datatable(tab, rownames = FALSE, options = list(dom = "t", pageLength = 200, ordering = FALSE),
                  caption = "Dans la CFA, chaque item ne charge que sur les facteurs prévus par la syntaxe; toutes les autres charges sont fixées à 0.") %>%
          formatRound(3:7, 3, dec.mark = ",")
      })

      output$cfa_phi <- renderUI({
        req(cfa_result())
        st <- lavaan::standardizedSolution(cfa_result())
        cv <- st[st$op == "~~" & st$lhs != st$rhs & st$lhs %in% lavaan::lavNames(cfa_result(), "lv"), ]
        if (nrow(cv) == 0) return(NULL)
        tags$p(class = "mt-3", tags$b("Corrélations entre facteurs (CFA) : "),
               paste0(paste0(cv$lhs, "–", cv$rhs, " = ", fr(cv$est.std, 3), collapse = "; "), "."))
      })

      # Afficher le résumé du modèle CFA
      output$cfa_summary <- renderPrint({
        req(cfa_result())
        summary(cfa_result(), fit.measures = TRUE, standardized = TRUE)
      })
      
      # Afficher les indices d'ajustement
      output$cfa_fit_indices <- renderPrint({
        req(cfa_result())
        cat("INDICES D'AJUSTEMENT DU MODÈLE CFA\n")
        cat("================================\n\n")
        fit_indices <- fitMeasures(cfa_result(), fit.measures = "all")
        print(fit_indices)
      })

      # Lecture brève, volontairement prudente, des résultats CFA. Les seuils
      # sont des repères; ils ne remplacent ni les charges, ni les résidus, ni
      # la plausibilité substantielle du modèle.
      output$cfa_interpretation <- renderUI({
        if (is.null(cfa_result())) {
          return(tags$div(
            class = "alert alert-secondary mb-0",
            "Les résultats apparaîtront ici après l'évaluation du modèle sur le nouvel échantillon."
          ))
        }

        fit <- tryCatch(lavaan::fitMeasures(cfa_result()), error = function(e) numeric())
        pick_measure <- function(candidates) {
          available <- candidates[candidates %in% names(fit)]
          values <- unname(fit[available])
          values <- values[is.finite(values)]
          if (length(values) == 0) NA_real_ else values[1]
        }
        cfi <- pick_measure(c("cfi.robust", "cfi.scaled", "cfi"))
        tli <- pick_measure(c("tli.robust", "tli.scaled", "tli"))
        rmsea <- pick_measure(c("rmsea.robust", "rmsea.scaled", "rmsea"))
        srmr <- pick_measure(c("srmr"))
        converged <- tryCatch(lavaan::lavInspect(cfa_result(), "converged"), error = function(e) FALSE)
        fmt <- function(x) if (is.finite(x)) fr(x, 3) else "non disponible"

        if (!isTRUE(converged)) {
          status_class <- "alert-danger"
          verdict <- "Le modèle n'a pas convergé : les indices ne doivent pas être interprétés."
        } else if (all(is.finite(c(cfi, rmsea, srmr))) && cfi >= .95 && rmsea <= .06 && srmr <= .08) {
          status_class <- "alert-success"
          verdict <- "L'ajustement est compatible avec une bonne reproduction des corrélations dans ce nouvel échantillon simulé."
        } else if (all(is.finite(c(cfi, rmsea, srmr))) && cfi >= .90 && rmsea <= .08 && srmr <= .10) {
          status_class <- "alert-warning"
          verdict <- "L'ajustement est utilisable, mais il reste des écarts à examiner avant de considérer la structure comme stable."
        } else {
          status_class <- "alert-danger"
          verdict <- "L'ajustement est faible : la structure proposée par l'AFE se transporte mal dans ce nouvel échantillon, ou le modèle doit être repensé."
        }

        couleurs <- c("alert-success" = "#2e7d32", "alert-warning" = "#b7791f", "alert-danger" = "#c0392b")
        fonds <- c("alert-success" = "#eef7ee", "alert-warning" = "#fff8e6", "alert-danger" = "#fdecea")
        tags$div(
          class = "p-3 mb-0",
          style = paste0("background:", fonds[[status_class]], "; border-left: 5px solid ",
                         couleurs[[status_class]], "; border-radius: 4px;"),
          tags$b("Lecture d'ensemble. "), verdict,
          tags$ul(
            tags$li(paste0("CFI = ", fmt(cfi), " ; TLI = ", fmt(tli), ". Plus ils sont près de 1, mieux c'est.")),
            tags$li(paste0("RMSEA = ", fmt(rmsea), " ; SRMR = ", fmt(srmr), ". Plus ils sont près de 0, mieux c'est."))
          ),
          tags$p(class = "mb-0",
                 "Même un bon ajustement ne prouve pas le modèle. Relisez les charges standardisées, le sens des facteurs et les éventuels indices de modification avant toute révision.")
        )
      })
      
      # Diagramme du modèle CFA
      output$cfa_diagram <- renderPlot({
        req(cfa_result())
        
        withProgress(message = "Génération du diagramme du modèle...", {
          # Générer un diagramme du modèle avec semPaths
          tryCatch({
            if (!requireNamespace("semPlot", quietly = TRUE)) {
              plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
              text(0, 0, "Package 'semPlot' requis pour le diagramme", cex = 1.2)
            } else {
              semPlot::semPaths(cfa_result(),
                                intercepts = FALSE,
                                whatLabels = "std",
                                rotation = 2,
                                sizeMan = 6,
                                sizeMan2 = 3,
                                sizeLat = 8,
                                style = "lisrel",
                                curvature = 1.8,
                                edge.label.cex = .8,
                                label.cex = .8,
                                residScale = 12)
            }
          }, error = function(e) {
            plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
            text(0, 0, paste("Erreur lors de la génération du diagramme:", e$message), cex = 1.2)
          })
        })
      })
      
      # Afficher les indices de modification
      output$modification_indices <- renderDT({
        req(cfa_result())
        
        # Obtenir les indices de modification
        mod_indices <- modificationIndices(cfa_result(), sort = TRUE)
        
        # Si aucun index de modification, afficher un message
        if (nrow(mod_indices) == 0) {
          return(datatable(data.frame(Message = "Aucun indice de modification disponible"),
                           options = list(dom = 't'),
                           rownames = FALSE))
        }
        
        # Filtrer et trier les indices
        mod_indices <- mod_indices[mod_indices$mi > 3.84, ]  # Seuil chi² à df=1, p<0.05
        mod_indices <- mod_indices[order(-mod_indices$mi), ]
        
        # Formater pour l'affichage
        mod_indices_formatted <- mod_indices[, c("lhs", "op", "rhs", "mi", "epc")]
        colnames(mod_indices_formatted) <- c("De", "Opération", "Vers", "IM", "CPC")
        
        # Créer le tableau
        datatable(
          mod_indices_formatted,
          options = list(pageLength = 10, scrollX = TRUE),
          rownames = FALSE
        ) %>%
          formatRound(columns = c("IM", "CPC"), digits = 2) %>%
          formatStyle(
            "IM",
            background = styleInterval(c(3.84, 10), c("white", "#ffff9c", "#ffcccc")),
            fontWeight = styleInterval(10, c("normal", "bold"))
          )
      })
      
      # Une AFC n'est pas proposée automatiquement après une AFE sur les mêmes
      # données. La validation confirmatoire requiert idéalement un nouvel
      # échantillon ou un partage apprentissage/validation.
      
      # Assigne un item à son facteur principal seulement. Un item dont la
      # seconde charge est trop proche de la première ne sert pas au calcul de
      # fiabilité : l'inclure dans deux échelles gonflerait artificiellement
      # l'impression de cohérence interne.
      factor_items_for_reliability <- function(loadings, factor_idx,
                                               threshold = .30, cross_gap = .15) {
        abs_loadings <- abs(loadings)
        primary <- apply(abs_loadings, 1, which.max)
        primary_value <- apply(abs_loadings, 1, max)
        secondary_value <- apply(abs_loadings, 1, function(x) {
          if (length(x) < 2) return(0)
          sort(x, decreasing = TRUE)[2]
        })
        rownames(loadings)[primary == factor_idx &
                            primary_value >= threshold &
                            (primary_value - secondary_value >= cross_gap)]
      }

      # Calculer les indices de fiabilité pour chaque facteur
      calculate_reliability <- function(data, loadings, threshold = 0.3) {
        req(data, loadings)
        
        # Nombre de facteurs
        n_factors <- ncol(loadings)
        
        # Initialiser le tableau de résultats
        reliability_results <- data.frame(
          Facteur = paste0("F", 1:n_factors),
          Alpha = numeric(n_factors),
          Omega = numeric(n_factors),
          N = integer(n_factors),
          Variables = character(n_factors)
        )
        
        # Pour chaque facteur
        for (f in 1:n_factors) {
          # Retenir seulement les indicateurs non ambigus du facteur principal.
          items <- factor_items_for_reliability(loadings, f, threshold)
          
          # Si au moins 2 variables sont associées au facteur
          if (length(items) >= 2) {
            # Extraire les données pour ces variables
            factor_data <- data[, items, drop = FALSE]
            factor_data <- factor_data[complete.cases(factor_data), , drop = FALSE]
            reliability_results$N[f] <- nrow(factor_data)
            
            # Calculer alpha de Cronbach
            alpha_result <- tryCatch(
              suppressWarnings(psych::alpha(factor_data, check.keys = FALSE, warnings = FALSE)$total$raw_alpha),
              error = function(e) NA)
            
            # Calculer oméga de McDonald
            omega_result <- omega_sous_echelle(factor_data)$omega
            
            # Stocker les résultats
            reliability_results$Alpha[f] <- alpha_result
            reliability_results$Omega[f] <- omega_result
            reliability_results$Variables[f] <- paste(items, collapse = ", ")
          } else {
            # Si moins de 2 variables, la fiabilité ne peut pas être calculée
            reliability_results$Alpha[f] <- NA
            reliability_results$Omega[f] <- NA
            reliability_results$N[f] <- 0
            reliability_results$Variables[f] <- ifelse(length(items) == 0, 
                                                       "Aucune variable (seuil trop élevé)",
                                                       paste(items, collapse = ", "))
          }
        }
        
        return(reliability_results)
      }
      
      # Mettre à jour le sélecteur de facteurs lorsque l'EFA est terminée
      observe({
        req(efa_result())
        efa <- efa_result()
        
        # Mettre à jour le sélecteur de facteurs avec les facteurs disponibles
        updateSelectInput(session, "factor_for_analysis",
                          choices = setNames(1:efa$nfactors, 
                                             paste0("Facteur ", 1:efa$nfactors)))
      })
      
      # Tableau des indices de fiabilité
      output$reliability_table <- renderDT({
        req(efa_result(), filtered_data())
        
        efa <- efa_result()
        data <- filtered_data()
        loadings <- efa$loadings
        threshold <- input$loadings_threshold
        
        # Calculer les indices de fiabilité
        reliability_data <- calculate_reliability(data, loadings, threshold)
        
        # Formatage pour l'affichage selon la méthode sélectionnée
        if (input$reliability_method == "alpha") {
          display_data <- reliability_data[, c("Facteur", "Alpha", "N", "Variables")]
          colnames(display_data) <- c("Facteur", "Alpha de Cronbach", "N", "Variables incluses")
          
          # Créer le tableau avec coloration selon la valeur - ALPHA SEULEMENT
          datatable(
            display_data,
            options = list(dom = "t", pageLength = 100, ordering = FALSE, scrollX = TRUE),
            rownames = FALSE
          ) %>% formatRound(columns = "Alpha de Cronbach", digits = 3, dec.mark = ",")
          
        } else if (input$reliability_method == "omega") {
          display_data <- reliability_data[, c("Facteur", "Omega", "N", "Variables")]
          colnames(display_data) <- c("Facteur", "Oméga de McDonald", "N", "Variables incluses")
          
          # Créer le tableau avec coloration selon la valeur - OMEGA SEULEMENT
          datatable(
            display_data,
            options = list(dom = "t", pageLength = 100, ordering = FALSE, scrollX = TRUE),
            rownames = FALSE
          ) %>% formatRound(columns = "Oméga de McDonald", digits = 3, dec.mark = ",")
          
        } else { # both
          display_data <- reliability_data[, c("Facteur", "Alpha", "Omega", "N", "Variables")]
          colnames(display_data) <- c("Facteur", "Alpha de Cronbach", "Oméga de McDonald", "N", "Variables incluses")
          
          # Créer le tableau avec coloration selon la valeur - LES DEUX
          datatable(
            display_data,
            options = list(dom = "t", pageLength = 100, ordering = FALSE, scrollX = TRUE),
            rownames = FALSE
          ) %>% formatRound(columns = c("Alpha de Cronbach", "Oméga de McDonald"), digits = 3, dec.mark = ",")
        }
      })
      
      # Analyse détaillée de fiabilité pour un facteur spécifique
      output$detailed_reliability <- renderUI({
        req(efa_result(), filtered_data(), input$factor_for_analysis)
        efa <- efa_result()
        f <- as.numeric(input$factor_for_analysis)
        items <- factor_items_for_reliability(efa$loadings, f, input$loadings_threshold)
        if (length(items) < 2) {
          return(tags$div(class = "alert alert-warning",
                          "Fiabilité non calculable : il faut au moins deux items nets (charge principale au-dessus du seuil ",
                          fr(input$loadings_threshold, 2), " et seconde charge nettement plus faible) sur ce facteur."))
        }
        d <- na.omit(filtered_data()[, items, drop = FALSE])
        a <- tryCatch(suppressWarnings(psych::alpha(d, check.keys = FALSE, warnings = FALSE)), error = function(e) NULL)
        o <- omega_sous_echelle(d)
        tab <- data.frame(Item = items, `Charge (AFE)` = efa$loadings[items, f],
                          `λ (modèle à 1 facteur)` = o$lambda[items],
                          `Corrélation item-reste` = if (!is.null(a)) a$item.stats[items, "r.drop"] else NA,
                          `Alpha si l'item est retiré` = if (!is.null(a)) a$alpha.drop[items, "raw_alpha"] else NA,
                          check.names = FALSE)
        tagList(
          tags$p(tags$b(paste0("Facteur ", f, " : ")), length(items), " items, ", nrow(d), " cas complets."),
          tags$ul(
            tags$li(HTML(paste0("Alpha de Cronbach = <b>", fr(if (!is.null(a)) a$total$raw_alpha else NA, 3),
                                "</b> : suppose des items de même charge (modèle tau-équivalent)."))),
            tags$li(HTML(paste0("Oméga total = <b>", fr(o$omega, 3),
                                "</b> = (Σλ)² / [(Σλ)² + Σ(1 − λ²)] : tient compte de charges inégales (modèle congénérique).")))
          ),
          renderDT(datatable(tab, rownames = FALSE, options = list(dom = "t", ordering = FALSE)) %>%
                     formatRound(2:5, 3, dec.mark = ",")),
          tags$p(class = "small text-muted mt-2",
                 "Un « alpha si l'item est retiré » plus élevé ne justifie pas à lui seul de retirer l'item : vérifiez d'abord ce que son contenu apporte à la mesure.")
        )
      })

      # Réactive pour stocker les résultats du bootstrap
      bootstrap_results <- reactiveVal(NULL)
      
      # Fonction pour exécuter l'analyse bootstrap
      observeEvent(input$run_bootstrap, {
        req(efa_result())
        efa <- efa_result()
        selected_vars <- efa$variables
        nfactors <- efa$nfactors
        rotation <- efa$rotation
        method <- efa$method
        n_samples <- input$bootstrap_samples
        ci_level <- input$bootstrap_ci / 100
        
        # Vérifier que les données sont suffisantes
        if (length(selected_vars) < nfactors + 1) {
          showNotification(
            "Erreur: Trop peu de variables pour le nombre de facteurs spécifié",
            type = "error",
            duration = 5
          )
          return()
        }
        
        withProgress(message = "Exécution du bootstrap...", {
          tryCatch({
            all_loadings <- list()
            # Les items et les cas sont ceux de la solution affichée (y compris
            # après un recalcul sans certains items).
            data_complete <- efa$data_complete
            n_obs <- nrow(data_complete)
            incProgress(0.1, detail = "Initialisation...")
            for (i in seq_len(n_samples)) {
              if (i %% 10 == 0) {
                incProgress(0.8 * (10 / n_samples),
                            detail = paste0("Échantillon ", i, " sur ", n_samples))
              }
              boot_data <- data_complete[sample.int(n_obs, n_obs, replace = TRUE), , drop = FALSE]
              res <- tryCatch(ajuster_afe(boot_data, nfactors, rotation, method, maxit = 500),
                              error = function(e) NULL)
              if (!is.null(res) && res$converged) {
                all_loadings[[i]] <- align_bootstrap_loadings(efa$loadings,
                                                              res$loadings[rownames(efa$loadings), , drop = FALSE])
              }
            }

            # Filtrer les échantillons qui ont échoué
            all_loadings <- all_loadings[!sapply(all_loadings, is.null)]
            n_successful <- length(all_loadings)

            if (n_successful == 0) {
              stop("Aucune réplication bootstrap n’a convergé; la stabilité ne peut pas être estimée.")
            }
            
            if (n_successful < n_samples * 0.5) {
              showNotification(
                paste0("Attention: Seulement ", n_successful, " échantillons sur ", 
                       n_samples, " ont convergé. Les résultats peuvent être instables."),
                type = "warning",
                duration = 10
              )
            }
            
            incProgress(0.9, detail = "Calcul des statistiques...")
            
            # Initialiser des tableaux pour les statistiques
            var_names <- rownames(efa$loadings)
            factor_names <- colnames(efa$loadings)
            n_vars <- length(var_names)
            
            # Tableaux pour stocker les résultats
            mean_loadings <- matrix(0, nrow = n_vars, ncol = nfactors)
            sd_loadings <- matrix(0, nrow = n_vars, ncol = nfactors)
            ci_lower <- matrix(0, nrow = n_vars, ncol = nfactors)
            ci_upper <- matrix(0, nrow = n_vars, ncol = nfactors)
            
            rownames(mean_loadings) <- rownames(sd_loadings) <- 
              rownames(ci_lower) <- rownames(ci_upper) <- var_names
            colnames(mean_loadings) <- colnames(sd_loadings) <- 
              colnames(ci_lower) <- colnames(ci_upper) <- factor_names
            
            # Pour chaque variable et chaque facteur
            for (v in 1:n_vars) {
              for (f in 1:nfactors) {
                # Extraire toutes les charges factorielles pour cette variable/facteur
                loadings_vf <- sapply(all_loadings, function(l) l[v, f])
                
                # Calculer la moyenne et l'écart-type
                mean_loadings[v, f] <- mean(loadings_vf, na.rm = TRUE)
                sd_loadings[v, f] <- sd(loadings_vf, na.rm = TRUE)
                
                # Calculer les intervalles de confiance
                alpha <- (1 - ci_level) / 2
                ci_lower[v, f] <- quantile(loadings_vf, alpha, na.rm = TRUE)
                ci_upper[v, f] <- quantile(loadings_vf, 1 - alpha, na.rm = TRUE)
              }
            }
            
            # Stocker tous les résultats
            bootstrap_results(list(
              mean_loadings = mean_loadings,
              sd_loadings = sd_loadings,
              ci_lower = ci_lower,
              ci_upper = ci_upper,
              raw_loadings = all_loadings,
              n_samples = n_successful,
              ci_level = ci_level
            ))
            
            # Mettre à jour les sélecteurs de facteurs et variables
            updateSelectInput(session, "bootstrap_display_factor", 
                              choices = setNames(1:nfactors, paste0("Facteur ", 1:nfactors)),
                              selected = 1)
            
            updateSelectInput(session, "bootstrap_stability_var",
                              choices = var_names,
                              selected = var_names[1])
            
            showNotification(paste0("Analyse bootstrap terminée (", n_successful, " échantillons)"), 
                             type = "message")
            
          }, error = function(e) {
            showNotification(
              paste("Erreur lors de l'analyse bootstrap:", e$message),
              type = "error",
              duration = 10
            )
          })
        })
      })
      
      # Tableau des intervalles de confiance des charges
      output$bootstrap_loadings_table <- renderDT({
        req(bootstrap_results(), input$bootstrap_display_factor, efa_result())
        res <- bootstrap_results()
        f <- as.numeric(input$bootstrap_display_factor)
        seuil <- input$loadings_threshold
        lo <- res$ci_lower[, f]; hi <- res$ci_upper[, f]
        lecture <- ifelse(lo > seuil | hi < -seuil, "au-delà du seuil",
                   ifelse(lo > -seuil & hi < seuil, "sous le seuil", "traverse le seuil"))
        tab <- data.frame(
          Item = rownames(res$mean_loadings),
          `Charge observée` = efa_result()$loadings[, f],
          `Moyenne bootstrap` = res$mean_loadings[, f],
          `Écart-type` = res$sd_loadings[, f],
          `Borne inf.` = lo, `Borne sup.` = hi, Largeur = hi - lo,
          `Par rapport au seuil` = lecture,
          check.names = FALSE)
        datatable(tab, rownames = FALSE,
                  options = list(dom = "t", pageLength = 200, ordering = FALSE, scrollX = TRUE),
                  caption = paste0(res$n_samples, " réplications alignées; intervalles percentiles à ",
                                   round(100 * res$ci_level), " %. « Traverse le seuil » (", fr(seuil, 2),
                                   ") : l'échantillon ne permet pas de dire si l'item atteint ce degré de force.")) %>%
          formatRound(2:7, 3, dec.mark = ",") %>%
          formatStyle("Charge observée", fontWeight = "bold") %>%
          formatStyle("Par rapport au seuil", backgroundColor = styleEqual(
            c("au-delà du seuil", "traverse le seuil", "sous le seuil"), c("#e3f1da", "#fff3cd", "white")))
      })

      # Graphique de la stabilité des charges
      output$bootstrap_stability_plot <- renderPlotly({
        req(bootstrap_results(), input$bootstrap_stability_var)
        
        results <- bootstrap_results()
        selected_var <- input$bootstrap_stability_var
        var_idx <- match(selected_var, rownames(results$mean_loadings))
        
        # S'assurer que la variable existe
        if (is.na(var_idx)) {
          return(plot_ly() %>%
            add_annotations(
              text = "Variable non trouvée dans les résultats",
              showarrow = FALSE,
              font = list(size = 14)
            ))
        }
        
        # Créer un dataframe avec les distributions des charges pour cette variable
        factor_names <- colnames(results$mean_loadings)
        n_factors <- length(factor_names)
        
        # Extraire toutes les charges bootstrappées pour cette variable
        all_samples <- results$raw_loadings
        n_samples <- length(all_samples)
        
        # Créer un dataframe pour stocker toutes les charges bootstrappées
        plot_data <- data.frame(
          Facteur = character(),
          Charge = numeric(),
          stringsAsFactors = FALSE
        )
        
        # Remplir le dataframe
        for (f in 1:n_factors) {
          factor_loadings <- sapply(all_samples, function(l) l[var_idx, f])
          plot_data <- rbind(
            plot_data,
            data.frame(
              Facteur = factor_names[f],
              Charge = factor_loadings
            )
          )
        }
        
        # Ajouter les valeurs originales pour référence
        original_loadings <- efa_result()$loadings[var_idx, ]
        
        # Créer le graphique violin plot
        p <- plot_ly() %>%
          add_trace(
            data = plot_data,
            x = ~Facteur,
            y = ~Charge,
            type = "violin",
            box = list(visible = TRUE),
            meanline = list(visible = TRUE),
            color = ~Facteur,
            points = "outliers"
          ) %>%
          add_trace(
            x = factor_names,
            y = original_loadings,
            type = "scatter",
            mode = "markers",
            marker = list(size = 10, symbol = "star", color = "black"),
            name = "Valeur originale"
          ) %>%
          layout(
            title = paste("Distribution bootstrap des charges de", selected_var, "(étoile : valeur observée)"),
            separators = ", ",
            xaxis = list(title = "Facteur"),
            yaxis = list(title = "Charge factorielle", range = c(-1, 1)),
            showlegend = FALSE
          )
        
        return(p)
      })
      
      # Mettre à jour les contrôles du bootstrap quand l'analyse factorielle est exécutée
      observe({
        req(efa_result())
        
        # Obtenir les noms des variables et facteurs
        var_names <- rownames(efa_result()$loadings)
        n_factors <- ncol(efa_result()$loadings)
        
        # Mettre à jour les sélecteurs
        updateSelectInput(
          session,
          "bootstrap_stability_var",
          choices = var_names,
          selected = if(length(var_names) > 0) var_names[1] else NULL
        )
        
        updateSelectInput(
          session,
          "bootstrap_display_factor",
          choices = setNames(1:n_factors, paste0("Facteur ", 1:n_factors)),
          selected = if(n_factors > 0) 1 else NULL
        )
      })
      
      # Mettre à jour les sélecteurs après le bootstrap également
      observeEvent(bootstrap_results(), {
        req(bootstrap_results(), efa_result())
        
        results <- bootstrap_results()
        
        # Mettre à jour les sélecteurs de facteurs et variables
        var_names <- rownames(efa_result()$loadings)
        n_factors <- ncol(efa_result()$loadings)
        
        updateSelectInput(
          session,
          "bootstrap_display_factor", 
          choices = setNames(1:n_factors, paste0("Facteur ", 1:n_factors)),
          selected = 1
        )
        
        updateSelectInput(
          session,
          "bootstrap_stability_var",
          choices = var_names,
          selected = var_names[1]
        )
      })
      
      observeEvent(suppressWarnings(event_data("plotly_click", source = "biplot", priority = "event")), {
        req(efa_result())
        click_data <- suppressWarnings(event_data("plotly_click", source = "biplot"))
        
        # Si on a des données de clic
        if (!is.null(click_data)) {
          # Récupérer la variable cliquée
          clicked_var <- click_data$customdata
          
          # Récupérer la liste actuelle des variables désactivées
          current_disabled <- input$disabled_variables
          if (is.null(current_disabled)) {
            current_disabled <- c()
          }
          
          # Basculer l'état de la variable (activer/désactiver)
          if (clicked_var %in% current_disabled) {
            # Si déjà désactivée, réactiver
            new_disabled <- setdiff(current_disabled, clicked_var)
          } else {
            # Si active, désactiver
            new_disabled <- c(current_disabled, clicked_var)
          }
          
          # Mettre à jour la liste des variables désactivées
          updateSelectizeInput(session, "disabled_variables", selected = new_disabled)
          
          # Si le recalcul automatique est activé, relancer l'analyse factorielle
          if (input$auto_recalculate) {
            # Déclencher le recalcul
            recalculate_efa_with_active_vars()
          }
        }
      })
      
      # Handler pour le bouton de recalcul manuel
      observeEvent(input$manual_recalculate, {
        recalculate_efa_with_active_vars()
      })
      
      # Ajouter cette fonction dans la partie serveur
      recalculate_efa_with_active_vars <- function() {
        req(efa_result(), filtered_data(), input$efa_variables)
        bootstrap_results(NULL)
        cfa_result(NULL)
        active_vars <- setdiff(input$efa_variables, input$disabled_variables)
        if (length(active_vars) <= input$efa_nfactors) {
          showNotification("Il doit rester plus d'items actifs que de facteurs.", type = "error")
          return()
        }
        withProgress(message = "Recalcul sans les items désactivés...", {
          tryCatch({
            data_subset <- filtered_data()[, active_vars, drop = FALSE]
            efa_result(ajuster_afe(na.omit(data_subset), input$efa_nfactors,
                                   input$dynamic_rotation, input$efa_method,
                                   n_total = nrow(data_subset)))
            showNotification(paste0("Solution recalculée avec ", length(active_vars), " items."),
                             type = "message")
          }, error = function(e) {
            showNotification(paste("Erreur lors du recalcul :", e$message), type = "error", duration = 10)
          })
        })
      }

      # Observer qui met à jour les sélecteurs de facteurs lorsque l'analyse factorielle est exécutée
      observe({
        req(efa_result())
        efa <- efa_result()
        n_factors <- efa$nfactors
        
        # Mettre à jour les sélecteurs de facteurs
        updateSelectInput(session, "plot_factor_x", 
                          choices = setNames(1:n_factors, paste0("Facteur ", 1:n_factors)),
                          selected = 1)
        
        updateSelectInput(session, "plot_factor_y", 
                          choices = setNames(1:n_factors, paste0("Facteur ", 1:n_factors)),
                          selected = min(2, n_factors))
        
        # Mettre à jour les variables désactivées
        updateSelectizeInput(session, "disabled_variables", 
                             choices = rownames(efa$loadings), 
                             selected = character(0))
      })
      
      # Ajouter cette fonction pour générer une sortie HTML avec les informations sur l'analyse recalculée
      output$recalculated_info <- renderUI({
        req(efa_result())
        efa <- efa_result()
        
        # Créer un affichage HTML des informations clés
        html_content <- tags$div(
          tags$p(tags$b("Méthode d'extraction:"), efa$method),
          tags$p(tags$b("Méthode de rotation:"), efa$rotation),
          tags$p(tags$b("Nombre de facteurs:"), efa$nfactors),
          tags$p(tags$b("Variables actives:"), length(rownames(efa$loadings))),
          tags$hr(),
          tags$p(tags$b("Variance expliquée (%):"), style = "margin-bottom: 5px;")
        )
        
        # Ajouter un tableau pour la variance expliquée
        variance_table <- tags$table(
          style = "width: 100%; border-collapse: collapse;",
          tags$thead(
            tags$tr(
              tags$th(style = "border: 1px solid #ddd; padding: 4px; text-align: left;", "Facteur"),
              tags$th(style = "border: 1px solid #ddd; padding: 4px; text-align: right;", "% Var"),
              tags$th(style = "border: 1px solid #ddd; padding: 4px; text-align: right;", "% Cum")
            )
          ),
          tags$tbody(
            lapply(1:efa$nfactors, function(i) {
              tags$tr(
                tags$td(style = "border: 1px solid #ddd; padding: 4px;", paste0("F", i)),
                tags$td(style = "border: 1px solid #ddd; padding: 4px; text-align: right;", 
                        sprintf("%.1f", efa$Vaccounted$Proportion[i] * 100)),
                tags$td(style = "border: 1px solid #ddd; padding: 4px; text-align: right;", 
                        sprintf("%.1f", efa$Vaccounted$Cumulative[i] * 100))
              )
            })
          )
        )
        
        # Ajouter un bouton pour afficher plus de détails
        detail_button <- tags$button(
          id = "show_detail_modal",
          class = "btn btn-sm btn-secondary mt-2",
          style = "width: 100%",
          onclick = "Shiny.setInputValue('show_detailed_results', Math.random(), {priority: 'event'});",
          "Afficher les détails complets"
        )
        
        return(tags$div(html_content, variance_table, detail_button))
      })
      
      # Ajouter ce gestionnaire d'événements pour afficher une fenêtre modale avec les résultats complets
      observeEvent(input$show_detailed_results, {
        req(efa_result())
        
        showModal(modalDialog(
          title = "Résultats détaillés de l'analyse factorielle",
          size = "l",
          easyClose = TRUE,
          
          tabsetPanel(
            tabPanel("Charges factorielles", 
                     DTOutput("modal_loadings_table")),
            tabPanel("Communalités", 
                     DTOutput("modal_communalities_table")),
            tabPanel("Indices d'ajustement", 
                     verbatimTextOutput("modal_fit_indices"))
          ),
          
          footer = modalButton("Fermer")
        ))
      })
      
      # Tableaux pour la fenêtre modale
      output$modal_loadings_table <- renderDT({
        req(efa_result())
        loadings <- efa_result()$loadings
        
        # Formater pour l'affichage
        loadings_formatted <- round(loadings, 3)
        
        # Créer un tableau avec coloration
        datatable(
          loadings_formatted,
          options = list(pageLength = 15, scrollX = TRUE),
          rownames = TRUE
        ) %>%
          formatStyle(
            columns = names(loadings_formatted),
            backgroundColor = styleInterval(
              c(-0.3, 0.3),
              c("#e15759", "white", "#a8d08d")
            ),
            color = styleInterval(
              c(-0.3, 0.3),
              c("#ffffff", "black", "#000000")
            )
          )
      })
      
      output$modal_communalities_table <- renderDT({
        req(efa_result())
        
        comm <- efa_result()$communality
        uniq <- efa_result()$uniquenesses
        
        # Préparer les données pour le tableau
        comm_data <- data.frame(
          Variable = names(comm),
          Communalité = round(comm, 3),
          Unicité = round(uniq, 3)
        )
        
        # Créer le tableau
        datatable(
          comm_data,
          options = list(pageLength = 15, scrollX = TRUE),
          rownames = FALSE
        ) %>%
          formatStyle(
            "Communalité",
            background = styleColorBar(c(0, 1), "lightsteelblue"),
            backgroundSize = "100% 90%",
            backgroundRepeat = "no-repeat",
            backgroundPosition = "center"
          )
      })
      
      output$modal_fit_indices <- renderPrint({
        req(efa_result())
        
        fit_indices <- efa_result()$fit_indices
        
        cat("INDICES D'AJUSTEMENT DU MODÈLE RECALCULÉ\n")
        cat("=====================================\n\n")
        
        # Sélectionner et afficher les indices les plus importants
        important_indices <- c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr", "aic", "bic")
        
        if (!is.null(fit_indices)) {
          for (idx in important_indices) {
            if (idx %in% names(fit_indices)) {
              cat(sprintf("%-10s: %8.4f\n", idx, fit_indices[idx]))
            }
          }
        } else {
          cat("Aucun indice d'ajustement disponible.")
        }
      })
      
      onSessionEnded(function() {
        # Nettoyer la mémoire lorsque l'utilisateur ferme l'application
        cleanup_reactive_objects()
        message("Session terminée - Mémoire libérée")
      })
      
      # Rendu des interprétations automatiques
      # Dans la partie rendu des interprétations automatiques, modifiez la section qui crée les cartes :
      
      output$auto_interpretation <- renderUI({
        req(efa_result())
        
        # Obtenir les charges factorielles et le seuil actuel
        loadings <- efa_result()$loadings
        threshold <- input$loadings_threshold
        
        # Générer les interprétations
        interpretations <- generate_factor_interpretation(loadings, threshold)
        
        # Créer des éléments d'interface pour chaque facteur
        factor_cards <- lapply(interpretations, function(interp) {
          # Définir la couleur en fonction de la qualité du facteur
          quality_colors <- list(
            "excellent" = "#28a745",  # vert
            "bon" = "#5cb85c",        # vert clair
            "modéré" = "#ffc107",     # jaune
            "faible" = "#dc3545",     # rouge
            "insuffisant" = "#6c757d" # gris
          )
          
          # Obtenir la couleur appropriée
          quality_color <- quality_colors[[interp$quality]]
          
          div(
            class = "card mb-3",
            div(
              class = "card-header",
              style = paste0("background-color: ", quality_color, "; color: white;"),
              h4(interp$name, style = "margin: 0;")
            ),
            div(
              class = "card-body",
              div(
                class = "alert", 
                style = paste0("background-color: ", quality_color, "20; color: black;"), # Ajout de "color: black"
                p(strong("Netteté du facteur : "), interp$quality),
                p(interp$clarity)
              ),
              h5("Variables associées"),
              p(interp$pos_summary),
              if (length(interp$negative_vars) > 0) p(interp$neg_summary) else NULL,
              hr(),
              if (interp$n_vars >= 3) div(
                class = "alert alert-info",
                style = "color: black;", # Ajout de cette ligne pour garantir le texte en noir
                h5("Suggestion d'interprétation"),
                p(interp$suggestion),
                if (interp$n_vars >= 5) {
                  top_vars <- interp$positive_vars[1:min(5, length(interp$positive_vars))]
                  p("Variables clés à examiner : ", strong(paste(top_vars, collapse=", ")))
                }
              ) else div(
                class = "alert alert-warning",
                style = "color: black;", # Ajout de cette ligne pour garantir le texte en noir
                p("Nombre insuffisant de variables avec des charges significatives pour une interprétation fiable.")
              )
            )
          )
        })
        
        # Renvoyer tous les éléments
        div(factor_cards)
      })
      
      # Description des jeux de données classiques
      dataset_descriptions <- list(
        hs1939 = list(title = "Holzinger et Swineford (1939)",
          desc = "Neuf tests d'aptitudes passés par 301 élèves de deux écoles. Trois facteurs attendus : visuel (x1-x3), verbal (x4-x6), vitesse (x7-x9).",
          variables = "x1 à x9", observations = "301", source = "lavaan"),
        sdt = list(title = "Besoins psychologiques fondamentaux (théorie de l'autodétermination)",
          desc = "Douze items, quatre par besoin : appartenance (R1-R4), compétence (C1-C4), autonomie (A1-A4). Trois facteurs corrélés attendus.",
          variables = "R1-R4, C1-C4, A1-A4", observations = "366", source = "EFA.dimensions (data_SDT)"),
        grips = list(title = "General Risk Propensity Scale (GRiPS)",
          desc = "Huit items sur l'attrait général pour le risque (0 à 5). Échelle conçue comme unidimensionnelle : un bon exemple de solution à un facteur.",
          variables = "fun, friends, enjoy, hurt, part, commonly, chances, attracted", observations = "810", source = "EFAtools (GRiPS_raw)"),
        rse = list(title = "Échelle d'estime de soi de Rosenberg (items non recodés)",
          desc = "Dix items (0 à 4). Les items Q3, Q5, Q8, Q9 et Q10 sont formulés négativement et n'ont PAS été inversés : leurs corrélations avec les autres sont négatives. À utiliser avec « Préparation des items ».",
          variables = "Q1 à Q10", observations = "300", source = "EFA.dimensions (data_RSE_not_recoded)"),
        bfi = list(title = "Big Five Inventory (25 items)",
          desc = "Cinq items par trait : agréabilité (A), conscience (C), extraversion (E), névrosisme (N), ouverture (O). Les items A1, C4, C5, E1, E2, O2 et O5 sont formulés à l'envers; le fichier contient des valeurs manquantes.",
          variables = "A1-A5, C1-C5, E1-E5, N1-N5, O1-O5", observations = "2 800", source = "psych (bfi)"),
        field = list(title = "Questionnaire d'anxiété face aux statistiques (SAQ)",
          desc = "Vingt-trois items de l'exemple d'Andy Field (anxiété face aux statistiques, aux ordinateurs, aux mathématiques, peur de l'évaluation). L'item Q03 est formulé à l'envers.",
          variables = "Q01 à Q23", observations = "2 571", source = "EFA.dimensions (data_Field)"),
        dospert = list(title = "Domain-Specific Risk-Taking (DOSPERT)",
          desc = "Trente items de probabilité de prise de risque (1 à 7), six par domaine : éthique, financier, santé et sécurité, récréatif, social.",
          variables = "ethR, finR, heaR, recR, socR (1 à 6)", observations = "3 123", source = "EFAtools (DOSPERT_raw)"),
        neo = list(title = "NEO-PI-R : trente facettes",
          desc = "Scores de six facettes pour chacun des cinq grands traits (n, e, o, a, c). Une structure à cinq facteurs est attendue, avec plusieurs facettes complexes.",
          variables = "n1-n6, e1-e6, o1-o6, a1-a6, c1-c6", observations = "1 000", source = "EFA.dimensions (data_NEOPIR)")
      )

      # Afficher les informations sur le jeu de données sélectionné
      output$dataset_info <- renderUI({
        req(input$classic_dataset)
        
        if (input$classic_dataset == "") {
          return(NULL)
        }
        
        info <- dataset_descriptions[[input$classic_dataset]]
        
        tagList(
          h5(info$title, style = "margin-top: 0;"),
          p(info$desc),
          p(tags$b("Variables :"), info$variables),
          p(tags$b("Observations :"), info$observations),
          p(tags$b("Source :"), info$source)
        )
      })
      
      # Chargement d'un jeu de données classique
      observeEvent(input$load_classic, {
        req(input$classic_dataset, input$classic_dataset != "")
        
        cleanup_reactive_objects()
        id <- showNotification("Chargement du jeu de données...", type = "message", duration = NULL)
        
        tryCatch({
          # Charger le jeu de données sélectionné
          lire <- function(nom, paquet) {
            if (!requireNamespace(paquet, quietly = TRUE)) stop("le paquet ", paquet, " n'est pas installé.")
            e <- new.env()
            utils::data(list = nom, package = paquet, envir = e)
            as.data.frame(get(nom, envir = e))
          }
          data <- switch(input$classic_dataset,
            hs1939  = lavaan::HolzingerSwineford1939[, paste0("x", 1:9)],
            sdt     = lire("data_SDT", "EFA.dimensions"),
            grips   = lire("GRiPS_raw", "EFAtools"),
            rse     = lire("data_RSE_not_recoded", "EFA.dimensions"),
            bfi     = lire("bfi", "psych")[, 1:25],
            field   = lire("data_Field", "EFA.dimensions"),
            dospert = lire("DOSPERT_raw", "EFAtools"),
            neo     = lire("data_NEOPIR", "EFA.dimensions"),
            stop("jeu de données inconnu."))

          # Stocker les données
          simulation_spec(NULL)
          validation_data(NULL)
          cfa_result(NULL)
          hideTab(inputId = "main_tabs", target = "Validation CFA")
          original_data(data)
          raw_data(data)
          filtered_data(data)
          preparation_log(preparation_log()[0, ])
          
          # Mettre à jour les sélecteurs
          updateSelectInput(session, "filter_var", 
                            choices = c("Aucun" = "", names(data)))
          
          numeric_cols <- names(data)[sapply(data, is.numeric)]
          updateSelectInput(session, "efa_variables", 
                            choices = numeric_cols,
                            selected = numeric_cols)
          updateSelectizeInput(session, "reverse_items", choices = numeric_cols,
                               selected = character(0))
          
          # Mettre à jour les sélecteurs pour la visualisation
          updateSelectInput(session, "viz_variables", 
                            choices = numeric_cols,
                            selected = numeric_cols[1:min(5, length(numeric_cols))])
          
          # Afficher l'onglet d'analyse
          showTab(inputId = "main_tabs", target = "Analyse factorielle")
          
          removeNotification(id)
          showNotification(paste("Données classiques chargées :", input$classic_dataset), 
                           type = "message")
          
        }, error = function(e) {
          # Gérer les erreurs
          removeNotification(id)
          showNotification(
            paste("Erreur lors du chargement du jeu de données:", e$message),
            type = "error",
            duration = 10
          )
        })
      })

} 

shinyApp(ui, server)

